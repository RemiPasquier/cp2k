!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

!> \brief Raw CD pivot products for the nonperiodic large-cell Gamma tensor route.
MODULE gw_cd_lcg
   USE cp_fm_diag,                      ONLY: cp_fm_syevd
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE dbt_api,                         ONLY: &
        dbt_clear, dbt_copy, dbt_create, dbt_filter, dbt_get_info, dbt_get_num_blocks_total, &
        dbt_get_nze_total, dbt_put_block, dbt_reserve_blocks, dbt_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE machine,                         ONLY: m_walltime
   USE mp2_cd_producer,                 ONLY: cd_column_consumer_type,&
                                              cd_product_columns,&
                                              cd_select_products
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: qs_environment_type
#include "base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: cd_lcg_select, cd_lcg_metrics, cd_lcg_build_T, cd_lcg_get_T, cd_lcg_matrix_pair
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_cd_lcg'
   ! Per-invocation state; arrays retain their original allocation and ownership.
   TYPE, EXTENDS(cd_column_consumer_type) :: cd_metric_consumer_type
      TYPE(post_scf_bandstructure_type), POINTER :: bs_env => NULL()
      INTEGER :: nrow
      INTEGER :: ncol
      INTEGER, POINTER :: rows(:) => NULL()
      INTEGER, POINTER :: cols(:) => NULL()
      INTEGER, ALLOCATABLE :: local_row(:)
      INTEGER, ALLOCATABLE :: local_col(:)
      REAL(dp), CONTIGUOUS, POINTER :: local_data(:, :) => NULL()
      REAL(dp), ALLOCATABLE :: chk(:)
      REAL(dp), ALLOCATABLE :: prj(:)
   CONTAINS
      PROCEDURE :: consume => store_metric_column
   END TYPE cd_metric_consumer_type

   ! Per-invocation state; arrays retain their original allocation and ownership.
   TYPE, EXTENDS(cd_column_consumer_type) :: cd_T_consumer_type
      TYPE(post_scf_bandstructure_type), POINTER :: bs_env => NULL()
      INTEGER :: coord(3)
      INTEGER :: nproc
      INTEGER :: me
      INTEGER, ALLOCATABLE :: dist_i(:)
      INTEGER, ALLOCATABLE :: dist_j(:)
      INTEGER, ALLOCATABLE :: dist_k(:)
      INTEGER, ALLOCATABLE :: copy_i(:, :)
      INTEGER, ALLOCATABLE :: copy_j(:, :)
      INTEGER, ALLOCATABLE :: copy_k(:, :)
      INTEGER, ALLOCATABLE :: rank_map(:, :, :, :)
      LOGICAL, ALLOCATABLE :: reserved(:, :, :)
      INTEGER(int_8) :: n_sent
      REAL(dp), ALLOCATABLE :: chk(:)
      REAL(dp), ALLOCATABLE :: prj(:)
   CONTAINS
      PROCEDURE :: consume => store_T_column
   END TYPE cd_T_consumer_type

CONTAINS

!> \brief Stable owner grouping, retaining original selection order as explicit metadata.
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE cd_lcg_select(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      INTEGER                                            :: atom, handle, n, p, q
      INTEGER, ALLOCATABLE                               :: mu(:), nu(:), owner(:)
      REAL(dp)                                           :: elapsed, tick

      CALL timeset('cd_lcg_select', handle)
      CPASSERT(ALL(bs_env%periodic == 0))
      tick = m_walltime()
      CALL cd_select_products(qs_env, bs_env%para_env, bs_env%cd_lcg_eps, bs_env%unit_nr, mu, nu, owner)
      elapsed = m_walltime() - tick
      CALL bs_env%para_env%max(elapsed)
      n = SIZE(mu)
      IF (n == 0) CPABORT("CD_LCG: selector returned an empty auxiliary space")
      ALLOCATE (bs_env%cd_lcg_mu(n), bs_env%cd_lcg_nu(n), bs_env%cd_lcg_owner(n), bs_env%cd_lcg_order(n))
      q = 0
      DO atom = 1, bs_env%n_atom
         bs_env%sizes_RI(atom) = COUNT(owner == atom)
         DO p = 1, n
            IF (owner(p) /= atom) CYCLE
            q = q + 1
            bs_env%cd_lcg_mu(q) = mu(p)
            bs_env%cd_lcg_nu(q) = nu(p)
            bs_env%cd_lcg_owner(q) = atom
            bs_env%cd_lcg_order(q) = p
         END DO
      END DO
      CPASSERT(q == n)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(A,I0)') 'CD_LCG| naux ', n
         WRITE (bs_env%unit_nr, '(A,3(I0,1X))') 'CD_LCG| owner_counts min max zero ', &
            MINVAL(bs_env%sizes_RI), MAXVAL(bs_env%sizes_RI), COUNT(bs_env%sizes_RI == 0)
         WRITE (bs_env%unit_nr, '(A,F15.6)') 'CD_LCG| selection_seconds ', elapsed
      END IF
      ! The source snapshot does not include DBCSR's zero-block contract. Fail
      ! before creating either its matrices or DBT's flattened block layouts.
      IF (ANY(bs_env%sizes_RI == 0)) THEN
         CPABORT("CD_LCG: zero-pivot owner atom; zero-size DBCSR/DBT blocks are not validated")
      END IF
      CALL timestop(handle)
   END SUBROUTINE cd_lcg_select

!> \brief Build both unregularized product metrics before any factorization.
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE cd_lcg_metrics(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      INTEGER                                            :: handle, operator_index, q
      REAL(dp)                                           :: cutoff, elapsed, tick
      TYPE(cd_metric_consumer_type)                      :: consumer

      consumer%bs_env => bs_env

      CALL timeset('cd_lcg_metrics', handle)
      ALLOCATE (consumer%chk(SIZE(bs_env%cd_lcg_mu)), consumer%prj(SIZE(bs_env%cd_lcg_mu)))
      CALL cp_fm_create(bs_env%cd_lcg_M, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(bs_env%cd_lcg_V, bs_env%fm_RI_RI%matrix_struct)
      ! M and V share this matrix structure. Invert its local index maps once.
      CALL cp_fm_get_info(bs_env%cd_lcg_M, nrow_local=consumer%nrow, ncol_local=consumer%ncol, &
                          row_indices=consumer%rows, col_indices=consumer%cols)
      ALLOCATE (consumer%local_row(SIZE(consumer%chk)), consumer%local_col(SIZE(consumer%chk)), source=0)
      DO q = 1, consumer%nrow
         consumer%local_row(consumer%rows(q)) = q
      END DO
      DO q = 1, consumer%ncol
         consumer%local_col(consumer%cols(q)) = q
      END DO
      DO operator_index = 1, 2
         IF (operator_index == 1) THEN
            CALL cp_fm_get_info(bs_env%cd_lcg_M, local_data=consumer%local_data)
         ELSE
            CALL cp_fm_get_info(bs_env%cd_lcg_V, local_data=consumer%local_data)
         END IF
         tick = m_walltime()
         consumer%chk = 0.0_dp
         consumer%prj = 0.0_dp
         cutoff = 0.0_dp
         IF (operator_index == 1) cutoff = bs_env%ri_metric%cutoff_radius
         CALL cd_product_columns(qs_env, bs_env%para_env, bs_env%cd_lcg_mu, bs_env%cd_lcg_nu, &
                                 cutoff, 0.0_dp, .TRUE., consumer)
         CALL cd_lcg_reduce_diagnostics(bs_env, consumer%chk, consumer%prj)
         elapsed = m_walltime() - tick
         ! Partial-column sums change the floating-point reduction order, not the metric entries.
         IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A,I0,1X,ES24.16)') &
            'CD_LCG| raw_checksum operator ', operator_index, SUM(consumer%chk)
         IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A,I0,1X,ES24.16)') &
            'CD_LCG| raw_projection operator ', operator_index, cd_lcg_combine(consumer%prj)
         CALL bs_env%para_env%max(elapsed)
         IF (operator_index == 1) THEN
            IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A,F15.6)') 'CD_LCG| M_Gamma_seconds ', elapsed
            CALL cd_lcg_spectrum(bs_env%cd_lcg_M, bs_env, 'M_Gamma')
         ELSE
            IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A,F15.6)') 'CD_LCG| V_seconds ', elapsed
            CALL cd_lcg_spectrum(bs_env%cd_lcg_V, bs_env, 'V')
         END IF
      END DO
      CALL timestop(handle)
   END SUBROUTINE cd_lcg_metrics

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param p ...
!> \param column ...
!> \param ao_shell ...
!> \param group_owner ...
! **************************************************************************************************
   SUBROUTINE store_metric_column(self, p, column, ao_shell, group_owner)
      CLASS(cd_metric_consumer_type), INTENT(INOUT) :: self
      TYPE(post_scf_bandstructure_type), POINTER :: bs_env
      INTEGER, INTENT(IN) :: p
      REAL(dp), INTENT(IN) :: column(:, :)
      INTEGER, INTENT(IN) :: ao_shell(:), group_owner(:)

      INTEGER :: c, q, r, reduce_handle
      REAL(dp) :: values(p)

      bs_env => self%bs_env
      CPASSERT(SIZE(ao_shell) == SIZE(column, 1) .AND. SIZE(group_owner) > 0)
      self%chk(p) = SUM(column*column)
      self%prj(p) = cd_lcg_project(column)
      DO q = 1, p
         values(q) = column(bs_env%cd_lcg_mu(q), bs_env%cd_lcg_nu(q))
      END DO
      CALL timeset('cd_raw_reduce', reduce_handle)
      CALL bs_env%para_env%sum(values)
      CALL timestop(reduce_handle)
      ! Mirror the same evaluated triangle, writing only locally owned matrix elements.
      c = self%local_col(p)
      IF (c > 0) THEN
         DO r = 1, self%nrow
            q = self%rows(r)
            IF (q <= p) self%local_data(r, c) = values(q)
         END DO
      END IF
      r = self%local_row(p)
      IF (r > 0) THEN
         DO c = 1, self%ncol
            q = self%cols(c)
            IF (q < p) self%local_data(r, c) = values(q)
         END DO
      END IF
   END SUBROUTINE store_metric_column

!> \brief Report spectrum on copies; never modify or regularize the supplied metric.
!> \param matrix ...
!> \param bs_env ...
!> \param label ...
! **************************************************************************************************
   SUBROUTINE cd_lcg_spectrum(matrix, bs_env, label)
      TYPE(cp_fm_type)                                   :: matrix
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      CHARACTER(LEN=*), INTENT(IN)                       :: label

      INTEGER                                            :: handle, info
      REAL(dp)                                           :: condition, high, low
      REAL(dp), ALLOCATABLE                              :: eigenvalues(:)
      TYPE(cp_fm_type)                                   :: vectors, work

      CALL timeset('cd_lcg_spectrum', handle)
      CALL cp_fm_create(work, matrix%matrix_struct)
      CALL cp_fm_create(vectors, matrix%matrix_struct)
      CALL cp_fm_to_fm(matrix, work)
      ALLOCATE (eigenvalues(bs_env%n_RI))
      CALL cp_fm_syevd(work, vectors, eigenvalues, info)
      IF (info /= 0) CPABORT("CD_LCG: metric diagnostic diagonalization failed")
      IF (.NOT. ALL(ieee_is_finite(eigenvalues))) CPABORT("CD_LCG: nonfinite metric spectrum")
      low = MINVAL(eigenvalues)
      high = MAXVAL(eigenvalues)
      condition = HUGE(1.0_dp)
      IF (low > 0.0_dp .AND. low >= high/HUGE(1.0_dp)) condition = high/low
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(3A,3(ES24.16,1X))') &
         'CD_LCG| ', label, ' eigen_min eigen_max condition ', low, high, condition
      IF (low <= 0.0_dp) CPABORT("CD_LCG: nonpositive metric eigenvalue")
      CALL cp_fm_release(work)
      CALL cp_fm_release(vectors)
      CALL timestop(handle)
   END SUBROUTINE cd_lcg_spectrum

!> \brief Store T once, distributed within each tensor subgroup, with no dense T replica.
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE cd_lcg_build_T(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      INTEGER                                            :: copy, d, handle, leader, rank
      INTEGER(int_8)                                     :: nblocks, nze
      INTEGER, ALLOCATABLE                               :: all_coord(:, :), colors(:)
      REAL(dp)                                           :: elapsed, occupancy, tick
      TYPE(cd_T_consumer_type)                           :: consumer

      consumer%bs_env => bs_env

      CALL timeset('cd_lcg_build_T', handle)
      ALLOCATE (consumer%chk(SIZE(bs_env%cd_lcg_mu)), consumer%prj(SIZE(bs_env%cd_lcg_mu)), source=0.0_dp)
      CPASSERT(.NOT. bs_env%cd_lcg_T_ready)
      tick = m_walltime()
      CALL dbt_create(bs_env%t_RI_AO__AO, bs_env%cd_lcg_T)
      ALLOCATE (consumer%dist_i(bs_env%n_atom), consumer%dist_j(bs_env%n_atom), &
                consumer%dist_k(bs_env%n_atom))
      ALLOCATE (consumer%reserved(bs_env%n_atom, bs_env%n_atom, bs_env%n_atom), source=.FALSE.)
      CALL dbt_get_info(bs_env%cd_lcg_T, my_ploc=consumer%coord, proc_dist_1=consumer%dist_i, &
                        proc_dist_2=consumer%dist_j, proc_dist_3=consumer%dist_k)
      ! Query each tensor's actual coordinates; never assume Cartesian ranks equal global ranks.
      consumer%nproc = bs_env%para_env%num_pe
      consumer%me = bs_env%para_env%mepos
      consumer%n_sent = 0_int_8
      ! Counts/displacements in the existing alltoallv interface are default integers.
      CPASSERT(INT(bs_env%n_ao, int_8)**2 <= INT(HUGE(0), int_8))
      ALLOCATE (all_coord(consumer%nproc, 3), colors(consumer%nproc))
      DO d = 1, 3
         CALL bs_env%para_env%allgather(consumer%coord(d), all_coord(:, d))
      END DO
      CALL bs_env%para_env%allgather(bs_env%tensor_group_color, colors)
      CPASSERT(bs_env%num_tensor_groups*bs_env%group_size_tensor == consumer%nproc)
      CPASSERT(ALL(colors >= 0) .AND. ALL(colors < bs_env%num_tensor_groups))
      CPASSERT(ALL(all_coord >= 0))
      ALLOCATE (consumer%rank_map(0:MAXVAL(all_coord(:, 1)), 0:MAXVAL(all_coord(:, 2)), &
                                  0:MAXVAL(all_coord(:, 3)), 0:bs_env%num_tensor_groups - 1), source=-1)
      DO rank = 1, consumer%nproc
         leader = consumer%rank_map(all_coord(rank, 1), all_coord(rank, 2), all_coord(rank, 3), colors(rank))
         CPASSERT(leader == -1)
         consumer%rank_map(all_coord(rank, 1), all_coord(rank, 2), all_coord(rank, 3), &
                           colors(rank)) = rank - 1
      END DO
      ALLOCATE (consumer%copy_i(bs_env%n_atom, bs_env%num_tensor_groups), &
                consumer%copy_j(bs_env%n_atom, bs_env%num_tensor_groups), consumer%copy_k(bs_env%n_atom, &
                                                                                          bs_env%num_tensor_groups))
      DO copy = 1, bs_env%num_tensor_groups
         CPASSERT(COUNT(colors == copy - 1) == bs_env%group_size_tensor)
         leader = MINLOC(colors, DIM=1, MASK=(colors == copy - 1)) - 1
         consumer%copy_i(:, copy) = consumer%dist_i
         consumer%copy_j(:, copy) = consumer%dist_j
         consumer%copy_k(:, copy) = consumer%dist_k
         CALL bs_env%para_env%bcast(consumer%copy_i(:, copy), leader)
         CALL bs_env%para_env%bcast(consumer%copy_j(:, copy), leader)
         CALL bs_env%para_env%bcast(consumer%copy_k(:, copy), leader)

         ASSOCIATE (ci => consumer%copy_i, cj => consumer%copy_j, ck => consumer%copy_k, &
                    rmap => consumer%rank_map)
            CPASSERT(ALL(ci(:, copy) >= 0) .AND. ALL(ci(:, copy) <= UBOUND(rmap, 1)))
            CPASSERT(ALL(cj(:, copy) >= 0) .AND. ALL(cj(:, copy) <= UBOUND(rmap, 2)))
            CPASSERT(ALL(ck(:, copy) >= 0) .AND. ALL(ck(:, copy) <= UBOUND(rmap, 3)))
         END ASSOCIATE
         IF (bs_env%tensor_group_color == copy - 1) THEN
            CPASSERT(ALL(consumer%copy_i(:, copy) == consumer%dist_i))
            CPASSERT(ALL(consumer%copy_j(:, copy) == consumer%dist_j))
            CPASSERT(ALL(consumer%copy_k(:, copy) == consumer%dist_k))
         END IF
      END DO
      CALL cd_product_columns(qs_env, bs_env%para_env, bs_env%cd_lcg_mu, bs_env%cd_lcg_nu, &
                              bs_env%ri_metric%cutoff_radius, bs_env%eps_filter, .FALSE., consumer)
      CALL cd_lcg_reduce_diagnostics(bs_env, consumer%chk, consumer%prj)
      CALL bs_env%para_env%sum(consumer%n_sent)
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A,ES24.16)') 'CD_LCG| raw_checksum T ', &
         SUM(consumer%chk)
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A,ES24.16)') 'CD_LCG| raw_projection T ', &
         cd_lcg_combine(consumer%prj)
      CALL dbt_filter(bs_env%cd_lcg_T, bs_env%eps_filter)
      nblocks = dbt_get_num_blocks_total(bs_env%cd_lcg_T)
      nze = dbt_get_nze_total(bs_env%cd_lcg_T)
      occupancy = REAL(nze, dp)/REAL(bs_env%n_RI, dp)/REAL(bs_env%n_ao, dp)**2
      elapsed = m_walltime() - tick
      CALL bs_env%para_env%max(elapsed)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(A,I0,A,I0)') 'CD_LCG| T_tensor_copies ', bs_env%num_tensor_groups, &
            ' group_size ', bs_env%group_size_tensor
         WRITE (bs_env%unit_nr, '(A,I0)') 'CD_LCG| T_raw_exchange_bytes ', 8_int_8*consumer%n_sent
         WRITE (bs_env%unit_nr, '(A,I0,A,ES24.16)') 'CD_LCG| T_blocks ', nblocks, ' occupancy ', occupancy
         WRITE (bs_env%unit_nr, '(A,I0)') 'CD_LCG| T_payload_bytes_per_subgroup ', 8_int_8*nze
         WRITE (bs_env%unit_nr, '(A,F15.6)') 'CD_LCG| T_seconds ', elapsed
      END IF
      IF (nze == 0) CPABORT("CD_LCG: empty filtered T; reduce EPS_FILTER")
      bs_env%cd_lcg_T_ready = .TRUE.
      CALL timestop(handle)
   END SUBROUTINE cd_lcg_build_T

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param p ...
!> \param column ...
!> \param ao_shell ...
!> \param group_owner ...
! **************************************************************************************************
   SUBROUTINE store_T_column(self, p, column, ao_shell, group_owner)
      CLASS(cd_T_consumer_type), INTENT(INOUT) :: self
      TYPE(post_scf_bandstructure_type), POINTER :: bs_env
      INTEGER, INTENT(IN) :: p
      REAL(dp), INTENT(IN) :: column(:, :)
      INTEGER, INTENT(IN) :: ao_shell(:), group_owner(:)

      INTEGER :: c, c_end, dest, dims(3), exchange_handle, i, icopy, ip, j, jp(2), k, kp(2), nvalue, peer, &
                 pos, r, r_end, recv_counts(self%nproc), recv_displs(self%nproc), recv_pos(self%nproc), &
                 send_counts(self%nproc), send_displs(self%nproc), send_pos(self%nproc), src
      INTEGER(int_8) :: n_col, n_pack
      INTEGER :: ao_end(SIZE(ao_shell))
      REAL(dp), ALLOCATABLE :: block(:, :, :), recv_buffer(:), send_buffer(:)

      bs_env => self%bs_env
      self%chk(p) = SUM(column*column)
      self%prj(p) = cd_lcg_project(column)
      i = bs_env%cd_lcg_owner(p)
      ip = p - bs_env%i_RI_start_from_atom(i) + 1
      ! Consecutive AOs in one shell share a source rank. Count shell rectangles and copy row
      ! spans, avoiding a source lookup per AO entry in the globally replicated packing scans.
      r_end = SIZE(ao_shell)
      DO r = SIZE(ao_shell), 1, -1
         IF (ao_shell(r) /= ao_shell(r_end)) r_end = r
         ao_end(r) = r_end
      END DO
      ! Process one tensor copy at a time, bounding exchange storage by one dense AO column.
      ! Keep the source column intact until every subgroup has received its copy.
      DO icopy = 1, bs_env%num_tensor_groups
         CALL timeset('cd_raw_exchange', exchange_handle)
         send_counts = 0
         recv_counts = 0
         DO k = 1, bs_env%n_atom
            kp = [bs_env%i_ao_start_from_atom(k), bs_env%i_ao_end_from_atom(k)]
            DO j = 1, bs_env%n_atom
               dest = tensor_owner(self, i, j, k, icopy)
               jp = [bs_env%i_ao_start_from_atom(j), bs_env%i_ao_end_from_atom(j)]
               c = kp(1)
               DO WHILE (c <= kp(2))
                  c_end = MIN(ao_end(c), kp(2))
                  r = jp(1)
                  DO WHILE (r <= jp(2))
                     r_end = MIN(ao_end(r), jp(2))
                     src = raw_owner(r, c, ao_shell, group_owner)
                     nvalue = (r_end - r + 1)*(c_end - c + 1)
                     IF (src >= 0 .AND. src /= dest) THEN
                        IF (src == self%me) send_counts(dest + 1) = send_counts(dest + 1) + nvalue
                        IF (dest == self%me) recv_counts(src + 1) = recv_counts(src + 1) + nvalue
                     END IF
                     r = r_end + 1
                  END DO
                  c = c_end + 1
               END DO
            END DO
         END DO
         send_displs(1) = 0
         recv_displs(1) = 0
         DO peer = 2, self%nproc
            send_displs(peer) = send_displs(peer - 1) + send_counts(peer - 1)
            recv_displs(peer) = recv_displs(peer - 1) + recv_counts(peer - 1)
         END DO
         n_pack = INT(SUM(send_counts), int_8) + INT(SUM(recv_counts), int_8)
         n_col = INT(SIZE(column, 1), int_8)*INT(SIZE(column, 2), int_8)
         CPASSERT(n_pack <= n_col)
         ALLOCATE (send_buffer(SUM(send_counts)), recv_buffer(SUM(recv_counts)))
         send_pos = send_displs
         DO k = 1, bs_env%n_atom
            kp = [bs_env%i_ao_start_from_atom(k), bs_env%i_ao_end_from_atom(k)]
            DO j = 1, bs_env%n_atom
               dest = tensor_owner(self, i, j, k, icopy)
               IF (dest == self%me) CYCLE
               jp = [bs_env%i_ao_start_from_atom(j), bs_env%i_ao_end_from_atom(j)]
               DO c = kp(1), kp(2)
                  r = jp(1)
                  DO WHILE (r <= jp(2))
                     r_end = MIN(ao_end(r), jp(2))
                     IF (raw_owner(r, c, ao_shell, group_owner) == self%me) THEN
                        nvalue = r_end - r + 1
                        pos = send_pos(dest + 1)
                        send_buffer(pos + 1:pos + nvalue) = column(r:r_end, c)
                        send_pos(dest + 1) = pos + nvalue
                     END IF
                     r = r_end + 1
                  END DO
               END DO
            END DO
         END DO
         CPASSERT(ALL(send_pos == send_displs + send_counts))
         self%n_sent = self%n_sent + INT(SIZE(send_buffer), int_8)
         ! Same typed alltoallv interface as cd_distribute_columns. Only owner payload travels;
         ! self entries are copied locally, and there is no arithmetic on received values.
         CALL bs_env%para_env%alltoall(send_buffer, send_counts, send_displs, recv_buffer, &
                                       recv_counts, recv_displs)
         CALL timestop(exchange_handle)
         recv_pos = recv_displs
         IF (bs_env%tensor_group_color == icopy - 1 .AND. self%dist_i(i) == self%coord(1)) THEN
            ! DBM reservation is an OpenMP collective, including in a one-thread run.
!$OMP PARALLEL NUM_THREADS(1) DEFAULT(NONE)  SHARED(self, bs_env, column, i, ip, ao_shell, ao_end, &
!$OMP group_owner, recv_pos, recv_buffer)  PRIVATE(j,k,jp,kp,dims,block,r,c,src,r_end,nvalue,pos)
            DO k = 1, bs_env%n_atom
               IF (self%dist_k(k) /= self%coord(3)) CYCLE
               kp = [bs_env%i_ao_start_from_atom(k), bs_env%i_ao_end_from_atom(k)]
               DO j = 1, bs_env%n_atom
                  IF (self%dist_j(j) /= self%coord(2)) CYCLE
                  jp = [bs_env%i_ao_start_from_atom(j), bs_env%i_ao_end_from_atom(j)]
                  dims = [bs_env%sizes_RI(i), bs_env%sizes_AO(j), bs_env%sizes_AO(k)]
                  ALLOCATE (block(dims(1), dims(2), dims(3)), source=0.0_dp)
                  DO c = kp(1), kp(2)
                     r = jp(1)
                     DO WHILE (r <= jp(2))
                        r_end = MIN(ao_end(r), jp(2))
                        src = raw_owner(r, c, ao_shell, group_owner)
                        nvalue = r_end - r + 1
                        IF (src == self%me) THEN
                           block(ip, r - jp(1) + 1:r_end - jp(1) + 1, c - kp(1) + 1) = column(r:r_end, c)
                        ELSE IF (src >= 0) THEN
                           pos = recv_pos(src + 1)
                           block(ip, r - jp(1) + 1:r_end - jp(1) + 1, c - kp(1) + 1) = &
                              recv_buffer(pos + 1:pos + nvalue)
                           recv_pos(src + 1) = pos + nvalue
                        END IF
                        r = r_end + 1
                     END DO
                  END DO
                  ! Match the old multi-rank sum's positive zero, including in mixed nonzero blocks.
                  IF (bs_env%para_env%num_pe > 1) THEN
                     WHERE (block(ip, :, :) == 0.0_dp) block(ip, :, :) = 0.0_dp
                  END IF
                  IF (ANY(block(ip, :, :) /= 0.0_dp)) THEN
                     IF (.NOT. self%reserved(i, j, k)) THEN
                        CALL dbt_reserve_blocks(bs_env%cd_lcg_T, [i], [j], [k])
                        self%reserved(i, j, k) = .TRUE.
                     END IF
                     CALL dbt_put_block(bs_env%cd_lcg_T, [i, j, k], dims, block, summation=.TRUE.)
                  END IF
                  DEALLOCATE (block)
               END DO
            END DO
!$OMP END PARALLEL
         END IF
         CPASSERT(ALL(recv_pos == recv_displs + recv_counts))
         DEALLOCATE (send_buffer, recv_buffer)
      END DO
   END SUBROUTINE store_T_column

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param i ...
!> \param j ...
!> \param k ...
!> \param icopy ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION tensor_owner(self, i, j, k, icopy) RESULT(owner)
      CLASS(cd_T_consumer_type), INTENT(IN) :: self
      INTEGER, INTENT(IN) :: i, j, k, icopy

      owner = self%rank_map(self%copy_i(i, icopy), self%copy_j(j, icopy), self%copy_k(k, icopy), icopy - 1)
      CPASSERT(owner >= 0)
   END FUNCTION tensor_owner

!> \brief Translate atom bounds to function bounds and copy from persistent T.
!> \param bs_env ...
!> \param tensor ...
!> \param atoms_RI ...
!> \param atoms_AO_1 ...
!> \param atoms_AO_2 ...
! **************************************************************************************************
   SUBROUTINE cd_lcg_get_T(bs_env, tensor, atoms_RI, atoms_AO_1, atoms_AO_2)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbt_type)                                     :: tensor
      INTEGER, INTENT(IN), OPTIONAL                      :: atoms_RI(2), atoms_AO_1(2), atoms_AO_2(2)

      INTEGER                                            :: atoms(2, 3), bounds(2, 3), d, handle

      CALL timeset('cd_lcg_get_T', handle)
      CPASSERT(bs_env%cd_lcg_T_ready)
      atoms(1, :) = 1
      atoms(2, :) = bs_env%n_atom
      IF (PRESENT(atoms_RI)) atoms(:, 1) = atoms_RI
      IF (PRESENT(atoms_AO_1)) atoms(:, 2) = atoms_AO_1
      IF (PRESENT(atoms_AO_2)) atoms(:, 3) = atoms_AO_2
      CALL dbt_clear(tensor)
      IF (ANY(atoms(1, :) > atoms(2, :))) THEN
         CALL timestop(handle)
         RETURN
      END IF
      CPASSERT(ALL(atoms >= 1) .AND. ALL(atoms <= bs_env%n_atom))
      bounds(:, 1) = [bs_env%i_RI_start_from_atom(atoms(1, 1)), bs_env%i_RI_end_from_atom(atoms(2, 1))]
      DO d = 2, 3
         bounds(:, d) = [bs_env%i_ao_start_from_atom(atoms(1, d)), bs_env%i_ao_end_from_atom(atoms(2, d))]
      END DO
      CALL dbt_copy(bs_env%cd_lcg_T, tensor, bounds=bounds)
      CALL timestop(handle)
   END SUBROUTINE cd_lcg_get_T

!> \brief Matrix-array contract of the Gaussian 2c provider, with zero imaginary part.
!> \param bs_env ...
!> \param matrices ...
!> \param full_coulomb ...
!> \param complex_pair ...
! **************************************************************************************************
   SUBROUTINE cd_lcg_matrix_pair(bs_env, matrices, full_coulomb, complex_pair)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, INTENT(OUT)         :: matrices(:, :)
      LOGICAL, INTENT(IN)                                :: full_coulomb, complex_pair

      INTEGER                                            :: handle, j, n

      CALL timeset('cd_lcg_matrix_pair', handle)
      CPASSERT(ALL(bs_env%periodic == 0))
      n = 1
      IF (complex_pair) n = 2
      ALLOCATE (matrices(1, n))
      DO j = 1, n
         CALL cp_fm_create(matrices(1, j), bs_env%fm_RI_RI%matrix_struct)
         CALL cp_fm_set_all(matrices(1, j), 0.0_dp)
      END DO
      IF (full_coulomb) THEN
         CALL cp_fm_to_fm(bs_env%cd_lcg_V, matrices(1, 1))
      ELSE
         CALL cp_fm_to_fm(bs_env%cd_lcg_M, matrices(1, 1))
      END IF
      CALL timestop(handle)
   END SUBROUTINE cd_lcg_matrix_pair
!> \brief Source rank of either orientation of an AO pair; -1 denotes an omitted shell group.
!> \param r ...
!> \param c ...
!> \param ao_shell ...
!> \param group_owner ...
!> \return ...
! **************************************************************************************************
   PURE INTEGER FUNCTION raw_owner(r, c, ao_shell, group_owner) RESULT(owner)
      INTEGER, INTENT(IN)                                :: r, c, ao_shell(:), group_owner(:)

      INTEGER                                            :: a, b, g

      a = MIN(ao_shell(r), ao_shell(c))
      b = MAX(ao_shell(r), ao_shell(c))
      g = INT(INT(b, int_8)*INT(b - 1, int_8)/2_int_8) + a
      owner = group_owner(g)
   END FUNCTION raw_owner

!> \brief Sum disjoint raw-column diagnostics once, without counting tensor-subgroup copies.
!> \param bs_env ...
!> \param chk ...
!> \param prj ...
! **************************************************************************************************
   SUBROUTINE cd_lcg_reduce_diagnostics(bs_env, chk, prj)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(dp), INTENT(INOUT)                            :: chk(:), prj(:)

      INTEGER                                            :: handle

      CALL timeset('cd_raw_reduce', handle)
      CALL bs_env%para_env%sum(chk)
      CALL bs_env%para_env%sum(prj)
      CALL timestop(handle)
   END SUBROUTINE cd_lcg_reduce_diagnostics

! **************************************************************************************************
!> \brief Fixed pseudo-random projection of one raw column; changes under sign flips and index permutations.
!> \param column ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION cd_lcg_project(column) RESULT(res)
      REAL(dp), INTENT(IN)                               :: column(:, :)
      REAL(dp)                                           :: res

      INTEGER                                            :: c, r

      res = 0.0_dp
      DO c = 1, SIZE(column, 2)
         DO r = 1, SIZE(column, 1)
            res = res + column(r, c)*COS(0.6180339887498949_dp*REAL(r, dp) + 0.4142135623730950_dp*REAL(c, dp))
         END DO
      END DO
   END FUNCTION cd_lcg_project

! **************************************************************************************************
!> \brief Combine per-pivot projections in pivot order with pivot-dependent weights.
!> \param prj ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION cd_lcg_combine(prj) RESULT(res)
      REAL(dp), INTENT(IN)                               :: prj(:)
      REAL(dp)                                           :: res

      INTEGER                                            :: p

      res = 0.0_dp
      DO p = 1, SIZE(prj)
         res = res + prj(p)*COS(0.7320508075688772_dp*REAL(p, dp))
      END DO
   END FUNCTION cd_lcg_combine

END MODULE gw_cd_lcg
