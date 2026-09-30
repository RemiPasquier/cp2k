!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Distributed full-Coulomb Cholesky of unique, unscaled contracted ORB pairs.
!>        Row-distributed producer and pure-L provider for the quartic RPA/GW/BSE path.
MODULE mp2_cd_producer
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_distribution_get, dbcsr_distribution_type, dbcsr_get_info, &
        dbcsr_get_stored_coordinates, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_reserve_blocks, &
        dbcsr_type, dbcsr_type_symmetric
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE gamma,                           ONLY: init_md_ftable
   USE hfx_energy_potential,            ONLY: coulomb4
   USE hfx_types,                       ONLY: &
        hfx_basis_type, hfx_create_basis_types, hfx_pgf_list, hfx_pgf_product_list, &
        hfx_potential_type, hfx_release_basis_types, hfx_screen_coeff_type, hfx_type, &
        init_t_c_g0_lmax, log_zero
   USE ieee_arithmetic,                 ONLY: ieee_is_finite,&
                                              ieee_positive_inf,&
                                              ieee_value
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_truncated
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE libint_wrapper,                  ONLY: cp_libint_cleanup_eri,&
                                              cp_libint_init_eri,&
                                              cp_libint_set_contrdepth,&
                                              cp_libint_static_init,&
                                              cp_libint_t,&
                                              libint_max_am
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_ri_libint,                   ONLY: prepare_integral_calc
   USE mp2_types,                       ONLY: mp2_type
   USE omp_lib,                         ONLY: omp_get_max_threads
   USE orbital_pointers,                ONLY: nco,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE t_c_g0,                          ONLY: get_lmax_init,&
                                              init
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: cd_build_aux, cd_distribute_columns, cd_get_rank, cd_fill_column, cd_release_aux, &
             cd_reserve_full_support, cd_select_products, cd_product_columns
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_cd_producer'
   INTEGER, PARAMETER :: cd_batch_size = 64
   ! Block pivoting qualifies residual diagonals above this fraction of the current maximum.
   ! RI-RPA/GW uses the Cholesky vectors themselves: the residual bound controls their error, and a
   ! near-greedy qualification gives the fewest vectors. Large-cell GW uses the pivot products as
   ! auxiliary functions in a truncated-Coulomb fit, which the residual bound does not control; the
   ! wider qualification keeps more products and is needed for converged quasiparticle energies.
   REAL(dp), PARAMETER :: cd_qualification_vectors = 0.3_dp, cd_qualification_products = 0.01_dp

   TYPE cd_controls
      REAL(dp) :: eps_cholesky = 1.0e-5_dp, eps_schwarz = 1.0e-10_dp, &
                  qualification = cd_qualification_vectors
   END TYPE cd_controls

   TYPE cd_shell
      INTEGER :: atom, kind, iset, l, n, first
   END TYPE cd_shell

   TYPE cd_pair
      ! group identifies two cd_shell records, hence both atom/set triples.
      INTEGER :: group, a, b, mu, nu
   END TYPE cd_pair

   TYPE cd_raw_column
      REAL(dp), ALLOCATABLE :: values(:)
   END TYPE cd_raw_column

   TYPE cd_l_basis
      INTEGER :: npgf_l = 0
      REAL(dp), ALLOCATABLE :: zet_l(:), sphi_ext(:, :, :)
   END TYPE cd_l_basis

   ! Row-distributed until cd_distribute_columns transposes it for the consumer.
   REAL(dp), ALLOCATABLE, SAVE :: retained_factor(:, :)
   INTEGER, ALLOCATABLE, SAVE :: retained_owner(:), retained_row(:)
   INTEGER, SAVE :: retained_first_column = 0, retained_rank = 0, retained_nao = 0

   ! Consumer of one raw product column (pivot ip) in cd_product_columns.
   ABSTRACT INTERFACE
! **************************************************************************************************
!> \brief ...
!> \param ip ...
!> \param column ...
!> \param ao_shell ...
!> \param group_owner ...
! **************************************************************************************************
      SUBROUTINE cd_column_consumer(ip, column, ao_shell, group_owner)
         IMPORT :: dp
         INTEGER, INTENT(IN) :: ip, ao_shell(:), group_owner(:)
         REAL(dp), INTENT(IN) :: column(:, :)
      END SUBROUTINE cd_column_consumer
   END INTERFACE

CONTAINS

!> \brief Molecular ORB scratch for the standalone selector and raw-product oracle.
!> \param qs_env ...
!> \param x ...
!> \param basis ...
!> \param particles ...
!> \param kind_of ...
!> \param max_set ...
!> \param natom ...
!> \param nsgf_max ...
!> \param max_pgf ...
!> \param integrals ...
!> \param work ...
!> \param work2 ...
!> \param buffer1 ...
!> \param buffer2 ...
!> \param primitives ...
!> \param pgf_ij ...
!> \param pgf_kl ...
!> \param products ...
!> \param nimages ...
! **************************************************************************************************
   SUBROUTINE cd_prepare_molecular(qs_env, x, basis, particles, kind_of, max_set, natom, &
                                   nsgf_max, max_pgf, integrals, work, work2, buffer1, buffer2, &
                                   primitives, pgf_ij, pgf_kl, products, nimages)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_type), POINTER                            :: x
      TYPE(hfx_basis_type), POINTER                      :: basis(:)
      TYPE(particle_type), POINTER                       :: particles(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: kind_of(:)
      INTEGER, INTENT(OUT)                               :: max_set, natom, nsgf_max, max_pgf
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: integrals(:), work(:), work2(:), &
                                                            buffer1(:), buffer2(:), primitives(:)
      TYPE(hfx_pgf_list), ALLOCATABLE, INTENT(OUT)       :: pgf_ij(:), pgf_kl(:)
      TYPE(hfx_pgf_product_list), ALLOCATABLE, &
         INTENT(OUT)                                     :: products(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: nimages(:)

      INTEGER                                            :: i, ikind, iset, ncos_max
      TYPE(atomic_kind_type), POINTER                    :: kinds(:)

      CALL get_qs_env(qs_env, atomic_kind_set=kinds, particle_set=particles)
      CALL get_atomic_kind_set(kinds, kind_of=kind_of)
      basis => x%basis_parameter
      max_set = x%basis_info%max_set
      natom = SIZE(particles)
      CALL init_md_ftable(4*x%basis_info%max_am)
      ncos_max = 0
      nsgf_max = 0
      max_pgf = 0
      DO ikind = 1, SIZE(basis)
         DO iset = 1, basis(ikind)%nset
            ncos_max = MAX(ncos_max, ncoset(basis(ikind)%lmax(iset)))
            nsgf_max = MAX(nsgf_max, basis(ikind)%nsgf(iset))
            max_pgf = MAX(max_pgf, basis(ikind)%npgf(iset))
         END DO
      END DO
      ALLOCATE (integrals(nsgf_max**4), primitives(nsgf_max**4))
      ALLOCATE (work(ncos_max**4), work2(ncos_max**4), buffer1(ncos_max**4), buffer2(ncos_max**4))
      ALLOCATE (pgf_ij(max_pgf**2), pgf_kl(max_pgf**2), products(1), nimages(max_pgf**2))
      DO i = 1, max_pgf**2
         ALLOCATE (pgf_ij(i)%image_list(1), pgf_kl(i)%image_list(1))
      END DO
   END SUBROUTINE cd_prepare_molecular

!> \brief Run exactly the existing full-Coulomb selector, returning products in pivot order.
!> \param qs_env ...
!> \param para_env ...
!> \param eps_cholesky ...
!> \param unit_nr ...
!> \param mu ...
!> \param nu ...
!> \param owner ...
! **************************************************************************************************
   SUBROUTINE cd_select_products(qs_env, para_env, eps_cholesky, unit_nr, mu, nu, owner)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(dp), INTENT(IN)                               :: eps_cholesky
      INTEGER, INTENT(IN)                                :: unit_nr
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: mu(:), nu(:), owner(:)

      INTEGER                                            :: handle
      TYPE(cd_controls)                                  :: controls

      CALL timeset('cd_select_products', handle)
      IF (omp_get_max_threads() /= 1) CPABORT("CD_LCG: requires one OpenMP thread")
      IF (.NOT. ieee_is_finite(eps_cholesky)) CPABORT("CD_LCG: nonfinite EPS_CHOLESKY")
      IF (eps_cholesky <= 0.0_dp .OR. eps_cholesky >= 1.0_dp) THEN
         CPABORT("CD_LCG: require 0 < EPS_CHOLESKY < 1")
      END IF
      controls%eps_cholesky = eps_cholesky
      controls%qualification = cd_qualification_products
      CALL cd_build_once(qs_env, controls, para_env, unit_nr, .FALSE., mu=mu, nu=nu, owner=owner)
      CALL timestop(handle)
   END SUBROUTINE cd_select_products

!> \brief Pack exact nonzero primitive support of one (kind, set, l), without renormalizing.
!> \param basis Original HFX basis and spherical contraction slices.
!> \param orb_basis Original normalized Cartesian contractions feeding sphi.
!> \param iset Original set index.
!> \param l Angular momentum to retain.
!> \param packed Private primitive data, shared by all atoms of this kind.
! **************************************************************************************************
   SUBROUTINE cd_pack_l_basis(basis, orb_basis, iset, l, packed)
      TYPE(hfx_basis_type), INTENT(IN)                   :: basis
      TYPE(gto_basis_set_type), INTENT(IN)               :: orb_basis
      INTEGER, INTENT(IN)                                :: iset, l
      TYPE(cd_l_basis), INTENT(OUT)                      :: packed

      INTEGER                                            :: first_cgf, ipgf, ishell, last_cgf, &
                                                            offset, p
      LOGICAL, ALLOCATABLE                               :: keep(:)

      ALLOCATE (keep(basis%npgf(iset)))
      keep = .FALSE.
      DO ishell = 1, basis%nshell(iset)
         IF (basis%nl(ishell, iset) /= l) CYCLE
         first_cgf = orb_basis%first_cgf(ishell, iset)
         last_cgf = orb_basis%last_cgf(ishell, iset)
         DO ipgf = 1, basis%npgf(iset)
            offset = (ipgf - 1)*ncoset(basis%lmax(iset)) + ncoset(l - 1)
            ! cphi = norm_cgf * gcc, before the Cartesian-to-spherical transform.
            ! Test every contraction of this l; never threshold transformed coefficients.
            keep(ipgf) = keep(ipgf) .OR. &
                         ANY(orb_basis%cphi(offset + 1:offset + nco(l), first_cgf:last_cgf) /= 0.0_dp)
         END DO
      END DO
      packed%npgf_l = COUNT(keep)
      ALLOCATE (packed%zet_l(packed%npgf_l))
      ! Preserve the l=0 origin: coulomb4 addresses this plane as l+1.
      ALLOCATE (packed%sphi_ext(SIZE(basis%sphi_ext, 1), 0:UBOUND(basis%sphi_ext, 2), packed%npgf_l))
      packed%sphi_ext = 0.0_dp
      p = 0
      DO ipgf = 1, basis%npgf(iset)
         IF (.NOT. keep(ipgf)) CYCLE
         p = p + 1
         packed%zet_l(p) = basis%zet(ipgf, iset)
         packed%sphi_ext(:, l, p) = basis%sphi_ext(:, l, ipgf, iset)
      END DO
      CPASSERT(p == packed%npgf_l)
   END SUBROUTINE cd_pack_l_basis

!> \brief Release the distributed provider after its AO columns have been consumed.
! **************************************************************************************************
   SUBROUTINE cd_release_aux()
      IF (ALLOCATED(retained_factor)) DEALLOCATE (retained_factor)
      IF (ALLOCATED(retained_owner)) DEALLOCATE (retained_owner)
      IF (ALLOCATED(retained_row)) DEALLOCATE (retained_row)
      retained_first_column = 0
      retained_rank = 0
      retained_nao = 0
   END SUBROUTINE cd_release_aux

!> \brief Auxiliary dimension; no Gaussian basis is involved.
!> \return ...
! **************************************************************************************************
   FUNCTION cd_get_rank() RESULT(k)
      INTEGER                                            :: k

      CPASSERT(ALLOCATED(retained_factor))
      k = retained_rank
      CPASSERT(k > 0)
   END FUNCTION cd_get_rank

!> \brief Unpack one unscaled unique-pair factor column into a full symmetric AO matrix.
!> \param column Global auxiliary index.
!> \param matrix Dense-block symmetric ORB template; only local DBCSR blocks are visited.
! **************************************************************************************************
   SUBROUTINE cd_fill_column(column, matrix)
      INTEGER, INTENT(IN)                                :: column
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CHARACTER                                          :: matrix_type
      INTEGER                                            :: col, col_offset, i, j, mu, nc, nr, nu, &
                                                            pair, row, row_offset, source
      LOGICAL                                            :: transposed
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iter

      CPASSERT(ALLOCATED(retained_factor))
      CPASSERT(column >= retained_first_column)
      CPASSERT(column < retained_first_column + SIZE(retained_factor, 2))
      CALL dbcsr_get_info(matrix, nfullrows_total=nr, nfullcols_total=nc, matrix_type=matrix_type)
      CPASSERT(nr == retained_nao .AND. nc == retained_nao)
      CPASSERT(matrix_type == dbcsr_type_symmetric)
      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, row_offset=row_offset, &
                                        col_offset=col_offset, transposed=transposed)
         IF (transposed) THEN
            i = row_offset
            row_offset = col_offset
            col_offset = i
         END IF
         DO j = 1, SIZE(block, 2)
            DO i = 1, SIZE(block, 1)
               mu = MIN(row_offset + i - 1, col_offset + j - 1)
               nu = MAX(row_offset + i - 1, col_offset + j - 1)
               pair = INT(INT(nu, int_8)*INT(nu - 1, int_8)/2_int_8 + INT(mu, int_8))
               source = retained_row(pair)
               block(i, j) = 0.0_dp
               IF (source > 0) block(i, j) = retained_factor(source, column - retained_first_column + 1)
            END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)
   END SUBROUTINE cd_fill_column

! **************************************************************************************************
!> \brief Transpose row-owned L to the auxiliary-column distribution used by consumer subgroups.
!>        Every rank in a subgroup receives the same columns; AO rows are reconstructed in global
!>        screened-pair order and cd_fill_column subsequently touches only locally owned blocks.
!> \param para_env Global communicator.
!> \param first_column First auxiliary column required by this rank's subgroup.
!> \param last_column Last auxiliary column required by this rank's subgroup.
! **************************************************************************************************
   SUBROUTINE cd_distribute_columns(para_env, first_column, last_column)
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      INTEGER, INTENT(IN)                                :: first_column, last_column

      INTEGER                                            :: dest, global_row, i, j, ncolumn, nlocal, &
                                                            nrow, pos, source
      INTEGER(int_8)                                     :: count8
      INTEGER, ALLOCATABLE                               :: ends(:), recv_counts(:), recv_displs(:), &
                                                            row_counts(:), send_counts(:), &
                                                            send_displs(:), starts(:)
      REAL(dp), ALLOCATABLE                              :: recv_buffer(:), send_buffer(:), &
                                                            transposed(:, :)

      CPASSERT(ALLOCATED(retained_factor))
      CPASSERT(ALLOCATED(retained_owner))
      CPASSERT(first_column >= 1 .AND. first_column <= retained_rank + 1)
      CPASSERT(last_column >= 0 .AND. last_column <= retained_rank)
      CPASSERT(last_column >= first_column - 1)
      CPASSERT(retained_first_column == 1)
      nlocal = SIZE(retained_factor, 1)
      nrow = SIZE(retained_owner)
      ncolumn = last_column - first_column + 1
      ALLOCATE (starts(para_env%num_pe), ends(para_env%num_pe), row_counts(para_env%num_pe))
      CALL para_env%allgather(first_column, starts)
      CALL para_env%allgather(last_column, ends)
      CALL para_env%allgather(nlocal, row_counts)
      CPASSERT(SUM(row_counts) == nrow)

      ALLOCATE (send_counts(para_env%num_pe), recv_counts(para_env%num_pe))
      ALLOCATE (send_displs(para_env%num_pe), recv_displs(para_env%num_pe))
      DO dest = 1, para_env%num_pe
         count8 = INT(nlocal, int_8)*INT(ends(dest) - starts(dest) + 1, int_8)
         CPASSERT(count8 <= INT(HUGE(0), int_8))
         send_counts(dest) = INT(count8)
         count8 = INT(row_counts(dest), int_8)*INT(ncolumn, int_8)
         CPASSERT(count8 <= INT(HUGE(0), int_8))
         recv_counts(dest) = INT(count8)
      END DO
      send_displs(1) = 0
      recv_displs(1) = 0
      DO dest = 2, para_env%num_pe
         send_displs(dest) = send_displs(dest - 1) + send_counts(dest - 1)
         recv_displs(dest) = recv_displs(dest - 1) + recv_counts(dest - 1)
      END DO
      ALLOCATE (send_buffer(SUM(send_counts)), recv_buffer(SUM(recv_counts)))
      pos = 0
      DO dest = 1, para_env%num_pe
         DO j = starts(dest), ends(dest)
            DO i = 1, nlocal
               pos = pos + 1
               send_buffer(pos) = retained_factor(i, j)
            END DO
         END DO
      END DO
      CPASSERT(pos == SIZE(send_buffer))
      CALL para_env%alltoall(send_buffer, send_counts, send_displs, recv_buffer, &
                             recv_counts, recv_displs)

      ALLOCATE (transposed(nrow, ncolumn))
      transposed = 0.0_dp
      DO source = 1, para_env%num_pe
         pos = recv_displs(source)
         DO j = 1, ncolumn
            DO global_row = 1, nrow
               IF (retained_owner(global_row) /= source - 1) CYCLE
               pos = pos + 1
               transposed(global_row, j) = recv_buffer(pos)
            END DO
         END DO
         CPASSERT(pos == recv_displs(source) + recv_counts(source))
      END DO
      CALL MOVE_ALLOC(transposed, retained_factor)
      retained_first_column = first_column
      DEALLOCATE (retained_owner, starts, ends, row_counts, send_counts, recv_counts)
      DEALLOCATE (send_displs, recv_displs, send_buffer, recv_buffer)
   END SUBROUTINE cd_distribute_columns

! **************************************************************************************************
!> \brief Reserve every atom-pair block a CD column can touch, honouring the symmetric storage
!>        convention (row <= col, as cp_dbcsr_alloc_block_from_nbl); reserving both triangles of a
!>        symmetric matrix stores each off-diagonal pair twice and double-counts it downstream.
!> \param matrix symmetric AO template
! **************************************************************************************************
   SUBROUTINE cd_reserve_full_support(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CHARACTER                                          :: matrix_type
      INTEGER                                            :: i, j, k, mynode, nblk, owner
      INTEGER, ALLOCATABLE                               :: cols(:), rows(:)
      TYPE(dbcsr_distribution_type)                      :: dist

      ! Reserve each upper-triangle block only on the rank DBCSR stores it on. Local row/col
      ! lists are not sufficient: symmetric storage may place block (i,j) by its transpose.
      CALL dbcsr_get_info(matrix, nblkrows_total=nblk, matrix_type=matrix_type, distribution=dist)
      CPASSERT(matrix_type == dbcsr_type_symmetric)
      CALL dbcsr_distribution_get(dist, mynode=mynode)
      ALLOCATE (rows(nblk*(nblk + 1)/2), cols(nblk*(nblk + 1)/2))
      k = 0
      DO j = 1, nblk
         DO i = 1, j
            CALL dbcsr_get_stored_coordinates(matrix, i, j, owner)
            IF (owner /= mynode) CYCLE
            k = k + 1
            rows(k) = i
            cols(k) = j
         END DO
      END DO
      CALL dbcsr_reserve_blocks(matrix, rows=rows(1:k), cols=cols(1:k))
      DEALLOCATE (rows, cols)
   END SUBROUTINE cd_reserve_full_support

! **************************************************************************************************
!> \brief Build and retain the Cholesky factor for the quartic correlation consumer.
!> \param qs_env Quickstep environment.
!> \param mp2_env Correlation controls.
!> \param para_env Global communicator.
!> \param unit_nr Report output unit.
! **************************************************************************************************
   SUBROUTINE cd_build_aux(qs_env, mp2_env, para_env, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cd_build_aux'

      INTEGER                                            :: handle
      TYPE(cd_controls)                                  :: controls

      CALL timeset(routineN, handle)
      CALL cd_release_aux()
      IF (omp_get_max_threads() /= 1) CPABORT("CD_PRODUCER: requires one OpenMP thread")
      IF (mp2_env%potential_parameter%potential_type /= do_potential_coulomb) THEN
         CPABORT("CD_PRODUCER: require full Coulomb INTERACTION_POTENTIAL")
      END IF
      controls = cd_controls(mp2_env%cholesky%eps_cholesky, mp2_env%cholesky%eps_schwarz, &
                             cd_qualification_vectors)
      CALL cd_build_once(qs_env, controls, para_env, unit_nr, .TRUE., mp2_env=mp2_env)
      IF (unit_nr > 0 .AND. para_env%is_source()) CALL m_flush(unit_nr)
      CALL timestop(handle)
   END SUBROUTINE cd_build_aux

! **************************************************************************************************
!> \brief Exact diagonal screening and batched scalar pivoting (Folkestad et al., 2019).
!> \param qs_env Quickstep environment.
!> \param controls Screening thresholds from correlation input or standalone selection.
!> \param para_env Global communicator.
!> \param unit_nr Report output unit.
!> \param retain Keep the final factor for the consumer.
!> \param mp2_env Optional legacy preparation context; absent for standalone selection.
!> \param mu Canonical first AO for each selected pivot, in selection order.
!> \param nu Canonical second AO for each selected pivot.
!> \param owner Atom owning the first AO.
! **************************************************************************************************
   SUBROUTINE cd_build_once(qs_env, controls, para_env, unit_nr, retain, mp2_env, mu, nu, owner)
      TYPE(qs_environment_type), POINTER :: qs_env
      TYPE(cd_controls), INTENT(IN) :: controls
      TYPE(mp2_type), OPTIONAL, INTENT(INOUT) :: mp2_env
      INTEGER, ALLOCATABLE, OPTIONAL, INTENT(OUT) :: mu(:), nu(:), owner(:)
      TYPE(mp_para_env_type), INTENT(IN) :: para_env
      INTEGER, INTENT(IN) :: unit_nr
      LOGICAL, INTENT(IN) :: retain

      TYPE(cell_type), POINTER :: cell
      TYPE(qs_kind_type), DIMENSION(:), POINTER :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
      TYPE(hfx_type), DIMENSION(:, :), POINTER :: saved_x, owned_x
      TYPE(hfx_type), POINTER :: actual_x_data
      TYPE(hfx_basis_type), DIMENSION(:), POINTER :: basis
      TYPE(gto_basis_set_type), POINTER :: orb_basis
      TYPE(cd_l_basis), ALLOCATABLE, TARGET :: packed_basis(:, :, :)
      TYPE(hfx_potential_type) :: potential
      TYPE(cp_libint_t) :: lib
      TYPE(hfx_pgf_list), ALLOCATABLE :: pgf_ij(:), pgf_kl(:)
      TYPE(hfx_pgf_product_list), ALLOCATABLE :: products(:)
      TYPE(hfx_screen_coeff_type), POINTER :: screen_set(:, :, :, :), screen_kind(:, :), &
                                              screen_pgf(:, :, :, :, :, :), radii(:, :, :, :, :, :)
      TYPE(hfx_screen_coeff_type), POINTER :: neutral(:, :)
      REAL(dp), POINTER :: p_work(:)
      REAL(dp), ALLOCATABLE :: integrals(:), work(:), work2(:), buffer1(:), buffer2(:), primitives(:), &
                               max_contraction(:, :), d(:), d_all(:), scores(:), factor(:, :), grown(:, :), &
                               panel(:, :), qfactor(:, :), qblock(:, :), pivot_row(:), qdiag(:), &
                               lsum(:), lsq(:)
      TYPE(cd_shell), ALLOCATABLE :: shells(:)
      TYPE(cd_pair), ALLOCATABLE :: pairs(:)
      TYPE(cd_raw_column), ALLOCATABLE :: cache(:)
      INTEGER, ALLOCATABLE :: kind_of(:), nimages(:), bra(:), ket(:), first(:), q(:), groups(:), &
                              pivots(:), missing(:), slots(:), group_owner(:), pair_owner(:), local_index(:)
      INTEGER(int_8), ALLOCATABLE :: group_weight(:)
      LOGICAL, ALLOCATABLE :: selected(:), touched(:), active(:), group_active(:)
      INTEGER :: max_set, natom, nsgf_max, nspins, max_pgf, periodic(3), nao, nshell, ngroup, ntotal, n, &
                 atom, ikind, iset, ishell, jshell, g, a, b, i, j, m, oldfirst, oldlast, k, capacity, &
                 newcap, nq, best, pivot, batch, nleft, ntouched, next_group, group_size, panel_limit, l, &
                 nmiss, cache_stat, nlocal, nlocal_full, max_nlocal, local_row, bad, rows_min, rows_max, &
                 rows_active_final, groups_active_final, low, high, mid
      ! Bound the local panel, gathered factor rows, Q block, and panel vectors to 2 GiB.
      INTEGER(int_8), PARAMETER :: panel_bytes = 2_int_8*1024_int_8**3
      ! Raw-column payload budget, separate from the panel and retained factor.
      INTEGER(int_8), PARAMETER :: raw_cache_budget = 2_int_8*1024_int_8**3
      ! Covers rank-dependent O(k*epsilon*dmax) Schur-update noise without merging distinct scales.
      REAL(dp), PARAMETER :: pivot_tie_tol = 1.0e-12_dp
      REAL(dp) :: tie_band
      INTEGER(int_8) :: npair8, quartets_unique, integrals_computed, columns_evaluated, &
                        cache_bytes, column_bytes, columns_cache_hits, panels, panel_words, linear_words, &
                        row_updates_skipped, quartets_skipped
      REAL(dp) :: helper_eps, helper_logeps, dmax_initial, dmax, bound, t_start, tick, t_diag, t_col, t_la, &
                  t_total, rel_final, roundoff, best_d, dmax_start, &
                  quartet_time, quartet_time_min, quartet_time_max, quartet_tick, prune_bound
      LOGICAL :: do_periodic, cache_full, do_prune

      t_start = m_walltime()
      t_col = 0.0_dp
      t_la = 0.0_dp
      quartets_unique = 0_int_8
      integrals_computed = 0_int_8
      columns_evaluated = 0_int_8
      cache_bytes = 0_int_8
      columns_cache_hits = 0_int_8
      panels = 0_int_8
      row_updates_skipped = 0_int_8
      quartets_skipped = 0_int_8
      quartet_time = 0.0_dp
      cache_full = .FALSE.
      ! The retained factor is consumed with the original row map.
      do_prune = .NOT. retain
      NULLIFY (cell, qs_kind_set, particle_set, saved_x, owned_x, actual_x_data, basis, p_work)
      CALL get_qs_env(qs_env, cell=cell, qs_kind_set=qs_kind_set)
      CALL get_cell(cell, periodic=periodic)
      IF (ANY(periodic /= 0)) CPABORT("CD_PRODUCER: only nonperiodic molecules are supported")

      ! prepare_integral_calc BORROWS x_data; it does not create it. Supply a private
      ! minimal ORB context even for HF/ADMM inputs, without reading an HF section.
      ! Restore the environment immediately after preparation. No SCF cache is touched.
      ALLOCATE (owned_x(1, 1))
      actual_x_data => owned_x(1, 1)
      CALL hfx_create_basis_types(actual_x_data%basis_parameter, actual_x_data%basis_info, qs_kind_set, 'ORB')
      IF (actual_x_data%basis_info%max_am > libint_max_am) THEN
         CPABORT("CD_PRODUCER: ORB angular momentum exceeds the libint build limit")
      END IF
      CALL cp_libint_static_init()
      CALL cp_libint_init_eri(actual_x_data%lib, actual_x_data%basis_info%max_am)
      CALL cp_libint_set_contrdepth(actual_x_data%lib, 1)
      ALLOCATE (actual_x_data%neighbor_cells(1))
      actual_x_data%periodic_parameter%do_periodic = .FALSE.
      actual_x_data%screening_parameter%eps_schwarz = 0.0_dp
      ! No fitted HFX bounds are needed. Initialized zero tables disable logarithmic
      ! screening; only their pointer associations are used by the preparation helper.
      ALLOCATE (actual_x_data%screen_funct_coeffs_set(1, 1, 1, 1))
      ALLOCATE (actual_x_data%screen_funct_coeffs_kind(1, 1))
      ALLOCATE (actual_x_data%screen_funct_coeffs_pgf(1, 1, 1, 1, 1, 1))
      ALLOCATE (actual_x_data%pair_dist_radii_pgf(1, 1, 1, 1, 1, 1))
      actual_x_data%screen_funct_is_initialized = .TRUE.
      IF (PRESENT(mp2_env)) THEN
         saved_x => qs_env%x_data
         qs_env%x_data => owned_x
         CALL prepare_integral_calc(cell, qs_env, mp2_env, para_env, potential, actual_x_data, &
                                    do_periodic, basis, max_set, particle_set, natom, kind_of, &
                                    nsgf_max, integrals, work, work2, buffer1, buffer2, primitives, &
                                    nspins, max_contraction, max_pgf, pgf_ij, pgf_kl, products, nimages, &
                                    helper_eps, helper_logeps, lib, p_work, screen_set, screen_kind, screen_pgf, radii)
         qs_env%x_data => saved_x
      ELSE
         CALL cd_prepare_molecular(qs_env, actual_x_data, basis, particle_set, kind_of, &
                                   max_set, natom, nsgf_max, max_pgf, integrals, work, work2, &
                                   buffer1, buffer2, primitives, pgf_ij, pgf_kl, products, nimages)
         potential%potential_type = do_potential_coulomb
         lib = actual_x_data%lib
         do_periodic = .FALSE.
         ALLOCATE (max_contraction(0, 0))
      END IF
      ALLOCATE (neutral(max_pgf, max_pgf))
      neutral(:, :)%x(1) = 0.0_dp
      neutral(:, :)%x(2) = 0.0_dp

      ALLOCATE (packed_basis(0:actual_x_data%basis_info%max_am, max_set, SIZE(basis)))
      DO ikind = 1, SIZE(basis)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis, basis_type='ORB')
         DO iset = 1, basis(ikind)%nset
            DO l = basis(ikind)%lmin(iset), basis(ikind)%lmax(iset)
               IF (basis(ikind)%nsgfl(l, iset) == 0) CYCLE
               CALL cd_pack_l_basis(basis(ikind), orb_basis, iset, l, packed_basis(l, iset, ikind))
            END DO
         END DO
      END DO

      nao = 0
      nshell = 0
      DO atom = 1, natom
         ikind = kind_of(atom)
         nao = nao + SUM(basis(ikind)%nsgf)
         DO iset = 1, basis(ikind)%nset
            nshell = nshell + COUNT(basis(ikind)%nsgfl(:, iset) > 0)
         END DO
      END DO
      npair8 = INT(nao, int_8)*(INT(nao, int_8) + 1_int_8)/2_int_8
      IF (npair8 > INT(HUGE(ntotal), int_8) - 1_int_8) THEN
         CPABORT("CD_PRODUCER: pair index overflow")
      END IF
      ntotal = INT(npair8)
      CPASSERT(ntotal > 0)
      ALLOCATE (shells(nshell), pairs(ntotal), d_all(ntotal))
      ishell = 0
      i = 1
      DO atom = 1, natom
         ikind = kind_of(atom)
         DO iset = 1, basis(ikind)%nset
            DO j = 2, basis(ikind)%nshell(iset)
               IF (basis(ikind)%nl(j, iset) < basis(ikind)%nl(j - 1, iset)) THEN
                  CPABORT("CD_PRODUCER: ORB shells must be ordered by angular momentum")
               END IF
            END DO
            ! hfx_create_basis_types packs each l independently in sphi_ext.
            ! Keep the original AO order; the oracle uses private packed primitives.
            DO l = basis(ikind)%lmin(iset), basis(ikind)%lmax(iset)
               IF (basis(ikind)%nsgfl(l, iset) == 0) CYCLE
               ishell = ishell + 1
               shells(ishell) = cd_shell(atom, ikind, iset, l, (2*l + 1)*basis(ikind)%nsgfl(l, iset), i)
               i = i + shells(ishell)%n
            END DO
         END DO
      END DO
      CPASSERT(i == nao + 1)
      ngroup = INT(INT(nshell, int_8)*(INT(nshell, int_8) + 1_int_8)/2_int_8)
      ALLOCATE (bra(ngroup), ket(ngroup), first(ngroup + 1), group_owner(ngroup))
      ALLOCATE (group_weight(ngroup), pair_owner(ntotal))
      g = 0
      n = 0
      DO jshell = 1, nshell
         DO ishell = 1, jshell
            g = g + 1
            bra(g) = ishell
            ket(g) = jshell
            first(g) = n + 1
            DO b = 1, shells(jshell)%n
               DO a = 1, shells(ishell)%n
                  IF (ishell == jshell .AND. a > b) CYCLE
                  n = n + 1
                  pairs(n) = cd_pair(g, a, b, shells(ishell)%first + a - 1, shells(jshell)%first + b - 1)
               END DO
            END DO
            group_weight(g) = INT(n - first(g) + 1, int_8)* &
                              INT(packed_basis(shells(ishell)%l, shells(ishell)%iset, &
                                               shells(ishell)%kind)%npgf_l, int_8)* &
                              INT(packed_basis(shells(jshell)%l, shells(jshell)%iset, &
                                               shells(jshell)%kind)%npgf_l, int_8)
         END DO
      END DO
      first(ngroup + 1) = n + 1
      CPASSERT(n == ntotal)
      CALL assign_groups(group_weight, group_owner)
      d_all = 0.0_dp
      DO g = 1, ngroup
         DO i = first(g), first(g + 1) - 1
            pair_owner(i) = group_owner(g)
         END DO
         IF (group_owner(g) /= para_env%mepos) CYCLE
         ! A diagonal quartet already has the canonical (max(group), min(group)) orientation.
         CALL quartet(g, g)
         ishell = bra(g)
         jshell = ket(g)
         DO i = first(g), first(g + 1) - 1
            a = pairs(i)%a + shells(ishell)%n*(pairs(i)%b - 1)
            j = shells(ishell)%n*shells(jshell)%n
            d_all(i) = integrals(a + j*(a - 1))
         END DO
      END DO
      CALL para_env%sum(d_all)
      IF (.NOT. ALL(ieee_is_finite(d_all))) CPABORT("CD_PRODUCER: nonfinite exact diagonal")
      dmax_initial = MAXVAL(d_all)
      roundoff = 100.0_dp*EPSILON(1.0_dp)*MAX(dmax_initial, TINY(1.0_dp))
      IF (MINVAL(d_all) < -roundoff) CPABORT("CD_PRODUCER: negative exact Coulomb diagonal")
      d_all = MAX(d_all, 0.0_dp)
      ! Stable compaction preserves the packed pair ordering and group contiguity.
      n = 0
      DO g = 1, ngroup
         oldfirst = first(g)
         oldlast = first(g + 1) - 1
         first(g) = n + 1
         DO i = oldfirst, oldlast
            IF (SQRT(d_all(i))*SQRT(MAX(dmax_initial, 0.0_dp)) <= &
                controls%eps_schwarz) CYCLE
            n = n + 1
            pairs(n) = pairs(i)
            d_all(n) = d_all(i)
            pair_owner(n) = pair_owner(i)
         END DO
      END DO
      first(ngroup + 1) = n + 1
      bound = controls%eps_cholesky*dmax_initial
      prune_bound = controls%qualification*bound
      ALLOCATE (active(n), group_active(ngroup))
      active = .TRUE.
      DO g = 1, ngroup
         group_active(g) = first(g) < first(g + 1)
         IF (do_prune .AND. group_active(g)) THEN
            group_active(g) = ANY(d_all(first(g):first(g + 1) - 1) > prune_bound)
         END IF
         IF (.NOT. group_active(g)) active(first(g):first(g + 1) - 1) = .FALSE.
      END DO
      nlocal_full = COUNT(pair_owner(1:n) == para_env%mepos)
      nlocal = COUNT(pair_owner(1:n) == para_env%mepos .AND. active)
      ALLOCATE (d(nlocal), local_index(n))
      local_index = 0
      local_row = 0
      DO i = 1, n
         IF (pair_owner(i) /= para_env%mepos .OR. .NOT. active(i)) CYCLE
         local_row = local_row + 1
         local_index(i) = local_row
         d(local_row) = d_all(i)
      END DO
      CPASSERT(local_row == nlocal)
      DEALLOCATE (d_all)
      t_diag = m_walltime() - t_start

      tick = m_walltime()
      k = 0
      batch = MIN(cd_batch_size, n)
      capacity = batch
      ALLOCATE (factor(nlocal, capacity), scores(n), selected(n), q(n), groups(batch), touched(ngroup))
      max_nlocal = nlocal
      CALL para_env%max(max_nlocal)
      panel_words = panel_bytes/INT(STORAGE_SIZE(0.0_dp)/8, int_8)
      linear_words = INT(max_nlocal, int_8) + INT(n, int_8) + 2_int_8
      high = INT(MIN(INT(n, int_8), panel_words/linear_words))
      low = 0
      DO WHILE (low < high)
         mid = low + (high - low + 1)/2
         IF (INT(mid, int_8) <= panel_words/(linear_words + INT(mid, int_8))) THEN
            low = mid
         ELSE
            high = mid - 1
         END IF
      END DO
      panel_limit = low
      ALLOCATE (cache(n), pivots(n), missing(n), slots(n))
      column_bytes = INT(nlocal, int_8)*INT(STORAGE_SIZE(0.0_dp)/8, int_8)
      selected = .FALSE.
      dmax = 0.0_dp
      IF (nlocal > 0) dmax = MAXVAL(d)
      CALL para_env%max(dmax)
      DO WHILE (dmax > bound .AND. k < n)
         dmax_start = dmax
         ! Disjoint support makes this an exact assembly of the global diagonal.
         scores = 0.0_dp
         DO i = 1, n
            IF (local_index(i) > 0) scores(i) = d(local_index(i))
         END DO
         CALL para_env%sum(scores)
         IF (do_prune) THEN
            ! Preserve complete expanded panels (and their tie ordering): retire a group only
            ! when none of its unselected rows can ever pass qualification again.
            active = .FALSE.
            DO g = 1, ngroup
               group_active(g) = first(g) < first(g + 1)
               IF (group_active(g)) THEN
                  group_active(g) = ANY(.NOT. selected(first(g):first(g + 1) - 1) .AND. &
                                        scores(first(g):first(g + 1) - 1) > prune_bound)
               END IF
               IF (group_active(g)) active(first(g):first(g + 1) - 1) = .TRUE.
            END DO
            CALL compact_active_rows()
         END IF
         WHERE (selected .OR. scores <= controls%qualification*dmax) scores = -1.0_dp
         nq = 0
         DO j = 1, MIN(batch, n - k)
            best_d = MAXVAL(scores)
            IF (best_d < 0.0_dp) EXIT
            ! Ties: relative to the current max, but never narrower than the accumulated Schur-update
            ! roundoff, which scales with k*eps*dmax_initial (not with the current max).
            tie_band = MAX(pivot_tie_tol*best_d, 10.0_dp*EPSILON(1.0_dp)*REAL(MAX(1, k), dp)*dmax_initial)
            ! The replicated maximum makes every rank choose the lowest canonical pair in the tie band.
            best = 0
            DO i = 1, n
               IF (scores(i) >= best_d - tie_band) THEN
                  IF (best == 0) THEN
                     best = i
                  ELSE IF (canonical_pair_index(i) < canonical_pair_index(best)) THEN
                     best = i
                  END IF
               END IF
            END DO
            CPASSERT(best > 0)
            nq = nq + 1
            q(nq) = best
            scores(best) = -1.0_dp
         END DO
         CPASSERT(nq > 0)
         ! Qualification is unchanged. Expand its groups to ALL retained pairs,
         ! including already selected rows; only scalar residual pivots are retained.
         ! Group order follows first occurrence in Q (largest score, stable ties).
         touched = .FALSE.
         ntouched = 0
         DO j = 1, nq
            g = pairs(q(j))%group
            IF (touched(g)) CYCLE
            ntouched = ntouched + 1
            groups(ntouched) = g
            touched(g) = .TRUE.
         END DO
         next_group = 1
         DO WHILE (next_group <= ntouched)
            ! Greedily pack whole groups. Separate chunks repeat cross-chunk
            ! quartets; exchange reuse remains exact within each chunk.
            nq = 0
            DO WHILE (next_group <= ntouched)
               g = groups(next_group)
               group_size = first(g + 1) - first(g)
               IF (group_size > panel_limit) THEN
                  CPABORT("CD_PRODUCER: one candidate group exceeds panel memory cap")
               END IF
               IF (group_size > panel_limit - nq) EXIT
               DO i = first(g), first(g + 1) - 1
                  nq = nq + 1
                  q(nq) = i
               END DO
               next_group = next_group + 1
            END DO
            CPASSERT(nq > 0)
            ALLOCATE (panel(nlocal, nq), pivot_row(nq), qdiag(nq))
            t_la = t_la + m_walltime() - tick
            tick = m_walltime()
            panels = panels + 1_int_8
            nmiss = 0
            DO j = 1, nq
               IF (ALLOCATED(cache(q(j))%values)) THEN
                  columns_cache_hits = columns_cache_hits + 1_int_8
               ELSE
                  nmiss = nmiss + 1
                  missing(nmiss) = q(j)
                  slots(nmiss) = j
               END IF
            END DO
            ! Scatter misses directly into the panel; no second panel-sized buffer.
            CALL columns(missing(1:nmiss), panel, slots(1:nmiss))
            columns_evaluated = columns_evaluated + INT(nmiss, int_8)
            DO j = 1, nq
               i = q(j)
               IF (ALLOCATED(cache(i)%values)) THEN
                  panel(:, j) = cache(i)%values
               ELSE IF (.NOT. cache_full) THEN
                  IF (column_bytes > raw_cache_budget - cache_bytes) THEN
                     cache_full = .TRUE.
                  ELSE
                     ALLOCATE (cache(i)%values(nlocal), STAT=cache_stat)
                     IF (cache_stat /= 0) THEN
                        cache_full = .TRUE.
                     ELSE
                        cache(i)%values(:) = panel(:, j)
                        cache_bytes = cache_bytes + column_bytes
                     END IF
                  END IF
               END IF
            END DO
            t_col = t_col + m_walltime() - tick
            tick = m_walltime()
            qdiag = 0.0_dp
            DO j = 1, nq
               IF (local_index(q(j)) > 0) qdiag(j) = d(local_index(q(j)))
            END DO
            CALL para_env%sum(qdiag)
            IF (k > 0) THEN
               ALLOCATE (qfactor(nq, k))
               qfactor = 0.0_dp
               DO j = 1, k
                  DO i = 1, nq
                     IF (local_index(q(i)) > 0) qfactor(i, j) = factor(local_index(q(i)), j)
                  END DO
               END DO
               CALL para_env%sum(qfactor)
               IF (nlocal > 0) THEN
                  CALL dgemm('N', 'T', nlocal, nq, k, -1.0_dp, factor, nlocal, &
                             qfactor, nq, 1.0_dp, panel, nlocal)
               END IF
               DEALLOCATE (qfactor)
            END IF
            ALLOCATE (qblock(nq, nq))
            qblock = 0.0_dp
            DO i = 1, nq
               IF (local_index(q(i)) > 0) qblock(i, :) = panel(local_index(q(i)), :)
            END DO
            ! Each Q row has exactly one owner, so the assembled Schur block is exact.
            CALL para_env%sum(qblock)
            nleft = MIN(nq, n - k)
            DO m = 1, nleft
               best_d = -1.0_dp
               DO j = 1, nq
                  i = q(j)
                  IF (selected(i)) CYCLE
                  best_d = MAX(best_d, qdiag(j))
               END DO
               IF (best_d < 0.0_dp) EXIT
               tie_band = MAX(pivot_tie_tol*best_d, 10.0_dp*EPSILON(1.0_dp)*REAL(MAX(1, k), dp)*dmax_initial)
               best = 0
               DO j = 1, nq
                  i = q(j)
                  IF (selected(i)) CYCLE
                  IF (qdiag(j) >= best_d - tie_band) THEN
                     IF (best == 0) THEN
                        best = j
                     ELSE IF (canonical_pair_index(i) < canonical_pair_index(q(best))) THEN
                        best = j
                     END IF
                  END IF
               END DO
               IF (best == 0) EXIT
               pivot = q(best)
               ! Delta may terminate the trajectory, never redirect it to a new panel.
               ! A capped qualification set can leave larger residuals outside this panel.
               IF (best_d <= bound) THEN
                  dmax = 0.0_dp
                  IF (nlocal > 0) dmax = MAXVAL(d)
                  CALL para_env%max(dmax)
                  IF (dmax <= bound) EXIT
               END IF
               IF (best_d <= controls%qualification*dmax_start) EXIT
               ! Grow only for selected pivots, not for the wider candidate count.
               IF (k == capacity) THEN
                  newcap = MIN(n, capacity + MAX(capacity/2, batch))
                  ALLOCATE (grown(nlocal, newcap))
                  grown(:, 1:k) = factor(:, 1:k)
                  CALL MOVE_ALLOC(grown, factor)
                  capacity = newcap
               END IF
               k = k + 1
               row_updates_skipped = row_updates_skipped + INT(nlocal_full - nlocal, int_8)
               ! Unscreened canonical AO-pair ID, independent of grouping/compaction.
               pivots(k) = canonical_pair_index(pivot)
               factor(:, k) = panel(:, best)/SQRT(qdiag(best))
               ! qblock is replicated: this is the broadcast pivot row for every rank.
               pivot_row(:) = qblock(:, best)/SQRT(qdiag(best))
               IF (nlocal > 0) THEN
                  CALL dger(nlocal, nq, -1.0_dp, factor(1, k), 1, pivot_row, 1, panel, nlocal)
               END IF
               CALL dger(nq, nq, -1.0_dp, pivot_row, 1, pivot_row, 1, qblock, nq)
               d(:) = d(:) - factor(:, k)**2
               selected(pivot) = .TRUE.
               ! A selected row cannot be chosen again because of roundoff.
               IF (local_index(pivot) > 0) d(local_index(pivot)) = 0.0_dp
               qdiag(:) = qdiag(:) - pivot_row(:)**2
               qdiag(best) = 0.0_dp
               roundoff = 100.0_dp*EPSILON(1.0_dp)*dmax_initial*REAL(MAX(1, k), dp)
               bad = 0
               IF (.NOT. ALL(ieee_is_finite(d))) bad = 1
               CALL para_env%max(bad)
               IF (bad /= 0) CPABORT("CD_PRODUCER: nonfinite residual diagonal")
               bad = 0
               IF (nlocal > 0) THEN
                  IF (MINVAL(d) < -roundoff) bad = 1
               END IF
               CALL para_env%max(bad)
               IF (bad /= 0) CPABORT("CD_PRODUCER: materially negative residual diagonal")
               d = MAX(d, 0.0_dp)
               qdiag = MAX(qdiag, 0.0_dp)
            END DO
            DEALLOCATE (panel, pivot_row, qdiag, qblock)
            dmax = 0.0_dp
            IF (nlocal > 0) dmax = MAXVAL(d)
            CALL para_env%max(dmax)
            IF (dmax <= bound .OR. k == n) EXIT
         END DO
      END DO
      t_la = t_la + m_walltime() - tick
      rel_final = 0.0_dp
      IF (dmax_initial > 0.0_dp) rel_final = dmax/dmax_initial
      rows_active_final = COUNT(active)
      groups_active_final = COUNT(group_active)

      ALLOCATE (lsum(k), lsq(k))
      DO j = 1, k
         lsum(j) = SUM(factor(:, j))
         lsq(j) = SUM(factor(:, j)**2)
      END DO
      CALL para_env%sum(lsum)
      CALL para_env%sum(lsq)

      IF (PRESENT(mu)) THEN
         CPASSERT(PRESENT(nu) .AND. PRESENT(owner))
         ALLOCATE (mu(k), nu(k), owner(k))
         DO j = 1, k
            ! Decode the canonical ID, without floating-point square roots.
            b = 1
            DO WHILE (INT(b, int_8)*INT(b + 1, int_8)/2_int_8 < INT(pivots(j), int_8))
               b = b + 1
            END DO
            mu(j) = pivots(j) - INT(INT(b, int_8)*INT(b - 1, int_8)/2_int_8)
            nu(j) = b
            DO i = 1, nshell
               IF (mu(j) < shells(i)%first .OR. mu(j) >= shells(i)%first + shells(i)%n) CYCLE
               owner(j) = shells(i)%atom
               EXIT
            END DO
            CPASSERT(i <= nshell)
         END DO
      END IF

      DO i = 1, SIZE(pgf_ij)
         DEALLOCATE (pgf_ij(i)%image_list, pgf_kl(i)%image_list)
      END DO
      DEALLOCATE (neutral, pgf_ij, pgf_kl, products, nimages, max_contraction, kind_of, packed_basis)
      DEALLOCATE (integrals, work, work2, buffer1, buffer2, primitives)
      CALL cp_libint_cleanup_eri(actual_x_data%lib)
      CALL hfx_release_basis_types(actual_x_data%basis_parameter)
      DEALLOCATE (actual_x_data%neighbor_cells, actual_x_data%screen_funct_coeffs_set, &
                  actual_x_data%screen_funct_coeffs_kind, actual_x_data%screen_funct_coeffs_pgf, &
                  actual_x_data%pair_dist_radii_pgf)
      DEALLOCATE (owned_x)
      IF (retain) THEN
         IF (k == 0) CPABORT("CD_CONSUMER: empty factor cannot define an auxiliary dimension")
         CALL MOVE_ALLOC(factor, retained_factor)
         retained_rank = k
         retained_nao = nao
         retained_first_column = 1
         ALLOCATE (retained_owner(n))
         retained_owner(:) = pair_owner(1:n)
         ALLOCATE (retained_row(ntotal))
         retained_row = 0
         DO i = 1, n
            j = INT(INT(pairs(i)%nu, int_8)*INT(pairs(i)%nu - 1, int_8)/2_int_8 + INT(pairs(i)%mu, int_8))
            retained_row(j) = i
         END DO
      END IF
      IF (ALLOCATED(factor)) DEALLOCATE (factor)
      DEALLOCATE (scores, selected, active, group_active, q, groups, touched, shells, pairs, d, bra, ket, first)
      DEALLOCATE (group_owner, pair_owner, local_index, group_weight)
      DEALLOCATE (cache, missing, slots)
      t_total = m_walltime() - t_start
      CALL para_env%sum(quartets_unique)
      CALL para_env%sum(integrals_computed)
      CALL para_env%max(columns_evaluated)
      CALL para_env%min(columns_cache_hits)
      CALL para_env%max(cache_bytes)
      CALL para_env%max(panels)
      CALL para_env%sum(row_updates_skipped)
      CALL para_env%sum(quartets_skipped)
      CALL para_env%max(t_diag)
      CALL para_env%max(t_col)
      CALL para_env%max(t_la)
      CALL para_env%max(t_total)
      quartet_time_min = quartet_time
      quartet_time_max = quartet_time
      CALL para_env%min(quartet_time_min)
      CALL para_env%max(quartet_time_max)
      rows_min = nlocal
      rows_max = nlocal
      CALL para_env%min(rows_min)
      CALL para_env%max(rows_max)
      IF (unit_nr > 0 .AND. para_env%is_source()) THEN
         WRITE (unit_nr, '(A,I0)') 'CD_PRODUCER| nao ', nao
         WRITE (unit_nr, '(A,I0)') 'CD_PRODUCER| npairs_total ', ntotal
         WRITE (unit_nr, '(A,I0)') 'CD_PRODUCER| npairs_screened ', n
         WRITE (unit_nr, '(A,ES24.16)') 'CD_PRODUCER| eps_cholesky ', controls%eps_cholesky
         WRITE (unit_nr, '(A,ES24.16)') 'CD_PRODUCER| eps_schwarz ', controls%eps_schwarz
         WRITE (unit_nr, '(A,I0)') 'CD_PRODUCER| batch_size ', cd_batch_size
         WRITE (unit_nr, '(A,ES24.16)') 'CD_PRODUCER| qualification ', controls%qualification
         WRITE (unit_nr, '(A,I0)') 'CD_PRODUCER| k_selected ', k
         WRITE (unit_nr, '(A,ES24.16)') 'CD_PRODUCER| dmax_initial ', dmax_initial
         WRITE (unit_nr, '(A,ES24.16)') 'CD_PRODUCER| dmax_final_rel ', rel_final
         WRITE (unit_nr, '(A,I0)') 'CD_PRODUCER| quartets_unique ', quartets_unique
         WRITE (unit_nr, '(A,I0)') 'CD_PRODUCER| integrals_computed ', integrals_computed
         WRITE (unit_nr, '(A,I0)') 'CD_PRODUCER| columns_evaluated ', columns_evaluated
         WRITE (unit_nr, '(A,F15.6)') 'CD_PRODUCER| t_diag ', t_diag
         WRITE (unit_nr, '(A,F15.6)') 'CD_PRODUCER| t_col ', t_col
         WRITE (unit_nr, '(A,F15.6)') 'CD_PRODUCER| t_la ', t_la
         WRITE (unit_nr, '(A,F15.6)') 'CD_PRODUCER| t_total ', t_total
         WRITE (unit_nr, '(A,F15.6)') 'CD_PRODUCER| t_quartet_rank_min ', quartet_time_min
         WRITE (unit_nr, '(A,F15.6)') 'CD_PRODUCER| t_quartet_rank_max ', quartet_time_max
         WRITE (unit_nr, '(A,I0)') 'CD_PRODUCER| rows_rank_min ', rows_min
         WRITE (unit_nr, '(A,I0)') 'CD_PRODUCER| rows_rank_max ', rows_max
         WRITE (unit_nr, '(A,I0)') 'CD_PRODUCER| panels ', panels
         WRITE (unit_nr, '(A,I0)') 'CD_PRODUCER| columns_cache_hits ', columns_cache_hits
         WRITE (unit_nr, '(A,I0)') 'CD_PRODUCER| cache_bytes ', cache_bytes
         WRITE (unit_nr, '(A,I0)') 'CD_PRODUCER| cache_budget_bytes ', raw_cache_budget
         DO i = 1, k, 20
            WRITE (unit_nr, '(A,20(I0,1X))') 'CD_PIVOTS| ', pivots(i:MIN(i + 19, k))
         END DO
         DO i = 1, MIN(5, k)
            WRITE (unit_nr, '(A,I0,A,ES24.16,A,ES24.16)') &
               'CD_LSUM| column ', i, ' sum ', lsum(i), ' sumsq ', lsq(i)
         END DO
         WRITE (unit_nr, '(A,ES24.16,A,ES24.16)') &
            'CD_LSUM| total sum ', SUM(lsum), ' sumsq ', SUM(lsq)
      END IF
      IF (unit_nr > 0 .AND. para_env%is_source()) THEN
         WRITE (unit_nr, '(A,I0,A,I0,A,I0,A,I0)') &
            'CD_PRODUCER| prune rows_active_final ', rows_active_final, &
            ' groups_active_final ', groups_active_final, &
            ' row_updates_skipped ', row_updates_skipped, ' quartets_skipped ', quartets_skipped
         CALL m_flush(unit_nr)
      END IF
      DEALLOCATE (lsum, lsq)

   CONTAINS

!> \brief Unscreened canonical AO-pair index for one packed pair.
!> \param ip ...
!> \return ...
! **************************************************************************************************
      INTEGER FUNCTION canonical_pair_index(ip) RESULT(index)
      INTEGER, INTENT(IN)                                :: ip

         index = INT(INT(pairs(ip)%nu, int_8)*INT(pairs(ip)%nu - 1, int_8)/2_int_8 &
                     + INT(pairs(ip)%mu, int_8))
      END FUNCTION canonical_pair_index

!> \brief Stably compact local rows after whole shell-pair groups become permanently inactive.
! **************************************************************************************************
      SUBROUTINE compact_active_rows()
      INTEGER                                            :: ip, new_nlocal, new_row, old_nlocal, &
                                                            old_row
      INTEGER, ALLOCATABLE                               :: row_map(:)
      REAL(dp), ALLOCATABLE                              :: compact_column(:), compact_d(:), &
                                                            compact_factor(:, :)

         old_nlocal = nlocal
         new_nlocal = COUNT(pair_owner(1:n) == para_env%mepos .AND. active)
         IF (new_nlocal == old_nlocal) RETURN
         CPASSERT(new_nlocal < old_nlocal)
         ALLOCATE (row_map(old_nlocal), compact_d(new_nlocal), compact_factor(new_nlocal, capacity))
         row_map = 0
         new_row = 0
         DO ip = 1, n
            old_row = local_index(ip)
            IF (old_row == 0) CYCLE
            IF (.NOT. active(ip)) THEN
               local_index(ip) = 0
               CYCLE
            END IF
            new_row = new_row + 1
            row_map(old_row) = new_row
            local_index(ip) = new_row
            compact_d(new_row) = d(old_row)
            IF (k > 0) compact_factor(new_row, 1:k) = factor(old_row, 1:k)
         END DO
         CPASSERT(new_row == new_nlocal)
         CALL MOVE_ALLOC(compact_d, d)
         CALL MOVE_ALLOC(compact_factor, factor)
         column_bytes = INT(new_nlocal, int_8)*INT(STORAGE_SIZE(0.0_dp)/8, int_8)
         cache_bytes = 0_int_8
         DO ip = 1, n
            IF (.NOT. ALLOCATED(cache(ip)%values)) CYCLE
            IF (.NOT. active(ip)) THEN
               DEALLOCATE (cache(ip)%values)
               CYCLE
            END IF
            ALLOCATE (compact_column(new_nlocal))
            DO old_row = 1, old_nlocal
               new_row = row_map(old_row)
               IF (new_row > 0) compact_column(new_row) = cache(ip)%values(old_row)
            END DO
            CALL MOVE_ALLOC(compact_column, cache(ip)%values)
            cache_bytes = cache_bytes + column_bytes
         END DO
         nlocal = new_nlocal
         DEALLOCATE (row_map)
      END SUBROUTINE compact_active_rows

!> \brief Evaluate one shell quartet. Fortran flattening has leg a fastest.
!> \param gbra ...
!> \param gket ...
! **************************************************************************************************
      SUBROUTINE quartet(gbra, gket)
      INTEGER, INTENT(IN)                                :: gbra, gket

      INTEGER                                            :: ia, ib, ic, id, sa, sb, sc, sd
      INTEGER(int_8)                                     :: neris
      REAL(dp)                                           :: estimate
      TYPE(cd_l_basis), POINTER                          :: pa, pb, pc, pd
      TYPE(hfx_basis_type), POINTER                      :: ba, bb, bc, bd

         ia = bra(gbra)
         ib = ket(gbra)
         ic = bra(gket)
         id = ket(gket)
         sa = shells(ia)%iset
         sb = shells(ib)%iset
         sc = shells(ic)%iset
         sd = shells(id)%iset
         ba => basis(shells(ia)%kind)
         bb => basis(shells(ib)%kind)
         bc => basis(shells(ic)%kind)
         bd => basis(shells(id)%kind)
         pa => packed_basis(shells(ia)%l, sa, shells(ia)%kind)
         pb => packed_basis(shells(ib)%l, sb, shells(ib)%kind)
         pc => packed_basis(shells(ic)%l, sc, shells(ic)%kind)
         pd => packed_basis(shells(id)%l, sd, shells(id)%kind)
         ! coulomb4 builds both PGF pair lists from these packed counts/exponents.
         ! The original max_pgf allocations bound every packed list and neutral table.
         quartet_tick = m_walltime()
         CALL coulomb4(lib, particle_set(shells(ia)%atom)%r, particle_set(shells(ib)%atom)%r, &
                       particle_set(shells(ic)%atom)%r, particle_set(shells(id)%atom)%r, &
                       pa%npgf_l, pb%npgf_l, pc%npgf_l, pd%npgf_l, &
                       shells(ia)%l, shells(ia)%l, shells(ib)%l, shells(ib)%l, &
                       shells(ic)%l, shells(ic)%l, shells(id)%l, shells(id)%l, &
                       shells(ia)%n, shells(ib)%n, shells(ic)%n, shells(id)%n, &
                       SIZE(pa%sphi_ext, 1), SIZE(pa%sphi_ext, 2), SIZE(pa%sphi_ext, 3), &
                       SIZE(pb%sphi_ext, 1), SIZE(pb%sphi_ext, 2), SIZE(pb%sphi_ext, 3), &
                       SIZE(pc%sphi_ext, 1), SIZE(pc%sphi_ext, 2), SIZE(pc%sphi_ext, 3), &
                       SIZE(pd%sphi_ext, 1), SIZE(pd%sphi_ext, 2), SIZE(pd%sphi_ext, 3), &
                       pa%zet_l, pb%zet_l, pc%zet_l, pd%zet_l, integrals, potential, &
                       actual_x_data%neighbor_cells, [0.0_dp, 0.0_dp], [0.0_dp, 0.0_dp], 0.0_dp, &
                       1.0_dp, estimate, cell, neris, 0.0_dp, log_zero, neutral, neutral, neutral, neutral, &
                       pgf_ij, pgf_kl, products, ba%nsgfl(:, sa), bb%nsgfl(:, sb), &
                       bc%nsgfl(:, sc), bd%nsgfl(:, sd), pa%sphi_ext, pb%sphi_ext, pc%sphi_ext, pd%sphi_ext, &
                       work, work2, buffer1, buffer2, primitives, nimages, do_periodic, p_work)
         quartet_time = quartet_time + m_walltime() - quartet_tick
         quartets_unique = quartets_unique + 1_int_8
         integrals_computed = integrals_computed + INT(shells(ia)%n, int_8)*INT(shells(ib)%n, int_8)* &
                              INT(shells(ic)%n, int_8)*INT(shells(id)%n, int_8)
      END SUBROUTINE quartet

!> \brief Evaluate owned bra rows in canonical group orientation, with local exchange reuse.
!> \param indices ...
!> \param values ...
!> \param destinations ...
! **************************************************************************************************
      SUBROUTINE columns(indices, values, destinations)
      INTEGER, INTENT(IN)                                :: indices(:)
      REAL(dp), INTENT(OUT)                              :: values(:, :)
      INTEGER, INTENT(IN)                                :: destinations(:)

      INTEGER                                            :: ab, cd, col, gbra, h, ip, na, nab, nc, &
                                                            row
      INTEGER, ALLOCATABLE                               :: head(:), next_col(:)
      LOGICAL                                            :: swapped

         ALLOCATE (head(ngroup), next_col(SIZE(indices)))
         head = 0
         DO col = SIZE(indices), 1, -1
            h = pairs(indices(col))%group
            next_col(col) = head(h)
            head(h) = col
         END DO
         values = 0.0_dp
         DO h = 1, ngroup
            IF (head(h) == 0) CYCLE
            DO gbra = 1, ngroup
               IF (first(gbra) == first(gbra + 1)) CYCLE
               IF (group_owner(gbra) /= para_env%mepos) CYCLE
               IF (.NOT. group_active(gbra)) THEN
                  quartets_skipped = quartets_skipped + 1_int_8
                  CYCLE
               END IF
               ! Reuse exchange symmetry only when this rank owns both row groups.
               IF (head(gbra) /= 0 .AND. gbra < h .AND. &
                   group_owner(h) == para_env%mepos) CYCLE
               ! The quartet arithmetic depends only on its groups, never on batch composition,
               ! cache state, or ownership. The larger group is always the bra group.
               swapped = gbra < h
               IF (swapped) THEN
                  CALL quartet(h, gbra)
                  na = shells(bra(h))%n
                  nc = shells(bra(gbra))%n
                  nab = na*shells(ket(h))%n
               ELSE
                  CALL quartet(gbra, h)
                  na = shells(bra(gbra))%n
                  nc = shells(bra(h))%n
                  nab = na*shells(ket(gbra))%n
               END IF
               col = head(h)
               DO WHILE (col /= 0)
                  ip = indices(col)
                  IF (swapped) THEN
                     ab = pairs(ip)%a + na*(pairs(ip)%b - 1)
                     DO row = first(gbra), first(gbra + 1) - 1
                        cd = pairs(row)%a + nc*(pairs(row)%b - 1)
                        values(local_index(row), destinations(col)) = integrals(ab + nab*(cd - 1))
                     END DO
                  ELSE
                     cd = pairs(ip)%a + nc*(pairs(ip)%b - 1)
                     DO row = first(gbra), first(gbra + 1) - 1
                        ab = pairs(row)%a + na*(pairs(row)%b - 1)
                        values(local_index(row), destinations(col)) = integrals(ab + nab*(cd - 1))
                     END DO
                  END IF
                  col = next_col(col)
               END DO
               IF (gbra == h .OR. group_owner(h) /= para_env%mepos) CYCLE
               ! (ab|cd) = (cd|ab): scatter the transposed block without another call.
               col = head(gbra)
               DO WHILE (col /= 0)
                  ip = indices(col)
                  ab = pairs(ip)%a + na*(pairs(ip)%b - 1)
                  DO row = first(h), first(h + 1) - 1
                     cd = pairs(row)%a + nc*(pairs(row)%b - 1)
                     values(local_index(row), destinations(col)) = integrals(ab + nab*(cd - 1))
                  END DO
                  col = next_col(col)
               END DO
            END DO
         END DO
      END SUBROUTINE columns

!> \brief Deterministic longest-processing-time assignment of whole shell-pair groups.
!> \param weight ...
!> \param owner ...
! **************************************************************************************************
      SUBROUTINE assign_groups(weight, owner)
      INTEGER(int_8), INTENT(IN)                         :: weight(:)
      INTEGER, INTENT(OUT)                               :: owner(:)

      INTEGER                                            :: ia, ib, idx, left, middle, out, rank, &
                                                            right, width
      INTEGER(int_8), ALLOCATABLE                        :: load(:)
      INTEGER, ALLOCATABLE                               :: order(:), work_order(:)

         ALLOCATE (order(SIZE(weight)), work_order(SIZE(weight)))
         order(:) = [(idx, idx=1, SIZE(weight))]
         width = 1
         DO WHILE (width < SIZE(weight))
            left = 1
            DO WHILE (left <= SIZE(weight))
               middle = MIN(left + width, SIZE(weight) + 1)
               right = MIN(left + 2*width, SIZE(weight) + 1)
               ia = left
               ib = middle
               DO out = left, right - 1
                  IF (ia < middle .AND. ib < right) THEN
                     IF (weight(order(ia)) > weight(order(ib)) .OR. &
                         (weight(order(ia)) == weight(order(ib)) .AND. order(ia) < order(ib))) THEN
                        work_order(out) = order(ia)
                        ia = ia + 1
                     ELSE
                        work_order(out) = order(ib)
                        ib = ib + 1
                     END IF
                  ELSE IF (ia < middle) THEN
                     work_order(out) = order(ia)
                     ia = ia + 1
                  ELSE
                     work_order(out) = order(ib)
                     ib = ib + 1
                  END IF
               END DO
               left = left + 2*width
            END DO
            order(:) = work_order(:)
            width = 2*width
         END DO
         ALLOCATE (load(para_env%num_pe))
         load = 0_int_8
         DO idx = 1, SIZE(order)
            rank = MINLOC(load, DIM=1)
            owner(order(idx)) = rank - 1
            load(rank) = load(rank) + weight(order(idx))
         END DO
         DEALLOCATE (order, work_order, load)
      END SUBROUTINE assign_groups
   END SUBROUTINE cd_build_once
!> \brief Stream raw four-centre product columns once, with one molecular image.
!> \param qs_env ...
!> \param para_env ...
!> \param mu ...
!> \param nu ...
!> \param cutoff ...
!> \param eps_filter ...
!> \param matrix_only Evaluate only pivot shell groups, for a product metric.
!> \param consume Called collectively with a rank-partial, unscaled AO-pair column, AO shell indices,
!>                and the global source rank of each shell-pair group (-1 for omitted groups).
!>                Consumers reduce pivot entries or exchange only tensor-owner payloads.
!> \note Exact Gaussian support is infinite. Never use zero HFX product radii here:
!>       those would truncate the product centres, in addition to the physical operator.
! **************************************************************************************************
   SUBROUTINE cd_product_columns(qs_env, para_env, mu, nu, cutoff, eps_filter, matrix_only, consume)
      TYPE(qs_environment_type), POINTER :: qs_env
      TYPE(mp_para_env_type), INTENT(IN) :: para_env
      INTEGER, INTENT(IN) :: mu(:), nu(:)
      REAL(dp), INTENT(IN) :: cutoff, eps_filter
      LOGICAL, INTENT(IN) :: matrix_only
      PROCEDURE(cd_column_consumer) :: consume
      TYPE(cell_type), POINTER :: cell
      TYPE(qs_kind_type), POINTER :: qs_kind_set(:)
      TYPE(particle_type), POINTER :: particle_set(:)
      TYPE(hfx_type), POINTER :: owned_x(:, :), actual_x_data
      TYPE(hfx_basis_type), POINTER :: basis(:)
      TYPE(gto_basis_set_type), POINTER :: orb_basis
      TYPE(cd_l_basis), ALLOCATABLE, TARGET :: packed_basis(:, :, :)
      TYPE(hfx_potential_type) :: potential
      TYPE(cp_libint_t) :: lib
      TYPE(hfx_pgf_list), ALLOCATABLE :: pgf_ij(:), pgf_kl(:)
      TYPE(hfx_pgf_product_list), ALLOCATABLE :: products(:)
      TYPE(hfx_screen_coeff_type), POINTER :: neutral(:, :), support(:, :)
      REAL(dp), POINTER :: p_work(:)
      ! Column-slab budget of the product-column batches (same 2 GiB policy as the selector panel/cache).
      INTEGER(int_8), PARAMETER :: column_slab_budget = 2_int_8*1024_int_8**3
      REAL(dp), ALLOCATABLE :: integrals(:), work(:), work2(:), buffer1(:), buffer2(:), primitives(:), &
                               diagonal(:), columns(:, :, :)
      TYPE(cd_shell), ALLOCATABLE :: shells(:)
      TYPE(cd_pair), ALLOCATABLE :: pairs(:)
      INTEGER, ALLOCATABLE :: kind_of(:), nimages(:), bra(:), ket(:), first(:), pivot_pair(:), &
                              count_g(:), start_g(:), order(:), ao_shell(:), group_owner(:)
      LOGICAL, ALLOCATABLE :: needed(:)
      INTEGER :: max_set, natom, nsgf_max, max_pgf, periodic(3), nao, nshell, ngroup, ntotal, n, &
                 atom, ikind, iset, ishell, jshell, g, a, b, i, j, l, ip, pair, gp, offset, &
                 handle, unit_id, n_owned, handle2, ibatch, b0, nb, nb_max, t
      INTEGER(int_8) :: npair8, n_done, n_surv, h_done, h_surv, hq
      REAL(dp) :: tolerance, margin
      LOGICAL :: do_periodic

      CALL timeset('cd_product_columns', handle)
      CPASSERT(SIZE(mu) == SIZE(nu) .AND. SIZE(mu) > 0)
      CPASSERT(ALL(mu <= nu) .AND. ALL(mu > 0))
      IF (.NOT. ieee_is_finite(cutoff) .OR. .NOT. ieee_is_finite(eps_filter)) THEN
         CPABORT("CD_LCG: nonfinite integral cutoff or filter")
      END IF
      NULLIFY (p_work)
      CALL get_qs_env(qs_env, cell=cell, qs_kind_set=qs_kind_set)
      CALL get_cell(cell, periodic=periodic)
      CPASSERT(ALL(periodic == 0))
      IF (omp_get_max_threads() /= 1) CPABORT("CD_LCG: requires one OpenMP thread")
      ALLOCATE (owned_x(1, 1))
      actual_x_data => owned_x(1, 1)
      CALL hfx_create_basis_types(actual_x_data%basis_parameter, actual_x_data%basis_info, qs_kind_set, 'ORB')
      IF (actual_x_data%basis_info%max_am > libint_max_am) THEN
         CPABORT("CD_PRODUCER: ORB angular momentum exceeds the libint build limit")
      END IF
      CALL cp_libint_static_init()
      CALL cp_libint_init_eri(actual_x_data%lib, actual_x_data%basis_info%max_am)
      CALL cp_libint_set_contrdepth(actual_x_data%lib, 1)
      ALLOCATE (actual_x_data%neighbor_cells(1))
      actual_x_data%periodic_parameter%do_periodic = .FALSE.
      actual_x_data%screening_parameter%eps_schwarz = 0.0_dp
      CALL cd_prepare_molecular(qs_env, actual_x_data, basis, particle_set, kind_of, &
                                max_set, natom, nsgf_max, max_pgf, integrals, work, work2, &
                                buffer1, buffer2, primitives, pgf_ij, pgf_kl, products, nimages)
      lib = actual_x_data%lib
      do_periodic = .FALSE.
      potential%potential_type = do_potential_coulomb
      IF (cutoff > 0.0_dp) THEN
         potential%potential_type = do_potential_truncated
         potential%cutoff_radius = cutoff
         IF (get_lmax_init() < 4*actual_x_data%basis_info%max_am) THEN
            unit_id = -1
            IF (para_env%is_source()) CALL open_file(unit_number=unit_id, file_name='t_c_g.dat')
            CALL init(4*actual_x_data%basis_info%max_am, unit_id, para_env%mepos, para_env)
            IF (para_env%is_source()) CALL close_file(unit_id)
            init_t_c_g0_lmax = 4*actual_x_data%basis_info%max_am
         END IF
      END IF
      ALLOCATE (neutral(max_pgf, max_pgf), support(max_pgf, max_pgf))
      neutral(:, :)%x(1) = 0.0_dp
      neutral(:, :)%x(2) = 0.0_dp
      support(:, :)%x(1) = 0.0_dp
      support(:, :)%x(2) = ieee_value(0.0_dp, ieee_positive_inf)
      ALLOCATE (packed_basis(0:actual_x_data%basis_info%max_am, max_set, SIZE(basis)))
      DO ikind = 1, SIZE(basis)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis, basis_type='ORB')
         DO iset = 1, basis(ikind)%nset
            DO l = basis(ikind)%lmin(iset), basis(ikind)%lmax(iset)
               IF (basis(ikind)%nsgfl(l, iset) == 0) CYCLE
               CALL cd_pack_l_basis(basis(ikind), orb_basis, iset, l, packed_basis(l, iset, ikind))
            END DO
         END DO
      END DO

      nao = 0
      nshell = 0
      DO atom = 1, natom
         ikind = kind_of(atom)
         nao = nao + SUM(basis(ikind)%nsgf)
         DO iset = 1, basis(ikind)%nset
            nshell = nshell + COUNT(basis(ikind)%nsgfl(:, iset) > 0)
         END DO
      END DO
      npair8 = INT(nao, int_8)*(INT(nao, int_8) + 1_int_8)/2_int_8
      IF (npair8 > INT(HUGE(ntotal), int_8) - 1_int_8) THEN
         CPABORT("CD_PRODUCER: pair index overflow")
      END IF
      ntotal = INT(npair8)
      CPASSERT(ntotal > 0)
      ALLOCATE (shells(nshell), pairs(ntotal), ao_shell(nao))
      ishell = 0
      i = 1
      DO atom = 1, natom
         ikind = kind_of(atom)
         DO iset = 1, basis(ikind)%nset
            DO j = 2, basis(ikind)%nshell(iset)
               IF (basis(ikind)%nl(j, iset) < basis(ikind)%nl(j - 1, iset)) THEN
                  CPABORT("CD_PRODUCER: ORB shells must be ordered by angular momentum")
               END IF
            END DO
            ! hfx_create_basis_types packs each l independently in sphi_ext.
            ! Keep the original AO order; the oracle uses private packed primitives.
            DO l = basis(ikind)%lmin(iset), basis(ikind)%lmax(iset)
               IF (basis(ikind)%nsgfl(l, iset) == 0) CYCLE
               ishell = ishell + 1
               shells(ishell) = cd_shell(atom, ikind, iset, l, (2*l + 1)*basis(ikind)%nsgfl(l, iset), i)
               ao_shell(i:i + shells(ishell)%n - 1) = ishell
               i = i + shells(ishell)%n
            END DO
         END DO
      END DO
      CPASSERT(i == nao + 1)
      ngroup = INT(INT(nshell, int_8)*(INT(nshell, int_8) + 1_int_8)/2_int_8)
      ALLOCATE (bra(ngroup), ket(ngroup), first(ngroup + 1))
      g = 0
      n = 0
      DO jshell = 1, nshell
         DO ishell = 1, jshell
            g = g + 1
            bra(g) = ishell
            ket(g) = jshell
            first(g) = n + 1
            DO b = 1, shells(jshell)%n
               DO a = 1, shells(ishell)%n
                  IF (ishell == jshell .AND. a > b) CYCLE
                  n = n + 1
                  pairs(n) = cd_pair(g, a, b, shells(ishell)%first + a - 1, shells(jshell)%first + b - 1)
               END DO
            END DO
         END DO
      END DO
      first(ngroup + 1) = n + 1
      CPASSERT(n == ntotal)
      ALLOCATE (pivot_pair(SIZE(mu)), needed(ngroup), diagonal(ngroup))
      pivot_pair = 0
      needed = .FALSE.
      DO ip = 1, SIZE(mu)
         DO i = 1, ntotal
            IF (pairs(i)%mu /= mu(ip) .OR. pairs(i)%nu /= nu(ip)) CYCLE
            pivot_pair(ip) = i
            needed(pairs(i)%group) = .TRUE.
            EXIT
         END DO
         CPASSERT(pivot_pair(ip) > 0)
      END DO
      ! Both full and spherical-cutoff Coulomb are positive semidefinite.
      ! Bound every entry of a shell-pair block by sqrt(d_g d_h), using the
      ! maximum EXACT contracted diagonal for this operator, not RI_AUX radii.
      diagonal = 0.0_dp
      tolerance = 0.0_dp
      IF (.NOT. matrix_only) THEN
         tolerance = eps_filter/SQRT(REAL(SIZE(mu), dp))/REAL(nao, dp)
         DO g = 1, ngroup
            IF (MOD(g - 1, para_env%num_pe) /= para_env%mepos) CYCLE
            CALL quartet(g, g)
            j = shells(bra(g))%n*shells(ket(g))%n
            DO i = first(g), first(g + 1) - 1
               a = pairs(i)%a + shells(bra(g))%n*(pairs(i)%b - 1)
               diagonal(g) = MAX(diagonal(g), ABS(integrals(a + j*(a - 1))))
            END DO
         END DO
         CALL para_env%sum(diagonal)
         IF (.NOT. ALL(ieee_is_finite(diagonal))) CPABORT("CD_LCG: nonfinite Schwarz diagonal")
         margin = 100.0_dp*EPSILON(1.0_dp)*MAX(MAXVAL(diagonal), TINY(1.0_dp))
         diagonal = diagonal + margin
      END IF
      ! Batch the pivots by ket shell-pair group gp: every surviving (g, gp) quartet is evaluated once per
      ! batch and all selected ket components are scattered from it. Raw columns are unchanged; only the
      ! consume order differs, and both consumers are order independent.
      ALLOCATE (count_g(ngroup), start_g(ngroup + 1), order(SIZE(mu)))
      count_g = 0
      DO ip = 1, SIZE(mu)
         gp = pairs(pivot_pair(ip))%group
         count_g(gp) = count_g(gp) + 1
      END DO
      start_g(1) = 1
      DO g = 1, ngroup
         start_g(g + 1) = start_g(g) + count_g(g)
      END DO
      count_g = 0
      DO ip = 1, SIZE(mu)
         gp = pairs(pivot_pair(ip))%group
         order(start_g(gp) + count_g(gp)) = ip
         count_g(gp) = count_g(gp) + 1
      END DO
      ! Size the column slab so that a whole pivot shell-pair group fits in one batch whenever the
      ! 2 GiB slab budget allows: a group split over several batches recomputes its shell quartets
      ! once per batch (measured 1.71x redundant integral work at 1392 AOs with the former 128 MB cap).
      nb_max = MAX(1, INT(MIN(REAL(MAXVAL(count_g), dp), &
                              REAL(column_slab_budget, dp)/(8.0_dp*REAL(nao, dp)**2))))
      DEALLOCATE (count_g)
      ALLOCATE (columns(nao, nao, nb_max), group_owner(ngroup))
      n_done = 0
      n_surv = 0
      h_done = 0
      h_surv = 0
      ibatch = 0
      DO gp = 1, ngroup
         DO b0 = start_g(gp), start_g(gp + 1) - 1, nb_max
            nb = MIN(nb_max, start_g(gp + 1) - b0)
            ibatch = ibatch + 1
            columns(:, :, 1:nb) = 0.0_dp
            group_owner = -1
            CALL timeset("cd_raw_quartets", handle2)
            ! Surviving bra groups are distributed round robin over ranks, rotated by the batch index.
            n_owned = ibatch - 1
            DO g = 1, ngroup
               IF (matrix_only .AND. .NOT. needed(g)) CYCLE
               IF (.NOT. matrix_only) THEN
                  IF (SQRT(diagonal(g))*SQRT(diagonal(gp)) < tolerance) CYCLE
               END IF
               n_owned = n_owned + 1
               ! Hash of this (bra group, batch) quartet: owners sum it, every rank sums all survivors.
               hq = MOD(INT(g, int_8)*2654435761_int_8 + INT(ibatch, int_8)*40503_int_8, 2147483647_int_8)
               h_surv = h_surv + hq
               group_owner(g) = MOD(n_owned - 1, para_env%num_pe)
               IF (group_owner(g) /= para_env%mepos) CYCLE
               n_done = n_done + 1
               h_done = h_done + hq
               CALL quartet(g, gp)
               j = shells(bra(g))%n*shells(ket(g))%n
               DO t = 1, nb
                  pair = pivot_pair(order(b0 + t - 1))
                  offset = j*(pairs(pair)%a - 1 + shells(bra(gp))%n*(pairs(pair)%b - 1))
                  DO i = first(g), first(g + 1) - 1
                     a = pairs(i)%a + shells(bra(g))%n*(pairs(i)%b - 1)
                     columns(pairs(i)%mu, pairs(i)%nu, t) = integrals(offset + a)
                     columns(pairs(i)%nu, pairs(i)%mu, t) = integrals(offset + a)
                  END DO
               END DO
            END DO
            n_surv = n_surv + INT(n_owned - (ibatch - 1), int_8)
            CALL timestop(handle2)
            IF (.NOT. ALL(ieee_is_finite(columns(:, :, 1:nb)))) THEN
               CPABORT("CD_LCG: nonfinite raw product integral")
            END IF
            CALL timeset("cd_raw_consume", handle2)
            DO t = 1, nb
               CALL consume(order(b0 + t - 1), columns(:, :, t), ao_shell, group_owner)
            END DO
            CALL timestop(handle2)
         END DO
      END DO
      ! Every surviving group must have been evaluated by exactly one rank.
      CALL para_env%sum(n_done)
      CALL para_env%sum(h_done)
      IF (n_done /= n_surv) THEN
         CPABORT("CD_LCG: raw-column group ownership mismatch between ranks")
      END IF
      IF (h_done /= h_surv) CPABORT("CD_LCG: raw-column quartet ownership hash mismatch")
      DO i = 1, SIZE(pgf_ij)
         DEALLOCATE (pgf_ij(i)%image_list, pgf_kl(i)%image_list)
      END DO
      DEALLOCATE (neutral, support)
      CALL cp_libint_cleanup_eri(actual_x_data%lib)
      CALL hfx_release_basis_types(actual_x_data%basis_parameter)
      DEALLOCATE (actual_x_data%neighbor_cells)
      DEALLOCATE (owned_x)
      CALL timestop(handle)

   CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param gbra ...
!> \param gket ...
! **************************************************************************************************
      SUBROUTINE quartet(gbra, gket)
      INTEGER, INTENT(IN)                                :: gbra, gket

      INTEGER                                            :: ia, ib, ic, id, sa, sb, sc, sd
      INTEGER(int_8)                                     :: neris
      REAL(dp)                                           :: estimate
      TYPE(cd_l_basis), POINTER                          :: pa, pb, pc, pd
      TYPE(hfx_basis_type), POINTER                      :: ba, bb, bc, bd

         ia = bra(gbra)
         ib = ket(gbra)
         ic = bra(gket)
         id = ket(gket)
         sa = shells(ia)%iset
         sb = shells(ib)%iset
         sc = shells(ic)%iset
         sd = shells(id)%iset
         ba => basis(shells(ia)%kind)
         bb => basis(shells(ib)%kind)
         bc => basis(shells(ic)%kind)
         bd => basis(shells(id)%kind)
         pa => packed_basis(shells(ia)%l, sa, shells(ia)%kind)
         pb => packed_basis(shells(ib)%l, sb, shells(ib)%kind)
         pc => packed_basis(shells(ic)%l, sc, shells(ic)%kind)
         pd => packed_basis(shells(id)%l, sd, shells(id)%kind)
         ! coulomb4 builds both PGF pair lists from these packed counts/exponents.
         ! The original max_pgf allocations bound every packed list and neutral table.
         CALL coulomb4(lib, particle_set(shells(ia)%atom)%r, particle_set(shells(ib)%atom)%r, &
                       particle_set(shells(ic)%atom)%r, particle_set(shells(id)%atom)%r, &
                       pa%npgf_l, pb%npgf_l, pc%npgf_l, pd%npgf_l, &
                       shells(ia)%l, shells(ia)%l, shells(ib)%l, shells(ib)%l, &
                       shells(ic)%l, shells(ic)%l, shells(id)%l, shells(id)%l, &
                       shells(ia)%n, shells(ib)%n, shells(ic)%n, shells(id)%n, &
                       SIZE(pa%sphi_ext, 1), SIZE(pa%sphi_ext, 2), SIZE(pa%sphi_ext, 3), &
                       SIZE(pb%sphi_ext, 1), SIZE(pb%sphi_ext, 2), SIZE(pb%sphi_ext, 3), &
                       SIZE(pc%sphi_ext, 1), SIZE(pc%sphi_ext, 2), SIZE(pc%sphi_ext, 3), &
                       SIZE(pd%sphi_ext, 1), SIZE(pd%sphi_ext, 2), SIZE(pd%sphi_ext, 3), &
                       pa%zet_l, pb%zet_l, pc%zet_l, pd%zet_l, integrals, potential, &
                       actual_x_data%neighbor_cells, [0.0_dp, 0.0_dp], [0.0_dp, 0.0_dp], 0.0_dp, &
                       1.0_dp, estimate, cell, neris, 0.0_dp, log_zero, support, support, neutral, neutral, &
                       pgf_ij, pgf_kl, products, ba%nsgfl(:, sa), bb%nsgfl(:, sb), &
                       bc%nsgfl(:, sc), bd%nsgfl(:, sd), pa%sphi_ext, pb%sphi_ext, pc%sphi_ext, pd%sphi_ext, &
                       work, work2, buffer1, buffer2, primitives, nimages, do_periodic, p_work)
      END SUBROUTINE quartet
   END SUBROUTINE cd_product_columns

END MODULE mp2_cd_producer
