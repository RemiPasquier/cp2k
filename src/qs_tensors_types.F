!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility methods to build 3-center integral tensors of various types.
! **************************************************************************************************

MODULE qs_tensors_types
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type
   USE cp_array_utils,                  ONLY: cp_1d_i_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE dbt_api,                         ONLY: dbt_create,&
                                              dbt_default_distvec,&
                                              dbt_distribution_destroy,&
                                              dbt_distribution_new,&
                                              dbt_distribution_type,&
                                              dbt_mp_environ_pgrid,&
                                              dbt_pgrid_type,&
                                              dbt_type
   USE distribution_2d_types,           ONLY: distribution_2d_create_prv => distribution_2d_create,&
                                              distribution_2d_release,&
                                              distribution_2d_type
   USE message_passing,                 ONLY: mp_cart_type,&
                                              mp_comm_type,&
                                              mp_para_env_release,&
                                              mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_iterator_p_type,&
                                              neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tensors_types'

   PUBLIC :: distribution_3d_type, neighbor_list_3c_type, neighbor_list_3c_iterator_type, &
             distribution_2d_create, distribution_3d_create, distribution_3d_destroy, &
             split_block_sizes, create_3c_tensor, create_2c_tensor, contiguous_tensor_dist, pgf_block_sizes, &
             create_tensor_batches

   INTEGER, PARAMETER, PUBLIC :: symmetric_none = 0, symmetric_ij = 1, symmetric_jk = 2, symmetrik_ik = 3, symmetric_ijk = 4

   INTEGER, PARAMETER, PUBLIC :: default_block_size = 64
   !! default block size for dense tensors, this block size should be covered by DBCSR/libcusmm

   TYPE distribution_3d_type
      TYPE(distribution_2d_type), POINTER :: dist_2d_1 => NULL(), dist_2d_2 => NULL()
      TYPE(mp_comm_type) :: comm_3d = mp_comm_type(), comm_2d_1 = mp_comm_type(), comm_2d_2 = mp_comm_type()
      LOGICAL :: owns_comm = .FALSE.
   END TYPE distribution_3d_type

   TYPE neighbor_list_3c_type
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: ij_list => NULL(), jk_list => NULL()
      INTEGER :: sym = symmetric_none
      TYPE(distribution_3d_type) :: dist_3d = distribution_3d_type()
      LOGICAL :: owns_dist = .FALSE.
   END TYPE

   TYPE neighbor_list_3c_iterator_type
      TYPE(neighbor_list_iterator_p_type), DIMENSION(:), POINTER :: iter_ij => NULL()
      TYPE(neighbor_list_iterator_p_type), DIMENSION(:), POINTER :: iter_jk => NULL()
      INTEGER                                                    :: iter_level = 0
      TYPE(neighbor_list_3c_type)                                :: ijk_nl = neighbor_list_3c_type()
      INTEGER, DIMENSION(2)                                      :: bounds_i = 0, bounds_j = 0, bounds_k = 0
   END TYPE

CONTAINS
! **************************************************************************************************
!> \brief Create a 3d distribution
!> \param dist_3d 3d distribution object
!> \param dist1 distribution vector along 1st process grid dimension
!> \param dist2 distribution vector along 2nd process grid dimension
!> \param dist3 distribution vector along 3rd process grid dimension
!> \param nkind ...
!> \param particle_set ...
!> \param mp_comm_3d MPI communicator with a 3d cartesian topology
!> \param own_comm Whether mp_comm_3d should be owned by dist_3d (default false)
! **************************************************************************************************
   SUBROUTINE distribution_3d_create(dist_3d, dist1, dist2, dist3, nkind, particle_set, mp_comm_3d, own_comm)
      TYPE(distribution_3d_type)                         :: dist_3d
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dist1, dist2, dist3
      INTEGER, INTENT(IN)                                :: nkind
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(mp_cart_type), INTENT(IN)                     :: mp_comm_3d
      LOGICAL, INTENT(IN), OPTIONAL                      :: own_comm

      CHARACTER(len=*), PARAMETER :: routineN = 'distribution_3d_create'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(2)                              :: mp_coor_1, mp_coor_2
      TYPE(mp_cart_type)                                 :: comm_2d_1, comm_2d_2

      CALL timeset(routineN, handle)

      IF (PRESENT(own_comm)) THEN
         IF (own_comm) dist_3d%comm_3d = mp_comm_3d
         dist_3d%owns_comm = own_comm
      ELSE
         dist_3d%owns_comm = .FALSE.
      END IF

      CALL comm_2d_1%from_sub(mp_comm_3d, [.TRUE., .TRUE., .FALSE.])
      CALL comm_2d_2%from_sub(mp_comm_3d, [.FALSE., .TRUE., .TRUE.])

      mp_coor_1 = comm_2d_1%mepos_cart
      mp_coor_2 = comm_2d_2%mepos_cart

      CPASSERT(mp_coor_1(2) == mp_coor_2(1))

      CALL distribution_2d_create(dist_3d%dist_2d_1, dist1, dist2, nkind, particle_set, comm_2d_1)
      CALL distribution_2d_create(dist_3d%dist_2d_2, dist2, dist3, nkind, particle_set, comm_2d_2)

      dist_3d%comm_2d_1 = comm_2d_1
      dist_3d%comm_2d_2 = comm_2d_2

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Destroy a 3d distribution
!> \param dist ...
! **************************************************************************************************
   SUBROUTINE distribution_3d_destroy(dist)
      TYPE(distribution_3d_type)                         :: dist

      CHARACTER(len=*), PARAMETER :: routineN = 'distribution_3d_destroy'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      CALL distribution_2d_release(dist%dist_2d_1)
      CALL distribution_2d_release(dist%dist_2d_2)
      CALL dist%comm_2d_1%free()
      CALL dist%comm_2d_2%free()
      IF (dist%owns_comm) CALL dist%comm_3d%free()

      NULLIFY (dist%dist_2d_1, dist%dist_2d_2)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Create a 2d distribution. This mainly wraps distribution_2d_create
!>        for consistency with distribution_3d_create.
!> \param dist_2d 2d distribution object
!> \param dist1 distribution vector along 1st process grid dimension
!> \param dist2 distribution vector along 2nd process grid dimension
!> \param nkind ...
!> \param particle_set ...
!> \param mp_comm_2d MPI communicator with a 3d cartesian topology
!> \param blacs_env_ext ...
! **************************************************************************************************
   SUBROUTINE distribution_2d_create(dist_2d, dist1, dist2, nkind, particle_set, mp_comm_2d, blacs_env_ext)
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dist1, dist2
      INTEGER, INTENT(IN)                                :: nkind
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(mp_cart_type), INTENT(IN), OPTIONAL           :: mp_comm_2d
      TYPE(cp_blacs_env_type), OPTIONAL, POINTER         :: blacs_env_ext

      CHARACTER(len=*), PARAMETER :: routineN = 'distribution_2d_create'

      INTEGER                                            :: handle, iatom, ikind, n, natom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nparticle_local_col, nparticle_local_row
      INTEGER, DIMENSION(2)                              :: mp_coor, mp_dims
      INTEGER, DIMENSION(:, :), POINTER                  :: dist1_prv, dist2_prv
      TYPE(cp_1d_i_p_type), DIMENSION(:), POINTER        :: local_particle_col, local_particle_row
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (blacs_env, local_particle_col, local_particle_row, para_env)

      CALL timeset(routineN, handle)

      CPASSERT(PRESENT(mp_comm_2d) .OR. PRESENT(blacs_env_ext))

      IF (PRESENT(mp_comm_2d)) THEN
         mp_dims = mp_comm_2d%num_pe_cart
         mp_coor = mp_comm_2d%mepos_cart
         ALLOCATE (para_env)
         para_env = mp_comm_2d
         CALL cp_blacs_env_create(blacs_env, para_env, &
                                  grid_2d=mp_dims)

         CPASSERT(blacs_env%mepos(1) == mp_coor(1))
         CPASSERT(blacs_env%mepos(2) == mp_coor(2))
         CALL mp_para_env_release(para_env)
      END IF

      IF (PRESENT(blacs_env_ext)) THEN
         blacs_env => blacs_env_ext
         mp_coor(1) = blacs_env%mepos(1)
         mp_coor(2) = blacs_env%mepos(2)
      END IF

      natom = SIZE(particle_set)
      ALLOCATE (dist1_prv(natom, 2), dist2_prv(natom, 2))
      dist1_prv(:, 1) = dist1
      dist2_prv(:, 1) = dist2

      ALLOCATE (local_particle_col(nkind), local_particle_row(nkind))
      ALLOCATE (nparticle_local_row(nkind), nparticle_local_col(nkind))
      nparticle_local_row = 0; nparticle_local_col = 0

      DO iatom = 1, natom
         ikind = particle_set(iatom)%atomic_kind%kind_number

         IF (dist1_prv(iatom, 1) == mp_coor(1)) nparticle_local_row(ikind) = nparticle_local_row(ikind) + 1
         IF (dist2_prv(iatom, 1) == mp_coor(2)) nparticle_local_col(ikind) = nparticle_local_col(ikind) + 1
      END DO

      DO ikind = 1, nkind
         n = nparticle_local_row(ikind)
         ALLOCATE (local_particle_row(ikind)%array(n))

         n = nparticle_local_col(ikind)
         ALLOCATE (local_particle_col(ikind)%array(n))
      END DO

      nparticle_local_row = 0; nparticle_local_col = 0
      DO iatom = 1, natom
         ikind = particle_set(iatom)%atomic_kind%kind_number

         IF (dist1_prv(iatom, 1) == mp_coor(1)) THEN
            nparticle_local_row(ikind) = nparticle_local_row(ikind) + 1
            local_particle_row(ikind)%array(nparticle_local_row(ikind)) = iatom
         END IF
         IF (dist2_prv(iatom, 1) == mp_coor(2)) THEN
            nparticle_local_col(ikind) = nparticle_local_col(ikind) + 1
            local_particle_col(ikind)%array(nparticle_local_col(ikind)) = iatom
         END IF
      END DO

      CALL distribution_2d_create_prv(dist_2d, row_distribution_ptr=dist1_prv, &
                                      col_distribution_ptr=dist2_prv, local_rows_ptr=local_particle_row, &
                                      local_cols_ptr=local_particle_col, blacs_env=blacs_env)

      IF (.NOT. PRESENT(blacs_env_ext)) THEN
         CALL cp_blacs_env_release(blacs_env)
      END IF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief contiguous distribution of weighted elements
!> \param nel ...
!> \param nbin ...
!> \param weights ...
!> \param limits_start ...
!> \param limits_end ...
!> \param dist ...
! **************************************************************************************************
   SUBROUTINE contiguous_tensor_dist(nel, nbin, weights, limits_start, limits_end, dist)
      INTEGER, INTENT(IN)                                :: nel
      INTEGER, INTENT(INOUT)                             :: nbin
      INTEGER, DIMENSION(nel), INTENT(IN)                :: weights
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT), &
         OPTIONAL                                        :: limits_start, limits_end
      INTEGER, DIMENSION(nel), INTENT(OUT), OPTIONAL     :: dist

      INTEGER                                            :: el_end, el_start, end_weight, ibin, &
                                                            nel_div, nel_rem, nel_split, nel_w, &
                                                            w_partialsum
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: lim_e, lim_s

      ALLOCATE (lim_s(nbin), lim_e(nbin))
      lim_s = 0; lim_e = 0

      nel_w = SUM(weights)
      nel_div = nel_w/nbin
      nel_rem = MOD(nel_w, nbin)

      w_partialsum = 0
      el_end = 0
      end_weight = 0
      DO ibin = 1, nbin
         nel_split = nel_div
         IF (ibin <= nel_rem) THEN
            nel_split = nel_split + 1
         END IF
         el_start = el_end + 1
         el_end = el_start
         w_partialsum = w_partialsum + weights(el_end)
         end_weight = end_weight + nel_split
         DO WHILE (w_partialsum < end_weight)
            !IF (ABS(w_partialsum + weights(el_end) - end_weight) > ABS(w_partialsum - end_weight)) EXIT
            el_end = el_end + 1
            w_partialsum = w_partialsum + weights(el_end)
            IF (el_end == nel) EXIT
         END DO

         IF (PRESENT(dist)) dist(el_start:el_end) = ibin - 1
         lim_s(ibin) = el_start
         lim_e(ibin) = el_end

         IF (el_end == nel) EXIT
      END DO

      IF (PRESENT(limits_start) .AND. PRESENT(limits_end)) THEN
         ALLOCATE (limits_start(ibin)); limits_start(:ibin) = lim_s(:ibin)
         ALLOCATE (limits_end(ibin)); limits_end(:ibin) = lim_e(:ibin)
      END IF

      nbin = ibin

   END SUBROUTINE contiguous_tensor_dist

! **************************************************************************************************
!> \brief ...
!> \param t3c Create 3-center tensor with load balanced default distribution.
!> \param dist_1 ...
!> \param dist_2 ...
!> \param dist_3 ...
!> \param pgrid ...
!> \param sizes_1 ...
!> \param sizes_2 ...
!> \param sizes_3 ...
!> \param map1 ...
!> \param map2 ...
!> \param name ...
! **************************************************************************************************
   SUBROUTINE create_3c_tensor(t3c, dist_1, dist_2, dist_3, pgrid, sizes_1, sizes_2, sizes_3, map1, map2, name)
      TYPE(dbt_type), INTENT(OUT)                        :: t3c
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: dist_1, dist_2, dist_3
      TYPE(dbt_pgrid_type), INTENT(IN)                   :: pgrid
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes_1, sizes_2, sizes_3, map1, map2
      CHARACTER(len=*), INTENT(IN)                       :: name

      CHARACTER(len=*), PARAMETER                        :: routineN = 'create_3c_tensor'

      INTEGER                                            :: handle, size_1, size_2, size_3
      INTEGER, DIMENSION(3)                              :: pcoord, pdims
      TYPE(dbt_distribution_type)                        :: dist

      CALL timeset(routineN, handle)

      CALL dbt_mp_environ_pgrid(pgrid, pdims, pcoord)

      size_1 = SIZE(sizes_1)
      size_2 = SIZE(sizes_2)
      size_3 = SIZE(sizes_3)

      ALLOCATE (dist_1(size_1))
      ALLOCATE (dist_2(size_2))
      ALLOCATE (dist_3(size_3))

      CALL dbt_default_distvec(size_1, pdims(1), sizes_1, dist_1)
      CALL dbt_default_distvec(size_2, pdims(2), sizes_2, dist_2)
      CALL dbt_default_distvec(size_3, pdims(3), sizes_3, dist_3)

      CALL dbt_distribution_new(dist, pgrid, dist_1, dist_2, dist_3)
      CALL dbt_create(t3c, name, dist, map1, map2, sizes_1, sizes_2, sizes_3)
      CALL dbt_distribution_destroy(dist)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param t2c ...
!> \param dist_1 ...
!> \param dist_2 ...
!> \param pgrid ...
!> \param sizes_1 ...
!> \param sizes_2 ...
!> \param order ...
!> \param name ...
! **************************************************************************************************
   SUBROUTINE create_2c_tensor(t2c, dist_1, dist_2, pgrid, sizes_1, sizes_2, order, name)
      TYPE(dbt_type), INTENT(OUT)                        :: t2c
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: dist_1, dist_2
      TYPE(dbt_pgrid_type), INTENT(IN)                   :: pgrid
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes_1, sizes_2
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: order
      CHARACTER(len=*), INTENT(IN)                       :: name

      CHARACTER(len=*), PARAMETER                        :: routineN = 'create_2c_tensor'

      INTEGER                                            :: handle, size_1, size_2
      INTEGER, DIMENSION(2)                              :: order_in, pcoord, pdims
      TYPE(dbt_distribution_type)                        :: dist

      CALL timeset(routineN, handle)

      IF (PRESENT(order)) THEN
         order_in = order
      ELSE
         order_in = [1, 2]
      END IF

      CALL dbt_mp_environ_pgrid(pgrid, pdims, pcoord)

      size_1 = SIZE(sizes_1)
      size_2 = SIZE(sizes_2)

      ALLOCATE (dist_1(size_1))
      ALLOCATE (dist_2(size_2))

      CALL dbt_default_distvec(size_1, pdims(1), sizes_1, dist_1)
      CALL dbt_default_distvec(size_2, pdims(2), sizes_2, dist_2)

      CALL dbt_distribution_new(dist, pgrid, dist_1, dist_2)
      CALL dbt_create(t2c, name, dist, [order_in(1)], [order_in(2)], sizes_1, sizes_2)
      CALL dbt_distribution_destroy(dist)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param blk_sizes ...
!> \param blk_sizes_split ...
!> \param max_size ...
! **************************************************************************************************
   SUBROUTINE split_block_sizes(blk_sizes, blk_sizes_split, max_size)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: blk_sizes
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: blk_sizes_split
      INTEGER, INTENT(IN)                                :: max_size

      INTEGER                                            :: blk_remainder, i, isplit, isplit_sum, &
                                                            nsplit

      isplit_sum = 0
      DO i = 1, SIZE(blk_sizes)
         nsplit = (blk_sizes(i) + max_size - 1)/max_size
         isplit_sum = isplit_sum + nsplit
      END DO

      ALLOCATE (blk_sizes_split(isplit_sum))

      isplit_sum = 0
      DO i = 1, SIZE(blk_sizes)
         nsplit = (blk_sizes(i) + max_size - 1)/max_size
         blk_remainder = blk_sizes(i)
         DO isplit = 1, nsplit
            isplit_sum = isplit_sum + 1
            blk_sizes_split(isplit_sum) = MIN(max_size, blk_remainder)
            blk_remainder = blk_remainder - max_size
         END DO
      END DO

   END SUBROUTINE split_block_sizes

! **************************************************************************************************
!> \brief ...
!> \param atomic_kind_set ...
!> \param basis ...
!> \param min_blk_size ...
!> \param pgf_blk_sizes ...
! **************************************************************************************************
   SUBROUTINE pgf_block_sizes(atomic_kind_set, basis, min_blk_size, pgf_blk_sizes)
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis
      INTEGER, INTENT(IN)                                :: min_blk_size
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: pgf_blk_sizes

      INTEGER                                            :: blk_count, blk_count_prev, blk_size, &
                                                            iatom, ikind, iset, natom, nblk, nset
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of, pgf_blk_sizes_tmp
      INTEGER, DIMENSION(:), POINTER                     :: nsgf_set

      CALL get_atomic_kind_set(atomic_kind_set, natom=natom, kind_of=kind_of)

      nblk = 0
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_gto_basis_set(basis(ikind)%gto_basis_set, nset=nset)
         nblk = nblk + nset
      END DO

      ALLOCATE (pgf_blk_sizes_tmp(nblk)); pgf_blk_sizes_tmp = 0

      blk_count = 0
      blk_size = 0
      DO iatom = 1, natom
         blk_count_prev = blk_count
         ikind = kind_of(iatom)
         CALL get_gto_basis_set(basis(ikind)%gto_basis_set, nset=nset, nsgf_set=nsgf_set)
         DO iset = 1, nset
            blk_size = blk_size + nsgf_set(iset)
            IF (blk_size >= min_blk_size) THEN
               blk_count = blk_count + 1
               pgf_blk_sizes_tmp(blk_count) = pgf_blk_sizes_tmp(blk_count) + blk_size
               blk_size = 0
            END IF
         END DO
         IF (blk_size > 0) THEN
            IF (blk_count == blk_count_prev) blk_count = blk_count + 1
            pgf_blk_sizes_tmp(blk_count) = pgf_blk_sizes_tmp(blk_count) + blk_size
            blk_size = 0
         END IF
      END DO

      ALLOCATE (pgf_blk_sizes(blk_count))
      pgf_blk_sizes(:) = pgf_blk_sizes_tmp(:blk_count)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param sizes ...
!> \param nbatches ...
!> \param starts_array ...
!> \param ends_array ...
!> \param starts_array_block ...
!> \param ends_array_block ...
! **************************************************************************************************
   SUBROUTINE create_tensor_batches(sizes, nbatches, starts_array, ends_array, &
                                    starts_array_block, ends_array_block)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes
      INTEGER, INTENT(INOUT)                             :: nbatches
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: starts_array, ends_array, &
                                                            starts_array_block, ends_array_block

      INTEGER                                            :: bsum, imem, nblocks

      nblocks = SIZE(sizes)

      CALL contiguous_tensor_dist(nblocks, nbatches, sizes, limits_start=starts_array_block, limits_end=ends_array_block)

      ALLOCATE (starts_array(nbatches))
      ALLOCATE (ends_array(nbatches))

      bsum = 0
      DO imem = 1, nbatches
         starts_array(imem) = bsum + 1
         bsum = bsum + SUM(sizes(starts_array_block(imem):ends_array_block(imem)))
         ends_array(imem) = bsum
      END DO
   END SUBROUTINE

END MODULE
