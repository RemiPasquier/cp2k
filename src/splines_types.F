!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines for handling splines_types
!> \par History
!>      2001-09-21-HAF added this doc entry and changed formatting
!> \author various
! **************************************************************************************************
MODULE splines_types

   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'splines_types'

   PUBLIC :: spline_env_release, spline_environment_type
   PUBLIC :: spline_env_create, spline_data_p_type
   PUBLIC :: spline_data_create, spline_data_p_copy
   PUBLIC :: spline_data_retain, spline_data_p_retain
   PUBLIC :: spline_data_release, spline_data_p_release
   PUBLIC :: spline_factor_copy, spline_factor_create, spline_factor_release
   PUBLIC :: spline_data_type ! the data structure for spline table
   PUBLIC :: spline_factor_type ! the multiplicative factors for splines

! **************************************************************************************************
!> \brief Data-structure that holds all needed information about
!>      a specific spline interpolation.
!> \par History
!>      2001-09-19-HAF added this doc entry and changed formatting
!> \author unknown
! **************************************************************************************************
   TYPE spline_data_type
      INTEGER :: ref_count = -1
      REAL(KIND=dp), POINTER :: y(:) => NULL() ! the function values y(x)
      REAL(KIND=dp), POINTER :: y2(:) => NULL() ! the 2nd derivative via interpolation
      INTEGER                 :: n = -1 ! dimension of above arrays
      ! not used if uniform increments
      REAL(KIND=dp)          :: h = -1.0_dp ! uniform increment of x if applicable
      REAL(KIND=dp)          :: invh = -1.0_dp ! inverse of h
      REAL(KIND=dp)          :: h26 = -1.0_dp ! 1/6 * h**2 if uniform increments
      ! 1/6 otherwise
      REAL(KIND=dp)          :: x1 = -1.0_dp ! starting x value if uniform incr.
      REAL(KIND=dp)          :: xn = -1.0_dp ! end x value if uniform incr.
   END TYPE spline_data_type

! **************************************************************************************************
   TYPE spline_data_p_type
      TYPE(spline_data_type), POINTER :: spline_data => NULL()
   END TYPE spline_data_p_type

! **************************************************************************************************
   TYPE spline_data_pp_type
      TYPE(spline_data_p_type), POINTER, DIMENSION(:)  :: spl_p => NULL()
   END TYPE spline_data_pp_type

! **************************************************************************************************
   TYPE spline_environment_type
      TYPE(spline_data_pp_type), POINTER, DIMENSION(:) :: spl_pp => NULL()
      INTEGER, POINTER, DIMENSION(:, :) :: spltab => NULL()
   END TYPE spline_environment_type

! **************************************************************************************************
   TYPE spline_factor_type
      REAL(KIND=dp) :: rcutsq_f = -1.0_dp, cutoff = -1.0_dp
      REAL(KIND=dp), DIMENSION(:), POINTER :: rscale => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER :: fscale => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER :: dscale => NULL()
   END TYPE spline_factor_type

CONTAINS

! **************************************************************************************************
!> \brief releases spline_env
!> \param spline_env ...
!> \author unknown
! **************************************************************************************************
   SUBROUTINE spline_env_release(spline_env)
      TYPE(spline_environment_type), INTENT(INOUT)       :: spline_env

      INTEGER                                            :: i
      TYPE(spline_data_p_type), DIMENSION(:), POINTER    :: spl_p

      DEALLOCATE (spline_env%spltab)
      DO i = 1, SIZE(spline_env%spl_pp)
         spl_p => spline_env%spl_pp(i)%spl_p
         CALL spline_data_p_release(spl_p)
      END DO
      DEALLOCATE (spline_env%spl_pp)

   END SUBROUTINE spline_env_release

! **************************************************************************************************
!> \brief releases spline_data
!> \param spline_data ...
!> \author CJM
! **************************************************************************************************
   SUBROUTINE spline_data_release(spline_data)
      TYPE(spline_data_type), POINTER                    :: spline_data

      IF (ASSOCIATED(spline_data)) THEN
         CPASSERT(spline_data%ref_count > 0)
         spline_data%ref_count = spline_data%ref_count - 1
         IF (spline_data%ref_count < 1) THEN
            IF (ASSOCIATED(spline_data%y)) THEN
               DEALLOCATE (spline_data%y)
            END IF
            IF (ASSOCIATED(spline_data%y2)) THEN
               DEALLOCATE (spline_data%y2)
            END IF
            DEALLOCATE (spline_data)
         END IF
      END IF
   END SUBROUTINE spline_data_release

! **************************************************************************************************
!> \brief releases spline_data_p
!> \param spl_p ...
!> \author CJM
! **************************************************************************************************
   SUBROUTINE spline_data_p_release(spl_p)
      TYPE(spline_data_p_type), DIMENSION(:), POINTER    :: spl_p

      INTEGER                                            :: i
      LOGICAL                                            :: release_kind

      IF (ASSOCIATED(spl_p)) THEN
         release_kind = .TRUE.
         DO i = 1, SIZE(spl_p)
            CALL spline_data_release(spl_p(i)%spline_data)
            release_kind = release_kind .AND. (.NOT. ASSOCIATED(spl_p(i)%spline_data))
         END DO
         IF (release_kind) THEN
            DEALLOCATE (spl_p)
         END IF
      END IF

   END SUBROUTINE spline_data_p_release

! **************************************************************************************************
!> \brief retains spline_env
!> \param spline_data ...
!> \author CJM
! **************************************************************************************************
   SUBROUTINE spline_data_retain(spline_data)
      TYPE(spline_data_type), POINTER                    :: spline_data

      CPASSERT(ASSOCIATED(spline_data))
      CPASSERT(spline_data%ref_count > 0)
      spline_data%ref_count = spline_data%ref_count + 1
   END SUBROUTINE spline_data_retain

! **************************************************************************************************
!> \brief retains spline_data_p_type
!> \param spl_p ...
!> \author CJM
! **************************************************************************************************
   SUBROUTINE spline_data_p_retain(spl_p)
      TYPE(spline_data_p_type), DIMENSION(:), POINTER    :: spl_p

      INTEGER                                            :: i

      CPASSERT(ASSOCIATED(spl_p))
      DO i = 1, SIZE(spl_p)
         CALL spline_data_retain(spl_p(i)%spline_data)
      END DO
   END SUBROUTINE spline_data_p_retain

! **************************************************************************************************
!> \brief Data-structure that holds all needed information about
!>      a specific spline interpolation.
!> \param spline_env ...
!> \param ntype ...
!> \param ntab_in ...
!> \par History
!>      2001-09-19-HAF added this doc entry and changed formatting
!> \author unknown
! **************************************************************************************************
   SUBROUTINE spline_env_create(spline_env, ntype, ntab_in)
      TYPE(spline_environment_type), INTENT(OUT)         :: spline_env
      INTEGER, INTENT(IN)                                :: ntype
      INTEGER, INTENT(IN), OPTIONAL                      :: ntab_in

      CHARACTER(len=*), PARAMETER                        :: routineN = 'spline_env_create'

      INTEGER                                            :: handle, i, isize, j, ntab

      CALL timeset(routineN, handle)

      NULLIFY (spline_env%spl_pp)
      NULLIFY (spline_env%spltab)
      ! Allocate the number of spline data tables (upper triangular)
      IF (PRESENT(ntab_in)) THEN
         ntab = ntab_in
      ELSE
         ntab = (ntype*ntype + ntype)/2
      END IF
      ALLOCATE (spline_env%spl_pp(ntab))

      ALLOCATE (spline_env%spltab(ntype, ntype))

      DO i = 1, ntab
         NULLIFY (spline_env%spl_pp(i)%spl_p)
         isize = 1
         ALLOCATE (spline_env%spl_pp(i)%spl_p(isize))
         DO j = 1, SIZE(spline_env%spl_pp(i)%spl_p)
            CALL spline_data_create(spline_env%spl_pp(i)%spl_p(j)%spline_data)
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE spline_env_create

! **************************************************************************************************
!> \brief Copy Data-structure of spline_data_p_type
!> \param spl_p_source ...
!> \param spl_p_dest ...
!> \author teo 06.2007
! **************************************************************************************************
   SUBROUTINE spline_data_p_copy(spl_p_source, spl_p_dest)
      TYPE(spline_data_p_type), DIMENSION(:), POINTER    :: spl_p_source, spl_p_dest

      INTEGER                                            :: i, nsized, nsizes

      CPASSERT(ASSOCIATED(spl_p_source))
      nsizes = SIZE(spl_p_source)
      IF (.NOT. ASSOCIATED(spl_p_dest)) THEN
         ALLOCATE (spl_p_dest(nsizes))
         DO i = 1, nsizes
            NULLIFY (spl_p_dest(i)%spline_data)
         END DO
      ELSE
         nsized = SIZE(spl_p_dest)
         CPASSERT(nsizes == nsized)
         DO i = 1, nsizes
            CALL spline_data_release(spl_p_dest(i)%spline_data)
         END DO
      END IF
      DO i = 1, nsizes
         CALL spline_data_copy(spl_p_source(i)%spline_data, spl_p_dest(i)%spline_data)
      END DO
   END SUBROUTINE spline_data_p_copy

! **************************************************************************************************
!> \brief Copy Data-structure that constains spline table
!> \param spline_data_source ...
!> \param spline_data_dest ...
!> \author teo 11.2005
! **************************************************************************************************
   SUBROUTINE spline_data_copy(spline_data_source, spline_data_dest)
      TYPE(spline_data_type), POINTER                    :: spline_data_source, spline_data_dest

      CPASSERT(ASSOCIATED(spline_data_source))
      IF (.NOT. ASSOCIATED(spline_data_dest)) CALL spline_data_create(spline_data_dest)

      spline_data_dest%ref_count = spline_data_source%ref_count
      spline_data_dest%n = spline_data_source%n
      spline_data_dest%h = spline_data_source%h
      spline_data_dest%invh = spline_data_source%invh
      spline_data_dest%h26 = spline_data_source%h26
      spline_data_dest%x1 = spline_data_source%x1
      spline_data_dest%xn = spline_data_source%xn
      IF (ASSOCIATED(spline_data_source%y)) THEN
         ALLOCATE (spline_data_dest%y(SIZE(spline_data_source%y)))
         spline_data_dest%y = spline_data_source%y
      END IF
      IF (ASSOCIATED(spline_data_source%y2)) THEN
         ALLOCATE (spline_data_dest%y2(SIZE(spline_data_source%y2)))
         spline_data_dest%y2 = spline_data_source%y2
      END IF
   END SUBROUTINE spline_data_copy

! **************************************************************************************************
!> \brief Data-structure that constains spline table
!> \param spline_data ...
!> \author unknown
! **************************************************************************************************
   SUBROUTINE spline_data_create(spline_data)
      TYPE(spline_data_type), POINTER                    :: spline_data

      ALLOCATE (spline_data)
      spline_data%ref_count = 1
      NULLIFY (spline_data%y)
      NULLIFY (spline_data%y2)
   END SUBROUTINE spline_data_create

! **************************************************************************************************
!> \brief releases spline_factor
!> \param spline_factor ...
!> \author teo
! **************************************************************************************************
   SUBROUTINE spline_factor_release(spline_factor)
      TYPE(spline_factor_type), POINTER                  :: spline_factor

      IF (ASSOCIATED(spline_factor)) THEN
         IF (ASSOCIATED(spline_factor%rscale)) THEN
            DEALLOCATE (spline_factor%rscale)
         END IF
         IF (ASSOCIATED(spline_factor%fscale)) THEN
            DEALLOCATE (spline_factor%fscale)
         END IF
         IF (ASSOCIATED(spline_factor%dscale)) THEN
            DEALLOCATE (spline_factor%dscale)
         END IF
         DEALLOCATE (spline_factor)
      END IF
   END SUBROUTINE spline_factor_release

! **************************************************************************************************
!> \brief releases spline_factor
!> \param spline_factor ...
!> \author teo
! **************************************************************************************************
   SUBROUTINE spline_factor_create(spline_factor)
      TYPE(spline_factor_type), POINTER                  :: spline_factor

      CPASSERT(.NOT. ASSOCIATED(spline_factor))
      ALLOCATE (spline_factor)
      ALLOCATE (spline_factor%rscale(1))
      ALLOCATE (spline_factor%fscale(1))
      ALLOCATE (spline_factor%dscale(1))
      spline_factor%rscale = 1.0_dp
      spline_factor%fscale = 1.0_dp
      spline_factor%dscale = 1.0_dp
      spline_factor%rcutsq_f = 1.0_dp
      spline_factor%cutoff = 0.0_dp
   END SUBROUTINE spline_factor_create

! **************************************************************************************************
!> \brief releases spline_factor
!> \param spline_factor_source ...
!> \param spline_factor_dest ...
!> \author teo
! **************************************************************************************************
   SUBROUTINE spline_factor_copy(spline_factor_source, spline_factor_dest)
      TYPE(spline_factor_type), POINTER                  :: spline_factor_source, spline_factor_dest

      INTEGER                                            :: isize, jsize, ksize

      IF (ASSOCIATED(spline_factor_dest)) CALL spline_factor_release(spline_factor_dest)
      IF (ASSOCIATED(spline_factor_source)) THEN
         isize = SIZE(spline_factor_source%rscale)
         jsize = SIZE(spline_factor_source%fscale)
         ksize = SIZE(spline_factor_source%dscale)
         CPASSERT(isize == jsize)
         CPASSERT(isize == ksize)
         CALL spline_factor_create(spline_factor_dest)
         spline_factor_dest%rscale = spline_factor_source%rscale
         spline_factor_dest%fscale = spline_factor_source%fscale
         spline_factor_dest%dscale = spline_factor_source%dscale
         spline_factor_dest%rcutsq_f = spline_factor_source%rcutsq_f
         spline_factor_dest%cutoff = spline_factor_source%cutoff
      END IF
   END SUBROUTINE spline_factor_copy

END MODULE splines_types
