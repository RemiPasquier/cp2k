!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief utils to manipulate splines on the regular grid of a pw
!> \par History
!>      01.2014 move routines related to input_section_types to separate file.
!> \author Ole Schuett
! **************************************************************************************************
MODULE cp_spline_utils
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE pw_methods,                      ONLY: do_standard_sum,&
                                              pw_axpy,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_spline_utils,                 ONLY: &
        add_coarse2fine, add_fine2coarse, find_coeffs, pw_spline_do_precond, &
        pw_spline_precond_create, pw_spline_precond_release, pw_spline_precond_set_kind, &
        pw_spline_precond_type, spl3_1d_transf_border1, spl3_1d_transf_coeffs, spl3_nopbc, &
        spl3_nopbct, spl3_pbc
   USE pw_types,                        ONLY: pw_r3d_rs_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_spline_utils'

   PUBLIC ::  pw_prolongate_s3, pw_restrict_s3

   ! input constants
   INTEGER, PARAMETER, PUBLIC               :: pw_interp = 1, &
                                               spline3_nopbc_interp = 2, &
                                               spline3_pbc_interp = 3

CONTAINS

! **************************************************************************************************
!> \brief restricts the function from a fine grid to a coarse one
!> \param pw_fine_in the fine grid
!> \param pw_coarse_out the coarse grid
!> \param coarse_pool ...
!> \param param_section ...
!> \author fawzi
!> \note
!>      extremely slow (but correct) version
! **************************************************************************************************
   SUBROUTINE pw_restrict_s3(pw_fine_in, pw_coarse_out, coarse_pool, param_section)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: pw_fine_in
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: pw_coarse_out
      TYPE(pw_pool_type), POINTER                        :: coarse_pool
      TYPE(section_vals_type), POINTER                   :: param_section

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_restrict_s3'

      INTEGER                                            :: aint_precond, handle, interp_kind, &
                                                            max_iter, precond_kind
      INTEGER, DIMENSION(2, 3)                           :: bo
      INTEGER, SAVE                                      :: ifile = 0
      LOGICAL                                            :: pbc, safe_computation, success
      REAL(kind=dp)                                      :: eps_r, eps_x
      TYPE(pw_r3d_rs_type)                               :: coeffs, values
      TYPE(pw_spline_precond_type)                       :: precond

      ifile = ifile + 1
      CALL timeset(routineN, handle)
      CALL section_vals_val_get(param_section, "safe_computation", &
                                l_val=safe_computation)
      CALL section_vals_val_get(param_section, "aint_precond", &
                                i_val=aint_precond)
      CALL section_vals_val_get(param_section, "precond", &
                                i_val=precond_kind)
      CALL section_vals_val_get(param_section, "max_iter", &
                                i_val=max_iter)
      CALL section_vals_val_get(param_section, "eps_r", &
                                r_val=eps_r)
      CALL section_vals_val_get(param_section, "eps_x", &
                                r_val=eps_x)
      CALL section_vals_val_get(param_section, "kind", &
                                i_val=interp_kind)

      pbc = (interp_kind == spline3_pbc_interp)
      CPASSERT(pbc .OR. interp_kind == spline3_nopbc_interp)
      bo = pw_coarse_out%pw_grid%bounds_local
      CALL coarse_pool%create_pw(values)
      CALL pw_zero(values)

      CALL add_fine2coarse(fine_values_pw=pw_fine_in, &
                           coarse_coeffs_pw=values, &
                           weights_1d=spl3_1d_transf_coeffs/2._dp, w_border0=0.5_dp, &
                           w_border1=spl3_1d_transf_border1/2._dp, pbc=pbc, &
                           safe_computation=safe_computation)

      CALL coarse_pool%create_pw(coeffs)
      CALL pw_spline_precond_create(precond, precond_kind=aint_precond, &
                                    pool=coarse_pool, pbc=pbc, transpose=.TRUE.)
      CALL pw_spline_do_precond(precond, values, coeffs)
      CALL pw_spline_precond_set_kind(precond, precond_kind)
      IF (pbc) THEN
         success = find_coeffs(values=values, coeffs=coeffs, &
                               linOp=spl3_pbc, preconditioner=precond, pool=coarse_pool, &
                               eps_r=eps_r, eps_x=eps_x, max_iter=max_iter, sumtype=do_standard_sum)
      ELSE
         success = find_coeffs(values=values, coeffs=coeffs, &
                               linOp=spl3_nopbct, preconditioner=precond, pool=coarse_pool, &
                               eps_r=eps_r, eps_x=eps_x, max_iter=max_iter, sumtype=do_standard_sum)
      END IF
      CALL pw_spline_precond_release(precond)

      CALL pw_zero(pw_coarse_out)
      CALL pw_axpy(coeffs, pw_coarse_out)

      CALL coarse_pool%give_back_pw(values)
      CALL coarse_pool%give_back_pw(coeffs)
      CALL timestop(handle)

   END SUBROUTINE pw_restrict_s3

! **************************************************************************************************
!> \brief prolongates a function from a coarse grid into a fine one
!> \param pw_coarse_in the coarse grid
!> \param pw_fine_out the fine grid
!> \param coarse_pool ...
!> \param param_section ...
!> \author fawzi
!> \note
!>      extremely slow (but correct) version
! **************************************************************************************************
   SUBROUTINE pw_prolongate_s3(pw_coarse_in, pw_fine_out, coarse_pool, &
                               param_section)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: pw_coarse_in, pw_fine_out
      TYPE(pw_pool_type), POINTER                        :: coarse_pool
      TYPE(section_vals_type), POINTER                   :: param_section

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_prolongate_s3'

      INTEGER                                            :: aint_precond, handle, interp_kind, &
                                                            max_iter, precond_kind
      INTEGER, DIMENSION(2, 3)                           :: bo
      INTEGER, SAVE                                      :: ifile = 0
      LOGICAL                                            :: pbc, safe_computation, success
      REAL(kind=dp)                                      :: eps_r, eps_x
      TYPE(pw_r3d_rs_type)                               :: coeffs
      TYPE(pw_spline_precond_type)                       :: precond

      ifile = ifile + 1
      CALL timeset(routineN, handle)
      CALL coarse_pool%create_pw(coeffs)
      bo = pw_coarse_in%pw_grid%bounds_local
      CALL section_vals_val_get(param_section, "safe_computation", &
                                l_val=safe_computation)
      CALL section_vals_val_get(param_section, "aint_precond", &
                                i_val=aint_precond)
      CALL section_vals_val_get(param_section, "precond", &
                                i_val=precond_kind)
      CALL section_vals_val_get(param_section, "max_iter", &
                                i_val=max_iter)
      CALL section_vals_val_get(param_section, "eps_r", &
                                r_val=eps_r)
      CALL section_vals_val_get(param_section, "eps_x", &
                                r_val=eps_x)
      CALL section_vals_val_get(param_section, "kind", &
                                i_val=interp_kind)

      pbc = (interp_kind == spline3_pbc_interp)
      CPASSERT(pbc .OR. interp_kind == spline3_nopbc_interp)
      CALL pw_spline_precond_create(precond, precond_kind=aint_precond, &
                                    pool=coarse_pool, pbc=pbc, transpose=.FALSE.)
      CALL pw_spline_do_precond(precond, pw_coarse_in, coeffs)
      CALL pw_spline_precond_set_kind(precond, precond_kind)
      IF (pbc) THEN
         success = find_coeffs(values=pw_coarse_in, coeffs=coeffs, &
                               linOp=spl3_pbc, preconditioner=precond, pool=coarse_pool, &
                               eps_r=eps_r, eps_x=eps_x, &
                               max_iter=max_iter, sumtype=do_standard_sum)
      ELSE
         success = find_coeffs(values=pw_coarse_in, coeffs=coeffs, &
                               linOp=spl3_nopbc, preconditioner=precond, pool=coarse_pool, &
                               eps_r=eps_r, eps_x=eps_x, &
                               max_iter=max_iter, sumtype=do_standard_sum)
      END IF
      CPASSERT(success)
      CALL pw_spline_precond_release(precond)

      CALL add_coarse2fine(coarse_coeffs_pw=coeffs, &
                           fine_values_pw=pw_fine_out, &
                           weights_1d=spl3_1d_transf_coeffs, &
                           w_border0=1._dp, &
                           w_border1=spl3_1d_transf_border1, &
                           pbc=pbc, safe_computation=safe_computation)

      CALL coarse_pool%give_back_pw(coeffs)

      CALL timestop(handle)

   END SUBROUTINE pw_prolongate_s3

END MODULE cp_spline_utils
