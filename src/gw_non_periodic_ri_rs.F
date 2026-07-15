!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief GW using RI-RS Approximation for molecules
!> \par History
!>      04.2026 created [Ritaj Tyagi]
! **************************************************************************************************
MODULE gw_non_periodic_ri_rs
   USE ai_contraction_sphi,             ONLY: ab_contract
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_binary_read, dbcsr_binary_write, dbcsr_copy, dbcsr_create, &
        dbcsr_deallocate_matrix, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_block_p, &
        dbcsr_get_info, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_deallocate_matrix_set,&
                                              max_elements_per_block
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                              cp_fm_uplo_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert,&
                                              cp_fm_cholesky_solve
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_should_output
   USE gamma,                           ONLY: init_md_ftable
   USE gw_integrals,                    ONLY: build_3c_integral_block_ctx,&
                                              gw_3c_ctx_create,&
                                              gw_3c_ctx_release,&
                                              gw_3c_ctx_type,&
                                              gw_3c_ws_create,&
                                              gw_3c_ws_release,&
                                              gw_3c_ws_type
   USE gw_large_cell_gamma,             ONLY: &
        Fourier_transform_w_to_t, G_occ_vir, compute_QP_energies, compute_fm_chi_Gamma_freq, &
        create_fm_W_MIC_time, delete_unnecessary_files, fill_fm_Sigma_c_Gamma_time, fm_write, &
        multiply_fm_W_MIC_time_with_Minv_Gamma
   USE gw_utils,                        ONLY: de_init_bs_env
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              rtp_method_bse
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE kpoint_coulomb_2c,               ONLY: build_2c_coulomb_matrix_kp
   USE lebedev,                         ONLY: get_number_of_lebedev_grid,&
                                              lebedev_grid
   USE libint_2c_3c,                    ONLY: eri_2center,&
                                              libint_potential_type
   USE libint_wrapper,                  ONLY: cp_libint_cleanup_2eri,&
                                              cp_libint_init_2eri,&
                                              cp_libint_set_contrdepth,&
                                              cp_libint_t
   USE machine,                         ONLY: m_walltime
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type,&
                                              rirs_grid_escalation_type,&
                                              rirs_grid_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE string_utilities,                ONLY: uppercase
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_non_periodic_ri_rs'

   INTEGER, PARAMETER, PRIVATE :: ri_rs_topk_size = 30, &
                                  ri_rs_topk_nlists = 4, &
                                  ri_rs_topk_e3c = 1, &
                                  ri_rs_topk_mbare = 2, &
                                  ri_rs_topk_mtrunc = 3, &
                                  ri_rs_topk_e3crel = 4
   REAL(KIND=dp), PARAMETER, PRIVATE :: ri_rs_j_full_max_bytes = 2.0E9_dp  ! ~2 GB/rank cap for replicated R_full
   ! Absolute |exact| floor, not a running/global max|v| scale: this single streaming pass over
   ! integrals ranks before the global max is known.
   REAL(KIND=dp), PARAMETER, PRIVATE :: ri_rs_relfloor = 1.0E-2_dp

   TYPE :: ri_rs_topk_type
      REAL(KIND=dp) :: val(ri_rs_topk_size, ri_rs_topk_nlists) = -1.0_dp
      REAL(KIND=dp) :: vref(ri_rs_topk_size, ri_rs_topk_nlists) = 0.0_dp
      INTEGER       :: idx(3, ri_rs_topk_size, ri_rs_topk_nlists) = 0
      INTEGER       :: atom(3, ri_rs_topk_size, ri_rs_topk_nlists) = 0
   END TYPE ri_rs_topk_type

   PUBLIC :: gw_calc_non_periodic_ri_rs, ri_rs_grid_assembler, &
             get_basis_offsets, precompute_ri_rs_radii, solve_D_lp_distributed, &
             ri_rs_e3c_accumulate, ri_rs_fit_error_finalize, ri_rs_j_full_compute, &
             ri_rs_topk_type

CONTAINS

! **************************************************************************************************
!> \brief Insert one absolute value into a fixed-size per-rank top-k diagnostic list.
! **************************************************************************************************
   SUBROUTINE ri_rs_topk_insert(topk, list_id, val, vref, i1, i2, i3, a1, a2, a3)

      TYPE(ri_rs_topk_type), INTENT(INOUT)               :: topk
      INTEGER, INTENT(IN)                                :: list_id, i1, i2, i3, a1, a2, a3
      REAL(KIND=dp), INTENT(IN)                          :: val, vref

      INTEGER                                            :: k, min_slot
      REAL(KIND=dp)                                      :: abs_val, min_val

      IF (list_id < 1 .OR. list_id > ri_rs_topk_nlists) RETURN

      abs_val = ABS(val)
      IF (abs_val <= 0.0_dp) RETURN

      min_slot = 1
      min_val = topk%val(1, list_id)
      DO k = 2, ri_rs_topk_size
         IF (topk%val(k, list_id) < min_val) THEN
            min_slot = k
            min_val = topk%val(k, list_id)
         END IF
      END DO

      IF (abs_val > min_val) THEN
         topk%val(min_slot, list_id) = abs_val
         topk%vref(min_slot, list_id) = vref
         topk%idx(:, min_slot, list_id) = [i1, i2, i3]
         topk%atom(:, min_slot, list_id) = [a1, a2, a3]
      END IF

   END SUBROUTINE ri_rs_topk_insert

! **************************************************************************************************
!> \brief Update one per-atom max diagnostic slot.
! **************************************************************************************************
   SUBROUTINE ri_rs_atommax_update(atommax, idxmax, iat, val, i1, i2, i3)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: atommax
      INTEGER, DIMENSION(:, :), INTENT(INOUT)            :: idxmax
      INTEGER, INTENT(IN)                                :: iat, i1, i2, i3
      REAL(KIND=dp), INTENT(IN)                          :: val

      IF (iat < 1 .OR. iat > SIZE(atommax)) RETURN

      IF (val > atommax(iat)) THEN
         atommax(iat) = val
         idxmax(:, iat) = [i1, i2, i3]
      END IF

   END SUBROUTINE ri_rs_atommax_update

! **************************************************************************************************
!> \brief Reduce one fixed-size local top-k list into global rank order and print it.
! **************************************************************************************************
   SUBROUTINE ri_rs_topk_reduce_print(bs_env, para_env, particle_set, atomic_kind_set, topk, &
                                      list_id, title, third_label)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(ri_rs_topk_type), INTENT(IN)                  :: topk
      INTEGER, INTENT(IN)                                :: list_id
      CHARACTER(LEN=*), INTENT(IN)                       :: title, third_label

      CHARACTER(LEN=2)                                   :: sym(3)
      INTEGER                                            :: ik, k, m, src, slot
      INTEGER, DIMENSION(6)                              :: loc_buf
      REAL(KIND=dp)                                      :: global_val, local_val, val_buf, vref_buf
      TYPE(ri_rs_topk_type)                              :: topk_work

      topk_work = topk

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A)') TRIM(title)
      END IF

      DO m = 1, ri_rs_topk_size
         local_val = MAXVAL(topk_work%val(:, list_id))
         global_val = local_val
         CALL para_env%max(global_val)
         IF (global_val <= 0.0_dp) EXIT

         src = -1
         IF (local_val == global_val) src = para_env%mepos
         CALL para_env%max(src)

         loc_buf(:) = 0
         val_buf = 0.0_dp
         vref_buf = 0.0_dp
         IF (para_env%mepos == src) THEN
            slot = 1
            DO k = 1, ri_rs_topk_size
               IF (topk_work%val(k, list_id) == global_val) THEN
                  slot = k
                  EXIT
               END IF
            END DO
            val_buf = topk_work%val(slot, list_id)
            vref_buf = topk_work%vref(slot, list_id)
            loc_buf(1:3) = topk_work%idx(:, slot, list_id)
            loc_buf(4:6) = topk_work%atom(:, slot, list_id)
            topk_work%val(slot, list_id) = -1.0_dp
         END IF
         CALL para_env%bcast(val_buf, src)
         CALL para_env%bcast(vref_buf, src)
         CALL para_env%bcast(loc_buf, src)

         sym(:) = "??"
         DO k = 1, 3
            IF (loc_buf(k + 3) >= 1) THEN
               ik = particle_set(loc_buf(k + 3))%atomic_kind%kind_number
               sym(k) = atomic_kind_set(ik)%element_symbol
            END IF
         END DO

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, &
                   '(T4,A,I0,A,ES11.3,A,I0,A,I0,A,A,A,I0,A,3(2X,I0,1X,A2),A,ES11.3)') &
               '#', m, ' val=', val_buf, ' (mu=', loc_buf(1), ' nu=', loc_buf(2), &
               ' ', TRIM(third_label), '=', loc_buf(3), ') atoms', &
               loc_buf(4), sym(1), loc_buf(5), sym(2), loc_buf(6), sym(3), ' |v|= ', vref_buf
         END IF
      END DO

   END SUBROUTINE ri_rs_topk_reduce_print

! **************************************************************************************************
!> \brief Reduce and print per-atom fit-error maxima.
! **************************************************************************************************
   SUBROUTINE ri_rs_atommax_reduce_print(bs_env, para_env, particle_set, atomic_kind_set, &
                                         atommax_e3c, atommax_e3c_idx, &
                                         atommax_e3crel, atommax_e3crel_idx, &
                                         atommax_j, atommax_j_idx)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: atommax_e3c, atommax_e3crel, &
                                                            atommax_j
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: atommax_e3c_idx, atommax_e3crel_idx, &
                                                            atommax_j_idx

      CHARACTER(LEN=2)                                   :: sym
      INTEGER                                            :: iat, ik, natom, src_e3c, src_e3crel, src_j
      INTEGER, DIMENSION(3)                              :: idx_e3c, idx_e3crel, idx_j
      REAL(KIND=dp)                                      :: e3c_glob, e3c_local, e3c_print, &
                                                            e3crel_glob, e3crel_local, e3crel_print, &
                                                            j_glob, j_local, j_print

      natom = SIZE(particle_set)

      DO iat = 1, natom
         e3c_local = atommax_e3c(iat)
         e3c_glob = e3c_local
         CALL para_env%max(e3c_glob)

         src_e3c = -1
         IF (e3c_local == e3c_glob) src_e3c = para_env%mepos
         CALL para_env%max(src_e3c)

         idx_e3c(:) = 0
         IF (para_env%mepos == src_e3c) idx_e3c(:) = atommax_e3c_idx(:, iat)
         CALL para_env%bcast(idx_e3c, src_e3c)

         IF (e3c_glob < 0.0_dp) THEN
            e3c_print = 0.0_dp
            idx_e3c(:) = 0
         ELSE
            e3c_print = e3c_glob
         END IF

         e3crel_local = atommax_e3crel(iat)
         e3crel_glob = e3crel_local
         CALL para_env%max(e3crel_glob)

         src_e3crel = -1
         IF (e3crel_local == e3crel_glob) src_e3crel = para_env%mepos
         CALL para_env%max(src_e3crel)

         idx_e3crel(:) = 0
         IF (para_env%mepos == src_e3crel) idx_e3crel(:) = atommax_e3crel_idx(:, iat)
         CALL para_env%bcast(idx_e3crel, src_e3crel)

         IF (e3crel_glob < 0.0_dp) THEN
            e3crel_print = 0.0_dp
            idx_e3crel(:) = 0
         ELSE
            e3crel_print = e3crel_glob
         END IF

         j_local = atommax_j(iat)
         j_glob = j_local
         CALL para_env%max(j_glob)

         src_j = -1
         IF (j_local == j_glob) src_j = para_env%mepos
         CALL para_env%max(src_j)

         idx_j(:) = 0
         IF (para_env%mepos == src_j) idx_j(:) = atommax_j_idx(:, iat)
         CALL para_env%bcast(idx_j, src_j)

         IF (j_glob < 0.0_dp) THEN
            j_print = 0.0_dp
            idx_j(:) = 0
         ELSE
            j_print = SQRT(MAX(j_glob, 0.0_dp))
         END IF

         sym = "??"
         ik = particle_set(iat)%atomic_kind%kind_number
         sym = atomic_kind_set(ik)%element_symbol

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, &
                   '(T2,A,I0,1X,A,A,ES11.3,A,I0,A,I0,A,I0,A,ES11.3,A,I0,A,I0,A,I0,A)') &
               'RI-RS fit error | atommax ', iat, TRIM(sym), ' | e3c ', e3c_print, &
               ' (mu= ', idx_e3c(1), ' nu= ', idx_e3c(2), ' P= ', idx_e3c(3), ') | J ', &
               j_print, ' (mu= ', idx_j(1), ' nu= ', idx_j(2), ' P= ', idx_j(3), ')'
            WRITE (bs_env%unit_nr, &
                   '(T2,A,I0,1X,A,A,ES11.3,A,I0,A,I0,A,I0,A)') &
               'RI-RS fit error | atommaxrel ', iat, TRIM(sym), ' | e3c_rel ', e3crel_print, &
               ' (mu= ', idx_e3crel(1), ' nu= ', idx_e3crel(2), ' P= ', idx_e3crel(3), ')'
         END IF
      END DO

   END SUBROUTINE ri_rs_atommax_reduce_print

! **************************************************************************************************
!> \brief GW calculation using RI-RS formalism for molecules
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************

   SUBROUTINE gw_calc_non_periodic_ri_rs(qs_env, bs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_calc_non_periodic_ri_rs'

      INTEGER                                            :: handle
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma, fm_W_time
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CALL timeset(routineN, handle)

       !!========================================================================
       !! 0. Precompute AO and RI Radii
       !!    Determine per-atom cutoff radii from the most diffuse Gaussian
       !!    primitives in the AO and RI auxiliary basis sets.
       !!    Equations:
       !!    a. α_min,ao = min_{i,j} { ζ_ao(i,j) | ζ_ao(i,j) > 10⁻³ }
       !!    b. α_min,ri = min_{i,j} { ζ_ri(i,j) | ζ_ri(i,j) > 10⁻³ }
       !!    c. r_ao     = sqrt( -ln(ε) / α_min,ao )  (radius_ao_per_atom)
       !!    d. r_ri     = sqrt( -ln(ε) / α_min,ri )  (radius_ri_per_atom)
       !!========================================================================
      CALL precompute_ri_rs_radii(qs_env, bs_env)

       !!========================================================================
       !! 1. Grid Generation for RI-RS
       !!    (Modified Lebedev grids from Ivan Duchemin and Xavier Blase)
       !!    Generate flattened 1D array of grid points for RI-RS.
       !!    Equation: r_g(k) = R_A + r_g(A)
       !!========================================================================
      CALL ri_rs_grid_assembler(qs_env, bs_env, bs_env%ri_rs%grid_points)

       !!========================================================================
       !! 2. Atomic Basis Evaluation
       !!    Compute values of spherical atomic basis functions at grid points.
       !!    Expression: Φ_μl = Φ_μ(r_l) (mat_phi_mu_l)
       !!========================================================================
      CALL atomic_basis_at_grid_point(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                                      bs_env%ri_rs%mat_phi_mu_l)

       !!========================================================================
       !! 3. Compute RI-RS Coefficients (Z_lp)
       !!    Solve the regularized system for each atom P, where the grid domain
       !!    is restricted to r_l within a cutoff distance of atom P:
       !!    a. D_ll' = [ Σ_μ Φ_μ(r_l) Φ_μ(r_l') ]^2 (Equation 13)
       !!    b. D_lP  = Σ_{μν} Φ_μ(r_l) Φ_ν(r_l) (μν|P) (Equation 15)
       !!    c. Conditioning:
       !!       Dvec_l   = 1 / sqrt(D_ll)  (Diagonal scaling vector)
       !!       D'_ll' = Dvec_l * D_ll' * Dvec_l' + λδ_ll'
       !!       D'_lP  = Dvec_l * D_lP
       !!    d. Solve: Σ_l' D'_ll' * Z'_l'P = D'_lP (Equation 14)
       !!    e. Rescale: Z_lP = Z'_lP * Dvec_l   (Z_lP stored in mat_Z_lP)
       !!========================================================================
      CALL compute_coeff_Z_lP(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                              bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

       !!========================================================================
       !! 4. Compute Independent-Particle Polarizability (χ)
       !!    G^occ_µλ(i|τ|)  = sum_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
       !!    G^vir_µλ(i|τ|)  = sum_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
       !!    G^occ_ll'(i|τ|) = sum_µν Φ_µ(r_l) G^occ_µν Φ_ν(r_l')
       !!    G^vir_ll'(i|τ|) = sum_µν Φ_µ(r_l) G^vir_µν Φ_ν(r_l')
       !!    χ_ll'(iτ)       = G^occ_ll'(i|τ|) * G^vir_ll'(i|τ|)
       !!    χ_PQ(iτ)        = sum_ll' Z_lP χ_ll'(iτ) Z_l'Q
       !!========================================================================
      CALL get_mat_chi_Gamma_tau(bs_env, bs_env%mat_chi_Gamma_tau, &
                                 bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

       !!========================================================================
       !! 5. Compute Screened Interaction (W)
       !!    χ_PQ(iτ) -> χ_PQ(iω) -> ε_PQ(iω) -> W_PQ(iω) -> W_PQ(iτ)
       !!========================================================================
      CALL compute_W(bs_env, qs_env, bs_env%mat_chi_Gamma_tau, fm_W_time)

       !!========================================================================
       !! 6. Compute Exact Exchange Self-Energy (Σ^x)
       !!    D_µν          = sum_n^occ C_µn C_νn
       !!    D_ll'         = sum_µν Φ_µ(r_l) D_µν Φ_ν(r_l')
       !!    V^trunc_ll'   = sum_PQ Z_lP V^trunc_PQ Z_l'Q
       !!    Σ^x_ll'       = D_ll' * V^trunc_ll'
       !!    Σ^x_λσ(k=0)   = -sum_ll' Φ_λ(r_l) Σ^x_ll' Φ_σ(r_l')
       !!========================================================================
      CALL compute_Sigma_x(bs_env, qs_env, bs_env%ri_rs%mat_phi_mu_l, &
                           bs_env%ri_rs%mat_Z_lP, fm_Sigma_x_Gamma)

       !!========================================================================
       !! 7. Compute Correlation Self-Energy (Σ^c)
       !!    W_ll'(iτ)     =  sum_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
       !!    Σ^c_ll'(iτ)   = -G^occ_ll'(i|τ|) * W_ll'(iτ), for τ < 0
       !!    Σ^c_ll'(iτ)   =  G^vir_ll'(i|τ|) * W_ll'(iτ), for τ > 0
       !!    Σ^c_λσ(iτ)    =  sum_ll' Φ_λ(r_l) Σ^c_ll'(iτ) Φ_σ(r_l')
       !!========================================================================
      CALL compute_Sigma_c(bs_env, fm_W_time, bs_env%ri_rs%mat_phi_mu_l, &
                           bs_env%ri_rs%mat_Z_lP, fm_Sigma_c_Gamma_time)

       !!========================================================================
       !! 8. Compute Quasiparticle Energies
       !!    Σ^c_λσ(iτ) -> Σ^c_nn(ϵ)
       !!    ϵ_nk^GW = ϵ_nk^DFT + Σ^c_nn(ϵ) + Σ^x_nn - v^xc_nn
       !!========================================================================
      CALL compute_QP_energies(bs_env, qs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

      CALL de_init_bs_env(bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_non_periodic_ri_rs

! **************************************************************************************************
!> \brief Compute per-atom AO and RI basis radii from the most diffuse Gaussian
!>        primitive in the AO ("ORB") and RI auxiliary ("RI_AUX") basis sets.
!>        Stores results in bs_env%ri_rs%radius_ao_per_atom(:) and
!>        bs_env%ri_rs%radius_ri_per_atom(:) and prints a per-atom table.
!>        Radius:  r_kind = sqrt(-log(eps) / alpha_min_kind)
!>        with eps = eps_filter.
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************

   SUBROUTINE precompute_ri_rs_radii(qs_env, bs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'precompute_ri_rs_radii'

      INTEGER                                            :: handle, i, iatom, ikind, j, natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      REAL(KIND=dp)                                      :: eps
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: alpha_min_ao_kind, alpha_min_ri_kind
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet_ao, zet_ri
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, nkind=nkind, atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set)
      natom = SIZE(particle_set)

      eps = bs_env%eps_filter

      ALLOCATE (alpha_min_ao_kind(nkind), alpha_min_ri_kind(nkind))
      alpha_min_ao_kind = HUGE(1.0_dp)
      alpha_min_ri_kind = HUGE(1.0_dp)

      DO ikind = 1, nkind
         zet_ao => bs_env%basis_set_AO(ikind)%gto_basis_set%zet
         zet_ri => bs_env%basis_set_RI(ikind)%gto_basis_set%zet

         DO i = 1, SIZE(zet_ao, 1)
            DO j = 1, SIZE(zet_ao, 2)
               IF (zet_ao(i, j) > 1.0E-3_dp) &
                  alpha_min_ao_kind(ikind) = MIN(alpha_min_ao_kind(ikind), zet_ao(i, j))
            END DO
         END DO
         DO i = 1, SIZE(zet_ri, 1)
            DO j = 1, SIZE(zet_ri, 2)
               IF (zet_ri(i, j) > 1.0E-3_dp) &
                  alpha_min_ri_kind(ikind) = MIN(alpha_min_ri_kind(ikind), zet_ri(i, j))
            END DO
         END DO
      END DO

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)

      ALLOCATE (bs_env%ri_rs%radius_ao_per_atom(natom))
      ALLOCATE (bs_env%ri_rs%radius_ri_per_atom(natom))
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         bs_env%ri_rs%radius_ao_per_atom(iatom) = SQRT(-LOG(eps)/alpha_min_ao_kind(ikind))
         bs_env%ri_rs%radius_ri_per_atom(iatom) = SQRT(-LOG(eps)/alpha_min_ri_kind(ikind))
      END DO

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A)') 'Per-kind RI-RS basis radii (Bohr):'
         WRITE (bs_env%unit_nr, '(T4,A6,2X,A4,2A14)') 'Kind', 'Elem', 'r_AO', 'r_RI'
         DO ikind = 1, nkind
            WRITE (bs_env%unit_nr, '(T4,I6,2X,A4,2F14.4)') &
               ikind, &
               atomic_kind_set(ikind)%element_symbol, &
               SQRT(-LOG(eps)/alpha_min_ao_kind(ikind)), &
               SQRT(-LOG(eps)/alpha_min_ri_kind(ikind))
         END DO
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      DEALLOCATE (alpha_min_ao_kind, alpha_min_ri_kind, kind_of)

      CALL timestop(handle)

   END SUBROUTINE precompute_ri_rs_radii

! **************************************************************************************************
!> \brief Compute grid points for RI-RS
!>        Right now based on Ivan and Xavier implementation
!>        JCP 150, 174120 (2019), JCTC 17, 2383 (2021)
!> \param qs_env ...
!> \param bs_env ...
!> \param ri_rs_grid_points ...
! **************************************************************************************************

   SUBROUTINE ri_rs_grid_assembler(qs_env, bs_env, ri_rs_grid_points)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: ri_rs_grid_points(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ri_rs_grid_assembler'

      INTEGER                                            :: atom_idx, end_idx, handle, ikind, j, &
                                                            natom, nkind, start_idx, &
                                                            total_grid_npts
      INTEGER, ALLOCATABLE                               :: atom_to_kind(:), ri_rs_grid_offsets(:)
      REAL(KIND=dp)                                      :: atomic_center(3)
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

       !! Get the information about the atoms in the system
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, particle_set=particle_set)

      nkind = SIZE(atomic_kind_set)
      natom = SIZE(particle_set)

       !! 1. Build the grid cache
      CALL build_grid_cache(bs_env, atomic_kind_set)

       !! 2. Calculate grid counts and offsets
      ALLOCATE (ri_rs_grid_offsets(natom + 1))
      ALLOCATE (atom_to_kind(natom))

      total_grid_npts = 0
      DO ikind = 1, nkind
         DO j = 1, SIZE(atomic_kind_set(ikind)%atom_list)
            atom_idx = atomic_kind_set(ikind)%atom_list(j)
            atom_to_kind(atom_idx) = ikind
            ri_rs_grid_offsets(atom_idx) = total_grid_npts + 1
            total_grid_npts = total_grid_npts + bs_env%ri_rs%grid_cache(ikind)%npts
         END DO
      END DO

      ri_rs_grid_offsets(natom + 1) = total_grid_npts + 1

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, FMT="(T2,A,T69,I12)") &
            'Total grid points used for RI-RS:', total_grid_npts
         WRITE (bs_env%unit_nr, "(A)") ' '
      END IF

       !! 3. Allocate the global ri_rs_grid arrays
      ALLOCATE (ri_rs_grid_points(3, total_grid_npts))

       !! 4. Parallelize the grid generation loop
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(ri_rs_grid_points, ri_rs_grid_offsets, atom_to_kind, &
      !$OMP        particle_set, bs_env, natom) &
      !$OMP PRIVATE(atom_idx, ikind, atomic_center, start_idx, end_idx) &
      !$OMP SCHEDULE(DYNAMIC, 1)
      DO atom_idx = 1, natom
         ikind = atom_to_kind(atom_idx)
         atomic_center(:) = particle_set(atom_idx)%r(:)

         start_idx = ri_rs_grid_offsets(atom_idx)
         end_idx = start_idx + bs_env%ri_rs%grid_cache(ikind)%npts - 1

          !! Shift the cached origin grid by the atom's center
         ri_rs_grid_points(1, start_idx:end_idx) = bs_env%ri_rs%grid_cache(ikind)%raw_points(1, :) + atomic_center(1)
         ri_rs_grid_points(2, start_idx:end_idx) = bs_env%ri_rs%grid_cache(ikind)%raw_points(2, :) + atomic_center(2)
         ri_rs_grid_points(3, start_idx:end_idx) = bs_env%ri_rs%grid_cache(ikind)%raw_points(3, :) + atomic_center(3)

      END DO
      !$OMP END PARALLEL DO

       !! 5. Cleanup memory
      IF (ALLOCATED(bs_env%ri_rs%grid_cache)) THEN
         DO ikind = 1, nkind
            IF (ALLOCATED(bs_env%ri_rs%grid_cache(ikind)%raw_points)) DEALLOCATE (bs_env%ri_rs%grid_cache(ikind)%raw_points)
         END DO
         DEALLOCATE (bs_env%ri_rs%grid_cache)
      END IF

      DEALLOCATE (atom_to_kind, ri_rs_grid_offsets)
      CALL timestop(handle)

   END SUBROUTINE ri_rs_grid_assembler

! **************************************************************************************************
!> \brief Reads grids from .ion files and stores them in memory based on grid_select
!> \param bs_env ...
!> \param atomic_kind_set ...
! **************************************************************************************************

   SUBROUTINE build_grid_cache(bs_env, atomic_kind_set)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

      CHARACTER(LEN=default_path_length)                 :: filename, full_path, line
      CHARACTER(LEN=default_string_length)               :: atom_sym, kind_name, suffix
      INTEGER                                            :: colon_idx, i, ierr, ikind, irec, iunit, &
                                                            nkind
      LOGICAL, ALLOCATABLE                               :: record_used(:)
      REAL(KIND=dp)                                      :: pt(3)

      ! Determine the file suffix based on the user's choice
      IF (bs_env%ri_rs%grid_select == 1) THEN
         suffix = "_def2-tzvp-rs.ion"
      ELSE IF (bs_env%ri_rs%grid_select == 2) THEN
         suffix = "_cc-pvtz-rs.ion"
      ELSE IF (bs_env%ri_rs%grid_select == 3) THEN
         IF (LEN_TRIM(bs_env%ri_rs%grid_file_suffix) > 0) THEN
            suffix = TRIM(bs_env%ri_rs%grid_file_suffix)
         ELSE
            suffix = "_rirs-tzvp.ion"
         END IF
      ELSE
         CPABORT("Unknown grid_select (1=def2, 2=cc-pVTZ, 3=RI-RS).")
      END IF

      nkind = SIZE(atomic_kind_set)
      IF (.NOT. ALLOCATED(bs_env%ri_rs%grid_cache)) ALLOCATE (bs_env%ri_rs%grid_cache(nkind))

      IF (ALLOCATED(bs_env%ri_rs%escalation)) THEN
         ALLOCATE (record_used(SIZE(bs_env%ri_rs%escalation)))
         record_used(:) = .FALSE.
      END IF

      DO ikind = 1, nkind
         atom_sym = TRIM(atomic_kind_set(ikind)%element_symbol)
         filename = TRIM(atom_sym)//TRIM(suffix)

         full_path = "ri_rs_grid/"//TRIM(filename)

         CALL open_file(file_name=TRIM(full_path), unit_number=iunit, &
                        file_action="READ", file_status="OLD")

         ! Parse the preamble for 'n points'
         bs_env%ri_rs%grid_cache(ikind)%npts = 0
         DO
            READ (iunit, '(A)', IOSTAT=ierr) line
            IF (ierr /= 0) EXIT
            IF (INDEX(line, 'n points') > 0) THEN
               colon_idx = INDEX(line, ':')
               READ (line(colon_idx + 1:), *) bs_env%ri_rs%grid_cache(ikind)%npts
               EXIT
            END IF
         END DO

         ! Allocate the cache array for this specific element
         ALLOCATE (bs_env%ri_rs%grid_cache(ikind)%raw_points(3, bs_env%ri_rs%grid_cache(ikind)%npts))

         REWIND (iunit)
         DO
            READ (iunit, '(A)', IOSTAT=ierr) line
            IF (ierr /= 0) EXIT
            IF (INDEX(line, '<grid_points>') > 0) EXIT
         END DO

         ! Read the raw grid coordinates
         DO i = 1, bs_env%ri_rs%grid_cache(ikind)%npts
            READ (iunit, *, IOSTAT=ierr) pt(1), pt(2), pt(3)
            IF (ierr /= 0) THEN
               CPABORT("Unexpected EOF in grid file ")
            END IF
            bs_env%ri_rs%grid_cache(ikind)%raw_points(:, i) = pt(:)
         END DO

         ! A-posteriori per-kind angular escalation (&RI_RS_GRID_ESCALATION)
         IF (ALLOCATED(bs_env%ri_rs%escalation)) THEN
            kind_name = TRIM(atomic_kind_set(ikind)%name)
            CALL uppercase(kind_name)
            DO irec = 1, SIZE(bs_env%ri_rs%escalation)
               IF (TRIM(bs_env%ri_rs%escalation(irec)%kind_name) == TRIM(kind_name)) THEN
                  record_used(irec) = .TRUE.
                  IF (bs_env%ri_rs%escalation(irec)%steps > 0 .OR. &
                      ALLOCATED(bs_env%ri_rs%escalation(irec)%shell_orders)) THEN
                     CALL escalate_grid_cache_entry(bs_env%ri_rs%grid_cache(ikind), &
                                                    bs_env%ri_rs%escalation(irec), iunit, &
                                                    TRIM(full_path), TRIM(kind_name), &
                                                    bs_env%unit_nr)
                  END IF
               END IF
            END DO
         END IF

         CALL close_file(unit_number=iunit)
      END DO

      IF (ALLOCATED(record_used)) THEN
         DO irec = 1, SIZE(record_used)
            IF (.NOT. record_used(irec)) THEN
               CALL cp_abort(__LOCATION__, &
                             "&RI_RS_GRID_ESCALATION: kind "// &
                             TRIM(bs_env%ri_rs%escalation(irec)%kind_name)// &
                             " does not match any atomic kind of the system.")
            END IF
         END DO
         DEALLOCATE (record_used)
      END IF

   END SUBROUTINE build_grid_cache

! **************************************************************************************************
!> \brief Applies one &RI_RS_GRID_ESCALATION record to a cached grid: parses the <ri_rs_shells>
!>        metadata block from the open .ion file and rebuilds the selected shells at higher
!>        Lebedev orders on frozen shell radii. Aborts with a clear message when the file
!>        carries no shell metadata (arbitrary point clouds cannot be escalated).
!> \param cache      the per-kind grid cache entry (points replaced in place)
!> \param esc        the escalation record (steps, radial window)
!> \param iunit      open unit of the .ion file (will be rewound)
!> \param grid_path  file path, for error messages
!> \param kind_name  kind label, for messages
!> \param unit_nr    output unit (> 0 on the ionode only)
! **************************************************************************************************
   SUBROUTINE escalate_grid_cache_entry(cache, esc, iunit, grid_path, kind_name, unit_nr)

      TYPE(rirs_grid_type), INTENT(INOUT)                :: cache
      TYPE(rirs_grid_escalation_type), INTENT(IN)        :: esc
      INTEGER, INTENT(IN)                                :: iunit
      CHARACTER(LEN=*), INTENT(IN)                       :: grid_path, kind_name
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=32)                                  :: word
      CHARACTER(LEN=default_path_length)                 :: line
      INTEGER :: colon_idx, i, idx_base, idx_new, ierr, ipt, ish, ish_read, n_ladder, n_lowered, &
         n_raised, n_shells, new_npts, npts_center, opos
      INTEGER, ALLOCATABLE                               :: shell_first(:), shell_idx_new(:), &
                                                            shell_npts(:), shell_order(:)
      LOGICAL                                            :: found, selected
      REAL(KIND=dp)                                      :: rnorm, tol
      REAL(KIND=dp), ALLOCATABLE                         :: new_points(:, :), shell_radius(:)

      ! trigger the lazy initialization of the Lebedev tables
      idx_base = get_number_of_lebedev_grid(l=3)
      n_ladder = SIZE(lebedev_grid)

      ! ---- parse the <ri_rs_shells> metadata block ----
      REWIND (iunit)
      found = .FALSE.
      DO
         READ (iunit, '(A)', IOSTAT=ierr) line
         IF (ierr /= 0) EXIT
         IF (INDEX(line, '<ri_rs_shells') > 0) THEN
            found = .TRUE.
            EXIT
         END IF
      END DO
      IF (.NOT. found) THEN
         CALL cp_abort(__LOCATION__, &
                       "&RI_RS_GRID_ESCALATION requested for kind "//kind_name// &
                       ", but the grid file "//grid_path// &
                       " carries no <ri_rs_shells> metadata block. Only grids shipped with "// &
                       "explicit shell metadata can be escalated; arbitrary point clouds cannot.")
      END IF

      n_shells = 0
      npts_center = -1
      ish_read = 0
      DO
         READ (iunit, '(A)', IOSTAT=ierr) line
         IF (ierr /= 0) THEN
            CALL cp_abort(__LOCATION__, "Truncated <ri_rs_shells> block in "//grid_path)
         END IF
         IF (INDEX(line, '</ri_rs_shells>') > 0) EXIT
         IF (INDEX(line, 'center_points') > 0) THEN
            colon_idx = INDEX(line, ':')
            READ (line(colon_idx + 1:), *) npts_center
         ELSE IF (INDEX(line, 'n_shells') > 0) THEN
            colon_idx = INDEX(line, ':')
            READ (line(colon_idx + 1:), *) n_shells
            ALLOCATE (shell_radius(n_shells), shell_order(n_shells), shell_npts(n_shells), &
                      shell_first(n_shells), shell_idx_new(n_shells))
         ELSE
            READ (line, *, IOSTAT=ierr) word, ish, rnorm, i, ipt, opos
            IF (ierr /= 0) CYCLE
            IF (TRIM(word) /= 'shell') CYCLE
            IF (n_shells == 0) THEN
               CALL cp_abort(__LOCATION__, "<ri_rs_shells> block in "//grid_path// &
                             ": shell record before n_shells.")
            END IF
            ish_read = ish_read + 1
            IF (ish /= ish_read .OR. ish > n_shells) THEN
               CALL cp_abort(__LOCATION__, "<ri_rs_shells> block in "//grid_path// &
                             ": shell indices not consecutive.")
            END IF
            shell_radius(ish) = rnorm
            shell_order(ish) = i
            shell_npts(ish) = ipt
            shell_first(ish) = opos
         END IF
      END DO
      IF (npts_center < 0 .OR. n_shells == 0 .OR. ish_read /= n_shells) THEN
         CALL cp_abort(__LOCATION__, "<ri_rs_shells> block in "//grid_path// &
                       " is incomplete (center_points/n_shells/shell records).")
      END IF

      ! ---- validate the metadata against the point payload ----
      opos = npts_center + 1
      DO ish = 1, n_shells
         IF (shell_first(ish) /= opos) THEN
            CALL cp_abort(__LOCATION__, "<ri_rs_shells> block in "//grid_path// &
                          ": shell point ranges are not contiguous after the center points.")
         END IF
         opos = opos + shell_npts(ish)
         IF (shell_radius(ish) <= 0.0_dp) THEN
            CALL cp_abort(__LOCATION__, "<ri_rs_shells> block in "//grid_path// &
                          ": non-positive shell radius.")
         END IF
      END DO
      IF (opos - 1 /= cache%npts) THEN
         CALL cp_abort(__LOCATION__, "<ri_rs_shells> block in "//grid_path// &
                       ": center_points + shell points do not add up to 'n points'.")
      END IF
      IF (ALLOCATED(esc%shell_orders)) THEN
         IF (SIZE(esc%shell_orders) /= n_shells) THEN
            WRITE (line, '(A,I0,A,I0,A)') &
               "&RI_RS_GRID_ESCALATION SHELL_ORDERS has ", SIZE(esc%shell_orders), &
               " entries but the grid has ", n_shells, " shells."
            CALL cp_abort(__LOCATION__, TRIM(line)//" File: "//TRIM(grid_path))
         END IF
      END IF
      DO ish = 1, n_shells
         idx_base = 0
         DO i = 1, n_ladder
            IF (lebedev_grid(i)%l == shell_order(ish)) idx_base = i
         END DO
         IF (idx_base == 0) THEN
            CALL cp_abort(__LOCATION__, "<ri_rs_shells> block in "//grid_path// &
                          ": shell Lebedev order is not implemented in CP2K "// &
                          "(implemented: l = 3, 5, 7, 9, 11, 15, ..., 53).")
         END IF
         IF (lebedev_grid(idx_base)%n /= shell_npts(ish)) THEN
            CALL cp_abort(__LOCATION__, "<ri_rs_shells> block in "//grid_path// &
                          ": shell point count does not match its declared Lebedev order.")
         END IF
         tol = 1.0E-8_dp*MAX(1.0_dp, shell_radius(ish))
         DO ipt = shell_first(ish), shell_first(ish) + shell_npts(ish) - 1
            rnorm = SQRT(SUM(cache%raw_points(:, ipt)**2))
            IF (ABS(rnorm - shell_radius(ish)) > tol) THEN
               CALL cp_abort(__LOCATION__, "<ri_rs_shells> block in "//grid_path// &
                             ": point radii deviate from the declared shell radius.")
            END IF
         END DO
         ! resolve the escalated ladder position of this shell
         IF (ALLOCATED(esc%shell_orders)) THEN
            ! explicit per-shell order list overrides STEPS and the radial window
            idx_new = 0
            DO i = 1, n_ladder
               IF (lebedev_grid(i)%l == esc%shell_orders(ish)) idx_new = i
            END DO
            IF (idx_new == 0) THEN
               WRITE (line, '(A,I0,A,I0,A)') &
                  "&RI_RS_GRID_ESCALATION SHELL_ORDERS entry ", ish, " requests Lebedev order ", &
                  esc%shell_orders(ish), " which is not implemented (l = 3, 5, 7, 9, 11, 15, ...)."
               CALL cp_abort(__LOCATION__, TRIM(line))
            END IF
            shell_idx_new(ish) = idx_new
         ELSE
            selected = (shell_radius(ish) >= esc%r_min) .AND. &
                       (esc%r_max < 0.0_dp .OR. shell_radius(ish) <= esc%r_max)
            IF (selected) THEN
               shell_idx_new(ish) = MIN(idx_base + esc%steps, n_ladder)
            ELSE
               shell_idx_new(ish) = idx_base
            END IF
         END IF
      END DO

      ! ---- rebuild the point payload with escalated shells (frozen radii) ----
      n_raised = 0
      n_lowered = 0
      new_npts = npts_center
      DO ish = 1, n_shells
         new_npts = new_npts + lebedev_grid(shell_idx_new(ish))%n
      END DO
      ALLOCATE (new_points(3, new_npts))
      new_points(:, 1:npts_center) = cache%raw_points(:, 1:npts_center)
      opos = npts_center
      DO ish = 1, n_shells
         idx_base = 0
         DO i = 1, n_ladder
            IF (lebedev_grid(i)%l == shell_order(ish)) idx_base = i
         END DO
         idx_new = shell_idx_new(ish)
         IF (idx_new /= idx_base) THEN
            IF (idx_new > idx_base) THEN
               n_raised = n_raised + 1
            ELSE
               n_lowered = n_lowered + 1
            END IF
            DO ipt = 1, lebedev_grid(idx_new)%n
               new_points(:, opos + ipt) = shell_radius(ish)*lebedev_grid(idx_new)%r(:, ipt)
            END DO
         ELSE
            new_points(:, opos + 1:opos + shell_npts(ish)) = &
               cache%raw_points(:, shell_first(ish):shell_first(ish) + shell_npts(ish) - 1)
         END IF
         opos = opos + lebedev_grid(idx_new)%n
      END DO

      IF (unit_nr > 0) THEN
         IF (ALLOCATED(esc%shell_orders)) THEN
            WRITE (unit_nr, FMT="(T2,A,A,A,I4,A)") &
               'RI-RS grid escalation for kind ', TRIM(kind_name), &
               ': explicit SHELL_ORDERS over ', n_shells, ' shells'
         ELSE
            WRITE (unit_nr, FMT="(T2,A,A,A,I3,A)") &
               'RI-RS grid escalation for kind ', TRIM(kind_name), ':', esc%steps, ' Lebedev step(s)'
         END IF
         WRITE (unit_nr, FMT="(T2,A,I4,A,I4,A)") &
            'shells raised:', n_raised, '  lowered:', n_lowered, '  (radii frozen)'
         WRITE (unit_nr, FMT="(T2,A,I8,A,I8)") &
            'grid points per atom of this kind:', cache%npts, ' ->', new_npts
         WRITE (unit_nr, "(A)") ' '
      END IF

      CALL MOVE_ALLOC(new_points, cache%raw_points)
      cache%npts = new_npts

      DEALLOCATE (shell_radius, shell_order, shell_npts, shell_first, shell_idx_new)

   END SUBROUTINE escalate_grid_cache_entry

! **************************************************************************************************
!> \brief Evaluates atomic basis functions on a real-space grid and builds a sparse DBCSR matrix.
!> \param qs_env ...
!> \param bs_env ...
!> \param ri_rs_grid_points ...
!> \param mat_phi_mu_l ...
! **************************************************************************************************

   SUBROUTINE atomic_basis_at_grid_point(qs_env, bs_env, ri_rs_grid_points, mat_phi_mu_l)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT)          :: ri_rs_grid_points(:, :)
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_phi_mu_l

      CHARACTER(LEN=*), PARAMETER :: routineN = 'atomic_basis_at_grid_point'

      INTEGER :: c_size, chunk_size, dimen_ORB, handle, i, i_blk, iatom, natom, npcol, nprow, &
         num_grid_chunks, r_end, r_start, total_grid_npts
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf
      INTEGER, DIMENSION(:), POINTER                     :: c_blk_sizes, col_dist, col_dist_ks, &
                                                            r_blk_sizes, row_dist, row_dist_ks
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atom_col_buffer
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist_ks
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      chunk_size = max_elements_per_block

      ! Extract environment variables
      CALL get_qs_env(qs_env, cell=cell, atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, particle_set=particle_set, &
                      para_env=para_env)

      natom = SIZE(particle_set)
      total_grid_npts = SIZE(ri_rs_grid_points, 2)

      ! Map the starting indices of spherical gaussian functions (SGF) for each atom
      ALLOCATE (first_sgf(natom + 1))
      CALL get_basis_offsets(particle_set, qs_kind_set, first_sgf, dimen_ORB)

      ! =========================================================================
      ! 1. SETUP DBCSR MATRIX TOPOLOGY
      ! =========================================================================

      ! A. Define Column Block Sizes (1 Block = 1 Atom's full basis set)
      ALLOCATE (c_blk_sizes(natom))
      DO iatom = 1, natom
         c_blk_sizes(iatom) = first_sgf(iatom + 1) - first_sgf(iatom)
      END DO

      ! B. Define Row Block Sizes (Grid chunks of max size 256)
      num_grid_chunks = CEILING(REAL(total_grid_npts, KIND=dp)/REAL(chunk_size, KIND=dp))
      ALLOCATE (r_blk_sizes(num_grid_chunks))
      r_blk_sizes = chunk_size
      IF (MOD(total_grid_npts, chunk_size) /= 0) THEN
         r_blk_sizes(num_grid_chunks) = MOD(total_grid_npts, chunk_size)
      END IF

      ! C. Fetch CP2K's Default Process Grid Configuration
      CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist_ks)
      CALL dbcsr_distribution_get(dbcsr_dist_ks, row_dist=row_dist_ks, col_dist=col_dist_ks)

      ! D. Build Custom Mappings using Round-Robin across the 2D process grid
      !    (MAXVAL + 1 accounts for 0-based process indexing in DBCSR)
      nprow = MAXVAL(row_dist_ks) + 1
      npcol = MAXVAL(col_dist_ks) + 1

      ALLOCATE (row_dist(num_grid_chunks))
      DO i = 1, num_grid_chunks
         row_dist(i) = MOD(i - 1, nprow)
      END DO

      ALLOCATE (col_dist(natom))
      DO i = 1, natom
         col_dist(i) = MOD(i - 1, npcol)
      END DO

      ! E. Create the DBCSR Distribution and Initialize the Matrix
      CALL dbcsr_distribution_new(dist, template=dbcsr_dist_ks, &
                                  row_dist=row_dist, col_dist=col_dist)

      CALL dbcsr_create(mat_phi_mu_l, name="phi_val_sparse", dist=dist, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=r_blk_sizes, col_blk_size=c_blk_sizes)

      ! =========================================================================
      ! 2. STREAM DATA DIRECTLY INTO SPARSE MATRIX
      ! =========================================================================
      ! Iterate over the atoms assigned to this specific MPI rank
      DO iatom = para_env%mepos + 1, natom, para_env%num_pe

         c_size = c_blk_sizes(iatom)

         ! Allocate a temporary dense buffer just for this specific atom
         ALLOCATE (atom_col_buffer(total_grid_npts, c_size))
         atom_col_buffer = 0.0_dp

         ! Evaluate the basis functions on the grid. Skip grid points outside
         ! the spatial extent of the most diffuse AO Gaussian on iatom; beyond
         ! that radius the contribution is guaranteed below eps_filter.
         CALL fill_phi_for_atom(atom_col_buffer, ri_rs_grid_points, total_grid_npts, &
                                iatom, particle_set, qs_kind_set, cell, &
                                bs_env%ri_rs%radius_ao_per_atom(iatom)**2)

         ! Slice the dense column into chunks and insert into DBCSR
         DO i_blk = 1, num_grid_chunks
            r_start = (i_blk - 1)*chunk_size + 1
            r_end = MIN(i_blk*chunk_size, total_grid_npts)

            ! Apply dynamic sparsity filtering: Only store blocks with physical significance
            IF (MAXVAL(ABS(atom_col_buffer(r_start:r_end, 1:c_size))) > bs_env%eps_filter) THEN
               CALL dbcsr_put_block(mat_phi_mu_l, row=i_blk, col=iatom, &
                                    block=atom_col_buffer(r_start:r_end, 1:c_size))
            END IF
         END DO

         DEALLOCATE (atom_col_buffer)

      END DO

      ! Finalize triggers internal MPI communication to route blocks to their correct 2D process owners
      CALL dbcsr_finalize(mat_phi_mu_l)

      ! -------------------------------------------------------------------------
      ! CLEANUP
      ! -------------------------------------------------------------------------
      DEALLOCATE (first_sgf, r_blk_sizes, c_blk_sizes, row_dist, col_dist)
      CALL dbcsr_distribution_release(dist)

      CALL timestop(handle)

   END SUBROUTINE atomic_basis_at_grid_point

! **************************************************************************************************
!> \brief Compute value of all basis function for single atom across all grid points
!> \param phi_val ...
!> \param ri_rs_grid ...
!> \param npts ...
!> \param iatom ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param r2_threshold ...
! **************************************************************************************************

   SUBROUTINE fill_phi_for_atom(phi_val, ri_rs_grid, npts, iatom, &
                                particle_set, qs_kind_set, cell, r2_threshold)

      REAL(KIND=dp), INTENT(INOUT)                       :: phi_val(:, :)
      INTEGER, INTENT(IN)                                :: npts
      REAL(KIND=dp), INTENT(IN)                          :: ri_rs_grid(3, npts)
      INTEGER, INTENT(IN)                                :: iatom
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(IN)                          :: r2_threshold

      INTEGER                                            :: first_sgf, i_pt, ico, iend_co, ikind, &
                                                            ipgf, iset, isgf, ishell, istart_co, &
                                                            l, last_sgf, lx, ly, lz, n_cart_total, &
                                                            row_idx
      REAL(KIND=dp)                                      :: alpha, dist_vec(3), exp_val, poly, r2, &
                                                            r_atom(3), weight
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      ! Get Atom Info
      ikind = particle_set(iatom)%atomic_kind%kind_number
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, basis_type="ORB")
      IF (.NOT. ASSOCIATED(orb_basis_set)) RETURN

      r_atom = particle_set(iatom)%r

      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(phi_val, ri_rs_grid, npts, orb_basis_set, r_atom, cell, ncoset, indco, &
      !$OMP        r2_threshold) &
      !$OMP PRIVATE(i_pt, dist_vec, r2, iset, n_cart_total, ishell, l, istart_co, iend_co, &
      !$OMP         first_sgf, last_sgf, ipgf, alpha, exp_val, isgf, ico, row_idx, weight, &
      !$OMP         lx, ly, lz, poly) &
      !$OMP SCHEDULE(STATIC)
      DO i_pt = 1, npts
         dist_vec = pbc(ri_rs_grid(:, i_pt) - r_atom, cell)
         r2 = DOT_PRODUCT(dist_vec, dist_vec)
         IF (r2 > r2_threshold) CYCLE

         DO iset = 1, orb_basis_set%nset
            n_cart_total = ncoset(orb_basis_set%lmax(iset))

            DO ishell = 1, orb_basis_set%nshell(iset)
               l = orb_basis_set%l(ishell, iset)
               istart_co = ncoset(l - 1) + 1
               iend_co = ncoset(l)

               first_sgf = orb_basis_set%first_sgf(ishell, iset)
               last_sgf = orb_basis_set%last_sgf(ishell, iset)

               DO ipgf = 1, orb_basis_set%npgf(iset)
                  alpha = orb_basis_set%zet(ipgf, iset)
                  exp_val = EXP(-alpha*r2)

                  DO isgf = first_sgf, last_sgf
                     DO ico = istart_co, iend_co
                        row_idx = (ipgf - 1)*n_cart_total + ico
                        weight = orb_basis_set%sphi(row_idx, isgf)
                        lx = indco(1, ico)
                        ly = indco(2, ico)
                        lz = indco(3, ico)
                        poly = (dist_vec(1)**lx)*(dist_vec(2)**ly)*(dist_vec(3)**lz)

                        phi_val(i_pt, isgf) = phi_val(i_pt, isgf) + (weight*poly*exp_val)

                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO

   END SUBROUTINE fill_phi_for_atom

! **************************************************************************************************
!> \brief Helper for OMP threads to fill phi_val column values
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param first_sgf ...
!> \param total_sgf ...
! **************************************************************************************************

   SUBROUTINE get_basis_offsets(particle_set, qs_kind_set, first_sgf, total_sgf)

      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, INTENT(OUT)                               :: first_sgf(:), total_sgf

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_basis_offsets'

      INTEGER                                            :: handle, iatom, ikind, nsgf

      CALL timeset(routineN, handle)

      total_sgf = 0
      DO iatom = 1, SIZE(particle_set)
         first_sgf(iatom) = total_sgf + 1
         ikind = particle_set(iatom)%atomic_kind%kind_number
         CALL get_qs_kind(qs_kind_set(ikind), nsgf=nsgf, basis_type="ORB")
         total_sgf = total_sgf + nsgf
      END DO
      first_sgf(SIZE(particle_set) + 1) = total_sgf + 1

      CALL timestop(handle)

   END SUBROUTINE get_basis_offsets

! **************************************************************************************************
!> \brief Compute RI-RS Coefficients (Z_lP)
!> \param qs_env ...
!> \param bs_env ...
!> \param ri_rs_grid_points ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
! **************************************************************************************************

   SUBROUTINE compute_coeff_Z_lP(qs_env, bs_env, ri_rs_grid_points, mat_phi_mu_l, mat_Z_lP)

      ! Arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT)          :: ri_rs_grid_points(:, :)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_Z_lP

      CHARACTER(LEN=*), PARAMETER :: key = 'PROPERTIES%BANDSTRUCTURE%GW%PRINT%RESTART', &
         routineN = 'compute_coeff_Z_lP'

      INTEGER :: atom_j_mepos, atom_j_stride, atom_P, atom_P_start, atom_P_stride, col_end, &
         col_start, current_chunk_size, g, group_handle, handle, i, i_blk, ikind, info, j, j_ri, &
         l, loc_idx, loc_ptr, max_ao_size, max_loc_ri, my_group, n_ao_total, n_grid_total, &
         n_groups, n_loc_ri, n_local_grid, n_procs_per_atom, natom, nkind, num_grid_chunks, &
         P_loop_atom, r_end, r_start, source_atom
      INTEGER                                            :: j_full_info, n_ao_pair_total, naux_total
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: atommax_e3c_idx, atommax_e3crel_idx, &
                                                            atommax_j_idx
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: local_grid_idx, row_offset
      INTEGER, DIMENSION(:), POINTER                     :: col_dist_phi, col_dist_ri, r_blk_sizes, &
                                                            ri_blk_sizes, row_dist_grid
      REAL(KIND=dp)                                      :: cutoff_ri, d_sP, dist, r2_threshold, &
                                                            r_c, t1
      REAL(KIND=dp)                                      :: e3c_abs_max, e3c_vmax, z_abs_max
      INTEGER, DIMENSION(3)                              :: e3c_atom_max, e3c_idx_max
      REAL(KIND=dp)                                      :: j_sq_sum, j_pair_max
      INTEGER, DIMENSION(3)                              :: j_atom_max, j_idx_max
      LOGICAL                                            :: do_j_full, j_full_ok
      REAL(KIND=dp)                                      :: j_full_sq_sum, j_full_pair_max
      INTEGER, DIMENSION(3)                              :: j_full_atom_max, j_full_idx_max
      TYPE(ri_rs_topk_type)                              :: fit_topk
      REAL(KIND=dp), ALLOCATABLE                         :: cutoff_ri_per_kind(:)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: atommax_e3c, atommax_e3crel, atommax_j, &
                                                            cutoff_ri_per_atom, d_vec_local
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: D_local, d_lp_local, phi_local, &
                                                            sphere_grid, Z_blk
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: R_full
      REAL(KIND=dp), DIMENSION(3)                        :: pos_P
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_b, fm_struct_D
      TYPE(cp_fm_type)                                   :: fm_b, fm_D
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type)                      :: dist_phi, dist_Z
      TYPE(gw_3c_ctx_type)                               :: ctx_3c
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      CALL get_qs_env(qs_env, para_env=para_env, particle_set=particle_set, input=input, &
                      qs_kind_set=qs_kind_set, cell=cell, atomic_kind_set=atomic_kind_set)

      ! ---------------------------------------------------------------------
      ! Subgroup setup. Default G=1 keeps the single-rank BLAS path; G>1 splits
      ! ranks into atom-groups so the Cholesky on D_local distributes across G
      ! ranks (memory ~1/G) and the compute_d_lp build also splits across the
      ! subgroup. G=1 leaves para_env_sub / blacs_env_sub NULL — no subgroup
      ! comms created, atom_P loop uses per-rank round-robin, compute_d_lp runs
      ! its full atom_j range on each rank, no allreduce.
      ! ---------------------------------------------------------------------
      n_procs_per_atom = MIN(bs_env%ri_rs%n_procs_per_atom_z_lp, para_env%num_pe)
      IF (n_procs_per_atom < 1) n_procs_per_atom = 1

      NULLIFY (para_env_sub, blacs_env_sub)
      IF (n_procs_per_atom > 1) THEN
         n_groups = para_env%num_pe/n_procs_per_atom
         my_group = MIN(para_env%mepos/n_procs_per_atom, n_groups - 1)
         ALLOCATE (para_env_sub)
         CALL para_env_sub%from_split(para_env, my_group)
         CALL cp_blacs_env_create(blacs_env=blacs_env_sub, para_env=para_env_sub)
         atom_P_start = my_group + 1
         atom_P_stride = n_groups
         atom_j_mepos = para_env_sub%mepos
         atom_j_stride = para_env_sub%num_pe
      ELSE
         atom_P_start = para_env%mepos + 1
         atom_P_stride = para_env%num_pe
         atom_j_mepos = 0
         atom_j_stride = 1
      END IF

      natom = SIZE(bs_env%i_RI_start_from_atom)
      n_ao_total = bs_env%i_ao_end_from_atom(natom)
      naux_total = bs_env%i_RI_end_from_atom(natom)
      n_ao_pair_total = n_ao_total*n_ao_total
      n_grid_total = SIZE(ri_rs_grid_points, 2)

      ! =========================================================================
      ! 1. SETUP DBCSR TOPOLOGY & EXACT OFFSETS
      ! =========================================================================
      CALL dbcsr_get_info(mat_phi_mu_l, row_blk_size=r_blk_sizes, distribution=dist_phi)
      CALL dbcsr_distribution_get(dist_phi, row_dist=row_dist_grid, col_dist=col_dist_phi, group=group_handle)

      num_grid_chunks = SIZE(r_blk_sizes)

      ALLOCATE (row_offset(num_grid_chunks))
      row_offset(1) = 0
      DO i_blk = 2, num_grid_chunks
         row_offset(i_blk) = row_offset(i_blk - 1) + r_blk_sizes(i_blk - 1)
      END DO

      ALLOCATE (ri_blk_sizes(natom), col_dist_ri(natom))
      DO atom_P = 1, natom
         ri_blk_sizes(atom_P) = bs_env%i_RI_end_from_atom(atom_P) - bs_env%i_RI_start_from_atom(atom_P) + 1
         col_dist_ri(atom_P) = MOD(atom_P - 1, MAXVAL(col_dist_phi) + 1)
      END DO

      CALL dbcsr_distribution_new(dist_Z, template=dist_phi, row_dist=row_dist_grid, col_dist=col_dist_ri)

      IF (bs_env%ri_rs%Z_lP_exists) THEN
         CALL dbcsr_binary_read(filepath=TRIM(bs_env%prefix)//"Z_lP.matrix", &
                                distribution=dist_Z, &
                                matrix_new=mat_Z_lP)
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T57,A,F7.1,A)') &
               'Read Z_lP from file ', ' Execution time', m_walltime() - t1, ' s'
            WRITE (bs_env%unit_nr, '(A)') ' '
         END IF
      ELSE

         CALL dbcsr_create(mat_Z_lP, name="mat_Z_lP", dist=dist_Z, &
                           matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=r_blk_sizes, col_blk_size=ri_blk_sizes)

         max_ao_size = 0
         DO j = 1, SIZE(bs_env%i_ao_start_from_atom)
            max_ao_size = MAX(max_ao_size, bs_env%i_ao_end_from_atom(j) - bs_env%i_ao_start_from_atom(j) + 1)
         END DO
         max_loc_ri = MAXVAL(ri_blk_sizes)

         ! Per-atom RI-RS integration sphere:
         !   cutoff_ri(P) = r_c + r_AO(P)
         ! where r_c is the truncated-Coulomb cutoff of the RI metric. The
         ! CUTOFF_RADIUS_RI_RS keyword (when > 0) overrides the entire cutoff calculation.
         nkind = SIZE(atomic_kind_set)
         ALLOCATE (cutoff_ri_per_atom(natom))

         IF (bs_env%ri_rs%cutoff_radius_ri_rs > 0.0_dp) THEN
            cutoff_ri_per_atom(:) = bs_env%ri_rs%cutoff_radius_ri_rs
         ELSE
            r_c = bs_env%ri_metric%cutoff_radius
            DO P_loop_atom = 1, natom
               cutoff_ri_per_atom(P_loop_atom) = &
                  r_c + bs_env%ri_rs%radius_ao_per_atom(P_loop_atom)
            END DO
         END IF

         ALLOCATE (cutoff_ri_per_kind(nkind))
         cutoff_ri_per_kind(:) = 0.0_dp
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A)') 'Per-kind maximum RI-RS sphere cutoff (Bohr):'
            WRITE (bs_env%unit_nr, '(T4,A4,A14)') 'Kind', 'max cutoff_ri'
            ! Walk atoms once; accumulate per-kind max using atomic_kind index
            DO P_loop_atom = 1, natom
               ikind = particle_set(P_loop_atom)%atomic_kind%kind_number
               cutoff_ri_per_kind(ikind) = MAX(cutoff_ri_per_kind(ikind), &
                                               cutoff_ri_per_atom(P_loop_atom))
            END DO
            DO ikind = 1, nkind
               WRITE (bs_env%unit_nr, '(T4,A4,F14.4)') &
                  atomic_kind_set(ikind)%element_symbol, &
                  cutoff_ri_per_kind(ikind)
            END DO
            WRITE (bs_env%unit_nr, '(A)') ' '
         END IF
         DEALLOCATE (cutoff_ri_per_kind)

         ! Build the shared 3c-integral context once, before the atom_P loop and outside
         ! any OMP region. This hoists out the libint / t_c_g0 / md_ftable inits and the
         ! contracted sphi tables, which would otherwise re-run on every per-block call.
         ! gw_3c_ctx_create is MPI-collective (truncated-Coulomb init reads a file + bcasts).
         CALL gw_3c_ctx_create(ctx_3c, qs_env, bs_env%ri_metric, &
                               basis_j=bs_env%basis_set_AO, basis_k=bs_env%basis_set_AO, &
                               basis_i=bs_env%basis_set_RI)

         ! =========================================================================
         ! 2. MPI LOOP OVER ATOMS (Fully independent, no MPI barriers inside)
         ! This distributes the N_atoms among the available MPI processors.
         ! phi_local for each atom_P's cutoff sphere is built on the fly via
         ! fill_phi_for_atom — no dense replicated phi_global, no allreduce.
         ! =========================================================================
         ! RI-RS fit-error accumulators (filled per atom below, reduced + printed after
         ! the loop). Always initialised; only consumed when PRINT_FIT_ERROR is set.
         e3c_abs_max = -1.0_dp
         e3c_vmax = 0.0_dp
         z_abs_max = 0.0_dp
         e3c_idx_max(:) = 0
         e3c_atom_max(:) = 0
         j_sq_sum = 0.0_dp
         j_pair_max = -1.0_dp
         j_idx_max(:) = 0
         j_atom_max(:) = 0
         j_full_ok = .FALSE.
         j_full_info = 0
         j_full_sq_sum = 0.0_dp
         j_full_pair_max = -1.0_dp
         j_full_idx_max(:) = 0
         j_full_atom_max(:) = 0

         do_j_full = bs_env%ri_rs%print_fit_error .AND. &
                     (REAL(n_ao_pair_total, dp)*REAL(naux_total, dp)*8.0_dp <= ri_rs_j_full_max_bytes)
         IF (bs_env%ri_rs%print_fit_error) THEN
            ALLOCATE (atommax_e3c(natom), atommax_e3crel(natom), atommax_j(natom))
            ALLOCATE (atommax_e3c_idx(3, natom), atommax_e3crel_idx(3, natom), atommax_j_idx(3, natom))
            atommax_e3c(:) = -1.0_dp
            atommax_e3crel(:) = -1.0_dp
            atommax_j(:) = -1.0_dp
            atommax_e3c_idx(:, :) = 0
            atommax_e3crel_idx(:, :) = 0
            atommax_j_idx(:, :) = 0
            j_full_ok = .FALSE.
            j_full_info = 0
            IF (do_j_full) THEN
               ALLOCATE (R_full(n_ao_pair_total, naux_total))
               R_full(:, :) = 0.0_dp
            ELSE
               ALLOCATE (R_full(1, 1))
               IF (bs_env%unit_nr > 0) THEN
                  WRITE (bs_env%unit_nr, '(T2,A,F0.1,A)') &
                     'RI-RS fit error | J_full skipped: replicated R_full ', &
                     REAL(n_ao_pair_total, dp)*REAL(naux_total, dp)*8.0_dp/1.0E9_dp, &
                     ' GB/rank > cap; e3c + top-30 still computed.'
               END IF
            END IF
         END IF

         DO atom_P = atom_P_start, natom, atom_P_stride

            n_loc_ri = ri_blk_sizes(atom_P)
            pos_P(:) = particle_set(atom_P)%r(:)

            cutoff_ri = cutoff_ri_per_atom(atom_P)

            ! ---------------------------------------------------------------------
            ! A. Determine Local Grid Domain based on cutoff_ri
            ! ---------------------------------------------------------------------
            n_local_grid = 0
            DO l = 1, n_grid_total
               dist = SQRT(SUM((ri_rs_grid_points(1:3, l) - pos_P(1:3))**2))
               IF (dist <= cutoff_ri) n_local_grid = n_local_grid + 1
            END DO

            ALLOCATE (local_grid_idx(n_local_grid))

            n_local_grid = 0
            DO l = 1, n_grid_total
               dist = SQRT(SUM((ri_rs_grid_points(1:3, l) - pos_P(1:3))**2))
               IF (dist <= cutoff_ri) THEN
                  n_local_grid = n_local_grid + 1
                  local_grid_idx(n_local_grid) = l
               END IF
            END DO

            ! ---------------------------------------------------------------------
            ! B. Build phi_local on the fly via fill_phi_for_atom (no replication)
            ! ---------------------------------------------------------------------
            ! Copy the local sphere coordinates into a contiguous buffer for
            ! fill_phi_for_atom. Only source atoms whose AO basis can reach the
            ! sphere of atom_P contribute; the others are pruned by the inequality
            ! using the per-atom radii precomputed in precompute_ri_rs_radii.
            ALLOCATE (sphere_grid(3, n_local_grid))
            DO loc_idx = 1, n_local_grid
               sphere_grid(:, loc_idx) = ri_rs_grid_points(:, local_grid_idx(loc_idx))
            END DO

            ALLOCATE (phi_local(n_local_grid, n_ao_total))
            phi_local = 0.0_dp

            DO source_atom = 1, natom
               d_sP = NORM2(particle_set(source_atom)%r(:) - pos_P(:))
               IF (d_sP > bs_env%ri_rs%radius_ao_per_atom(source_atom) + cutoff_ri) CYCLE

               col_start = bs_env%i_ao_start_from_atom(source_atom)
               col_end = bs_env%i_ao_end_from_atom(source_atom)
               r2_threshold = bs_env%ri_rs%radius_ao_per_atom(source_atom)**2

               CALL fill_phi_for_atom(phi_local(:, col_start:col_end), sphere_grid, &
                                      n_local_grid, source_atom, particle_set, qs_kind_set, &
                                      cell, r2_threshold)
            END DO

            DEALLOCATE (sphere_grid)

            ! ---------------------------------------------------------------------
            ! C. Build Local RHS Matrix (d_lp_local)
            !    Done first so the subgroup-distributed compute_d_lp + allreduce
            !    is not entangled with the LHS build. compute_d_lp does not
            !    depend on D_local or d_vec_local.
            ! ---------------------------------------------------------------------
            ALLOCATE (d_lp_local(n_local_grid, n_loc_ri))
            d_lp_local = 0.0_dp

            CALL compute_d_lp(bs_env, ctx_3c, phi_local, d_lp_local, n_local_grid, &
                              n_loc_ri, atom_P, max_ao_size, atom_j_mepos, atom_j_stride)

            ! Reduce per-subgroup-rank partials into the replicated d_lp_local.
            ! Skipped for G=1 (BLAS path): each rank has the full sum locally.
            IF (n_procs_per_atom > 1) THEN
               CALL para_env_sub%sum(d_lp_local)
            END IF

            ! ---------------------------------------------------------------------
            ! D. Build d_vec_local (Jacobi diagonal) + LHS — BLAS or ScaLAPACK
            ! ---------------------------------------------------------------------
            ALLOCATE (d_vec_local(n_local_grid))

            IF (n_procs_per_atom == 1) THEN
               ! BLAS path: build D_local densely, compute d_vec as side-effect
               ! of the Jacobi step (preserves bit-identical arithmetic with the
               ! previous branch).
               ALLOCATE (D_local(n_local_grid, n_local_grid))
               D_local = 0.0_dp

               CALL dsyrk("L", "N", n_local_grid, n_ao_total, 1.0_dp, phi_local, &
                          n_local_grid, 0.0_dp, D_local, n_local_grid)

               !$OMP PARALLEL DO DEFAULT(NONE) &
               !$OMP SHARED(n_local_grid, D_local, d_vec_local, bs_env) &
               !$OMP PRIVATE(i) &
               !$OMP SCHEDULE(STATIC)
               DO i = 1, n_local_grid
                  D_local(i, i) = D_local(i, i)**2
                  d_vec_local(i) = 1.0_dp/SQRT(MAX(D_local(i, i), 1.0E-16_dp))
                  D_local(i, i) = (D_local(i, i)*d_vec_local(i)**2) + bs_env%ri_rs%tikhonov
               END DO
               !$OMP END PARALLEL DO

               !$OMP PARALLEL DO DEFAULT(NONE) &
               !$OMP SHARED(n_local_grid, D_local, d_vec_local) &
               !$OMP PRIVATE(j, i) &
               !$OMP SCHEDULE(DYNAMIC)
               DO j = 1, n_local_grid
                  DO i = j + 1, n_local_grid
                     D_local(i, j) = D_local(i, j)**2
                     D_local(i, j) = D_local(i, j)*d_vec_local(i)*d_vec_local(j)
                     D_local(j, i) = D_local(i, j)
                  END DO
               END DO
               !$OMP END PARALLEL DO
            ELSE
               ! ScaLAPACK path: d_vec computed directly from phi (= 1/||phi_i||^2);
               ! solve_D_lp_distributed builds D block-cyclic internally with
               ! the squared+scaled values, so no dense D_local on this rank.
               !$OMP PARALLEL DO DEFAULT(NONE) &
               !$OMP SHARED(n_local_grid, n_ao_total, phi_local, d_vec_local) &
               !$OMP PRIVATE(i, j) &
               !$OMP SCHEDULE(STATIC)
               DO i = 1, n_local_grid
                  d_vec_local(i) = 0.0_dp
                  DO j = 1, n_ao_total
                     d_vec_local(i) = d_vec_local(i) + phi_local(i, j)*phi_local(i, j)
                  END DO
                  d_vec_local(i) = 1.0_dp/MAX(d_vec_local(i), 1.0E-16_dp)
               END DO
               !$OMP END PARALLEL DO
            END IF

            ! ---------------------------------------------------------------------
            ! E. Pre-scale d_lp by d_vec
            ! ---------------------------------------------------------------------
            !$OMP PARALLEL DO DEFAULT(NONE) &
            !$OMP SHARED(n_loc_ri, n_local_grid, d_lp_local, d_vec_local) &
            !$OMP PRIVATE(j_ri, i) &
            !$OMP SCHEDULE(STATIC)
            DO j_ri = 1, n_loc_ri
               DO i = 1, n_local_grid
                  d_lp_local(i, j_ri) = d_lp_local(i, j_ri)*d_vec_local(i)
               END DO
            END DO
            !$OMP END PARALLEL DO

            ! ---------------------------------------------------------------------
            ! F. Solve — BLAS dpotrf/dpotrs or ScaLAPACK pdpotrf/pdpotrs
            ! ---------------------------------------------------------------------
            IF (n_procs_per_atom == 1) THEN
               CALL dpotrf('L', n_local_grid, D_local, n_local_grid, info)
               CALL dpotrs('L', n_local_grid, n_loc_ri, D_local, n_local_grid, &
                           d_lp_local, n_local_grid, info)
               DEALLOCATE (D_local)
            ELSE
               CALL solve_D_lp_distributed(phi_local, d_vec_local, d_lp_local, &
                                           n_local_grid, n_ao_total, n_loc_ri, &
                                           bs_env%ri_rs%tikhonov, &
                                           para_env_sub, blacs_env_sub, &
                                           fm_struct_D, fm_struct_b, fm_D, fm_b, info)
            END IF

            ! ---------------------------------------------------------------------
            ! G. Post-scale solution by d_vec (common to both paths)
            ! ---------------------------------------------------------------------
            !$OMP PARALLEL DO DEFAULT(NONE) &
            !$OMP SHARED(n_loc_ri, n_local_grid, d_lp_local, d_vec_local) &
            !$OMP PRIVATE(j_ri, i) &
            !$OMP SCHEDULE(STATIC)
            DO j_ri = 1, n_loc_ri
               DO i = 1, n_local_grid
                  d_lp_local(i, j_ri) = d_lp_local(i, j_ri)*d_vec_local(i)
               END DO
            END DO
            !$OMP END PARALLEL DO

            ! RI-RS fit-error diagnostic (opt-in): reconstruct (mu nu|P) from Phi and the
            ! just-solved Z (= d_lp_local after post-scaling) and accumulate the max error.
            IF (bs_env%ri_rs%print_fit_error) THEN
               CALL ri_rs_e3c_accumulate(bs_env, ctx_3c, cell, phi_local, d_lp_local, &
                                         n_local_grid, n_loc_ri, atom_P, max_ao_size, &
                                         atom_j_mepos, atom_j_stride, &
                                         e3c_abs_max, e3c_idx_max, e3c_atom_max, &
                                         e3c_vmax, z_abs_max, &
                                         j_sq_sum, j_pair_max, j_idx_max, j_atom_max, &
                                         atommax_e3c, atommax_e3c_idx, &
                                         atommax_e3crel, atommax_e3crel_idx, &
                                         atommax_j, atommax_j_idx, &
                                         R_full, n_ao_total, fit_topk, do_j_full)
            END IF

            ! ---------------------------------------------------------------------
            ! H. Scatter Local Solution Back to Global DBCSR Matrix
            ! local_grid_idx is ascending (built by the ordered scan above), so a
            ! single walking pointer over it sweeps the chunks in one pass.
            ! Under ScaLAPACK (G>1) the d_lp_local solution is identical on all
            ! G subgroup ranks (gathered via cp_fm_get_submatrix); only the
            ! subgroup root writes to mat_Z_lP so each atom column is emitted
            ! exactly once. DBCSR routes blocks to their global owner on finalize.
            ! ---------------------------------------------------------------------
            IF (n_procs_per_atom == 1 .OR. para_env_sub%mepos == 0) THEN
               ALLOCATE (Z_blk(MAXVAL(r_blk_sizes), n_loc_ri))
               loc_ptr = 1

               DO i_blk = 1, num_grid_chunks
                  r_start = row_offset(i_blk) + 1
                  r_end = row_offset(i_blk) + r_blk_sizes(i_blk)
                  current_chunk_size = r_blk_sizes(i_blk)

                  Z_blk = 0.0_dp

                  DO WHILE (loc_ptr <= n_local_grid)
                     g = local_grid_idx(loc_ptr)
                     IF (g > r_end) EXIT
                     Z_blk(g - r_start + 1, 1:n_loc_ri) = d_lp_local(loc_ptr, 1:n_loc_ri)
                     loc_ptr = loc_ptr + 1
                  END DO

                  IF (MAXVAL(ABS(Z_blk(1:current_chunk_size, 1:n_loc_ri))) > bs_env%eps_filter) THEN
                     CALL dbcsr_put_block(mat_Z_lP, row=i_blk, col=atom_P, &
                                          block=Z_blk(1:current_chunk_size, 1:n_loc_ri))
                  END IF
               END DO

               DEALLOCATE (Z_blk)
            END IF

            DEALLOCATE (d_vec_local, d_lp_local)
            DEALLOCATE (local_grid_idx, phi_local)

         END DO

         IF (bs_env%ri_rs%print_fit_error .AND. do_j_full) THEN
            CALL para_env%sum(R_full)
            ! NOTE: J_full costs O(naux^2) memory and O(naux^3) Cholesky -- diagnostic only.
            CALL ri_rs_j_full_compute(bs_env, ctx_3c, para_env, particle_set, R_full, &
                                      n_ao_total, j_full_ok, j_full_info, &
                                      j_full_sq_sum, j_full_pair_max, &
                                      j_full_idx_max, j_full_atom_max)
         END IF

         DEALLOCATE (cutoff_ri_per_atom)

         CALL gw_3c_ctx_release(ctx_3c)

         IF (bs_env%ri_rs%print_fit_error) THEN
            CALL ri_rs_fit_error_finalize(bs_env, para_env, particle_set, atomic_kind_set, &
                                          e3c_abs_max, e3c_idx_max, e3c_atom_max, &
                                          e3c_vmax, z_abs_max, &
                                          j_sq_sum, j_pair_max, j_idx_max, j_atom_max, &
                                          atommax_e3c, atommax_e3c_idx, &
                                          atommax_e3crel, atommax_e3crel_idx, &
                                          atommax_j, atommax_j_idx, &
                                          j_full_ok, j_full_info, j_full_sq_sum, &
                                          j_full_pair_max, j_full_idx_max, j_full_atom_max, &
                                          fit_topk)
            DEALLOCATE (R_full)
            DEALLOCATE (atommax_e3c, atommax_e3crel, atommax_j)
            DEALLOCATE (atommax_e3c_idx, atommax_e3crel_idx, atommax_j_idx)
         END IF

         CALL dbcsr_finalize(mat_Z_lP)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T57,A,F7.1,A)') &
               'Computed Z_lP ', ' Execution time', m_walltime() - t1, ' s'
            WRITE (bs_env%unit_nr, '(A)') ' '
         END IF

         logger => cp_get_default_logger()

         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, key), cp_p_file)) THEN
            CALL dbcsr_binary_write(matrix=mat_Z_lP, filepath=TRIM(bs_env%prefix)//"Z_lP.matrix")
         END IF

      END IF

      DEALLOCATE (row_offset, ri_blk_sizes, col_dist_ri)
      CALL dbcsr_distribution_release(dist_Z)

      IF (n_procs_per_atom > 1) THEN
         CALL cp_blacs_env_release(blacs_env_sub)
         CALL para_env_sub%free()
         DEALLOCATE (para_env_sub)
      END IF

      DEALLOCATE (ri_rs_grid_points)

      CALL timestop(handle)

   END SUBROUTINE compute_coeff_Z_lP

! **************************************************************************************************
!> \brief Computes the dense localized RHS d_lp(l,P) = Σ_{μν} Φ_μ(r_l)·Φ_ν(r_l)·(μν|P) for one
!>        RI atom P, OMP-threaded over (atom_j, atom_k) AO-pair blocks: per thread, build the 3c
!>        block, then contract grid-chunked pair densities into a private d_lp partial; partials
!>        are reduced into d_lp at the end.
!>        Pair screening is handled inside build_3c_integral_block_ctx via the `screened` output,
!> \param bs_env ...
!> \param ctx shared 3c-integral context (gw_3c_ctx_create)
!> \param phi_val Φ_μ(r_l) on the local-sphere grid (n_grid_total × n_ao)
!> \param d_lp output (n_grid_total × n_loc_ri), zeroed by the caller, accumulated here
!> \param n_grid_total number of local-sphere grid rows
!> \param n_loc_ri number of RI functions of atom_P
!> \param atom_P ...
!> \param max_ao_size ...
!> \param atom_j_mepos ...
!> \param atom_j_stride ...
! **************************************************************************************************

   SUBROUTINE compute_d_lp(bs_env, ctx, phi_val, d_lp, n_grid_total, n_loc_ri, atom_P, &
                           max_ao_size, atom_j_mepos, atom_j_stride)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gw_3c_ctx_type), INTENT(IN)                   :: ctx
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi_val
      INTEGER, INTENT(IN)                                :: n_grid_total, n_loc_ri
      REAL(KIND=dp), INTENT(INOUT)                       :: d_lp(n_grid_total, n_loc_ri)
      INTEGER, INTENT(IN)                                :: atom_P, max_ao_size, atom_j_mepos, &
                                                            atom_j_stride

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_d_lp'
      INTEGER, PARAMETER                                 :: grid_chunk = 1024

      INTEGER                                            :: atom_j, atom_k, c, handle, j, jk_idx, &
                                                            jsize, jstart, k, ksize, kstart, l, &
                                                            l0, ri
      LOGICAL                                            :: screened
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: int_2d_prv, rho_chunk
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: int_3c_prv
      TYPE(gw_3c_ws_type)                                :: ws

      CALL timeset(routineN, handle)

      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(bs_env, ctx, phi_val, d_lp, n_grid_total, n_loc_ri, atom_P, max_ao_size, &
      !$OMP        atom_j_mepos, atom_j_stride) &
      !$OMP PRIVATE(atom_j, atom_k, c, j, jk_idx, jsize, jstart, k, ksize, kstart, l, l0, ri, &
      !$OMP         screened, int_2d_prv, rho_chunk, int_3c_prv, ws)

      CALL gw_3c_ws_create(ws, ctx)
      ALLOCATE (int_3c_prv(max_ao_size, max_ao_size, n_loc_ri))
      ALLOCATE (int_2d_prv(max_ao_size*max_ao_size, n_loc_ri))
      ALLOCATE (rho_chunk(grid_chunk, max_ao_size*max_ao_size))

      ! MPI-stride atom_j over the subgroup (atom_j_stride = 1 for the BLAS
      ! path, > 1 for the ScaLAPACK path). The OMP DO parallelizes the inner
      ! atom_k while the outer atom_j carries the MPI stride. COLLAPSE(2) is
      ! dropped because the outer stride is non-unit under ScaLAPACK; for
      ! typical natom the inner loop has plenty of work for DYNAMIC.
      !$OMP DO SCHEDULE(DYNAMIC) REDUCTION(+:d_lp)
      DO atom_j = atom_j_mepos + 1, SIZE(bs_env%i_ao_start_from_atom), atom_j_stride
         DO atom_k = 1, SIZE(bs_env%i_ao_start_from_atom)
            jstart = bs_env%i_ao_start_from_atom(atom_j)
            jsize = bs_env%i_ao_end_from_atom(atom_j) - jstart + 1
            kstart = bs_env%i_ao_start_from_atom(atom_k)
            ksize = bs_env%i_ao_end_from_atom(atom_k) - kstart + 1

            int_3c_prv(1:jsize, 1:ksize, 1:n_loc_ri) = 0.0_dp

            ! Compute B_{μν,P} = (μν|P); ctx-internal triangle-inequality screening on
            ! kind_radius sets `screened=.TRUE.` for negligible triples.
            CALL build_3c_integral_block_ctx(int_3c_prv(1:jsize, 1:ksize, 1:n_loc_ri), &
                                             ctx, ws, atom_j=atom_j, atom_k=atom_k, atom_i=atom_P, &
                                             screened=screened)

            IF (screened) CYCLE

            ! Flatten 3D B_{μν, P} tensor to 2D B_{(μν), P} matrix for BLAS
            DO ri = 1, n_loc_ri
               DO k = 1, ksize
                  DO j = 1, jsize
                     jk_idx = (k - 1)*jsize + j
                     int_2d_prv(jk_idx, ri) = int_3c_prv(j, k, ri)
                  END DO
               END DO
            END DO

            ! Pair density ρ(l, μν) = Φ_μ(r_l) Φ_ν(r_l) in grid chunks, contracted on the fly:
            ! d_{l,P} += ρ(l, μν) B_{(μν),P}  (dgemm runs serially inside the parallel region)
            DO l0 = 1, n_grid_total, grid_chunk
               c = MIN(grid_chunk, n_grid_total - l0 + 1)
               DO k = 1, ksize
                  DO j = 1, jsize
                     jk_idx = (k - 1)*jsize + j
                     DO l = 1, c
                        rho_chunk(l, jk_idx) = phi_val(l0 + l - 1, jstart + j - 1)* &
                                               phi_val(l0 + l - 1, kstart + k - 1)
                     END DO
                  END DO
               END DO
               CALL dgemm("N", "N", c, n_loc_ri, jsize*ksize, &
                          1.0_dp, rho_chunk, grid_chunk, &
                          int_2d_prv, max_ao_size*max_ao_size, &
                          1.0_dp, d_lp(l0, 1), n_grid_total)
            END DO
         END DO
      END DO
      !$OMP END DO

      DEALLOCATE (int_3c_prv, int_2d_prv, rho_chunk)
      CALL gw_3c_ws_release(ws)

      !$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE compute_d_lp

! **************************************************************************************************
!> \brief Accumulate the RI-RS 3-centre reconstruction error e3c for one RI atom P. Rebuilds the
!>        exact (image-summed) 3c block (mu nu|P) with build_3c_integral_block_ctx, mirroring
!>        compute_d_lp: the molecular (perd=0) path collapses to the single cell-(0,0,0) block and
!>        the periodic-Gamma path sums the same (cell_R,cell_S) images, so the "exact" reference is
!>        exactly the integral the solve targeted. It then forms the reconstruction
!>        int3c_approx(mu,nu,P) = sum_l Phi_mu(r_l) Phi_nu(r_l) Z_lP and updates, in place, the
!>        running max |approx - exact| with its global (mu,nu,P) / (atom_j,atom_k,atom_P) location,
!>        the max |exact| (rel denominator) and the max |Z|. Serial (opt-in diagnostic). MPI: the
!>        atom_j loop is strided over the subgroup exactly as compute_d_lp, so calling on every
!>        rank covers each block once; reduce afterwards with ri_rs_fit_error_finalize.
!> \param bs_env ...
!> \param ctx shared 3c-integral context (gw_3c_ctx_create)
!> \param cell ...
!> \param phi_val Phi_mu(r_l) on the local-sphere grid (n_grid x n_ao_total)
!> \param Z_lP solved coefficients on the local-sphere grid (n_grid x n_loc_ri)
!> \param n_grid number of local-sphere grid rows
!> \param n_loc_ri number of RI functions of atom_P
!> \param atom_P ...
!> \param max_ao_size ...
!> \param atom_j_mepos ...
!> \param atom_j_stride ...
!> \param e_abs_max running max abs reconstruction error (INOUT)
!> \param idx_max global (mu,nu,P) at the running max (INOUT)
!> \param atom_max (atom_j,atom_k,atom_P) at the running max (INOUT)
!> \param vmax_loc running max |exact 3c| (INOUT)
!> \param zmax_loc running max |Z_lP| (INOUT)
!> \param j_sq_sum running sum over (mu,nu) pairs of J_pair = R^T V_PP^-1 R, the Coulomb-norm
!>        squared of the per-pair fit-error density (INOUT). R(P) = approx(mu,nu,P)-exact(mu,nu,P)
!>        and V_PP = on-atom 2c bare-Coulomb (P|Q). Block-diagonal in atom_P (see body).
!> \param j_pair_max running max over (mu,nu) of J_pair (squared; sqrt taken at print) (INOUT)
!> \param j_idx_max global (mu,nu,P_start) at the running J_pair max (INOUT)
!> \param j_atom_max (atom_j,atom_k,atom_P) at the running J_pair max (INOUT)
!> \param atommax_e3c per-rank per-atom max abs 3c residual (INOUT)
!> \param atommax_e3c_idx global (mu,nu,P) at each per-atom 3c max (INOUT)
!> \param atommax_e3crel per-rank per-atom max relative 3c residual (INOUT)
!> \param atommax_e3crel_idx global (mu,nu,P) at each per-atom relative 3c max (INOUT)
!> \param atommax_j per-rank per-atom max J_pair, squared (INOUT)
!> \param atommax_j_idx global (mu,nu,P_start) at each per-atom J max (INOUT)
! **************************************************************************************************
   SUBROUTINE ri_rs_e3c_accumulate(bs_env, ctx, cell, phi_val, Z_lP, n_grid, n_loc_ri, &
                                   atom_P, max_ao_size, atom_j_mepos, atom_j_stride, &
                                   e_abs_max, idx_max, atom_max, vmax_loc, zmax_loc, &
                                   j_sq_sum, j_pair_max, j_idx_max, j_atom_max, &
                                   atommax_e3c, atommax_e3c_idx, &
                                   atommax_e3crel, atommax_e3crel_idx, atommax_j, atommax_j_idx, &
                                   R_full, n_ao_total, fit_topk, do_j_full)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gw_3c_ctx_type), INTENT(IN)                   :: ctx
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi_val
      INTEGER, INTENT(IN)                                :: n_grid, n_loc_ri, atom_P, max_ao_size, &
                                                            atom_j_mepos, atom_j_stride
      REAL(KIND=dp), INTENT(IN)                          :: Z_lP(n_grid, n_loc_ri)
      REAL(KIND=dp), INTENT(INOUT)                       :: e_abs_max, vmax_loc, zmax_loc
      INTEGER, DIMENSION(3), INTENT(INOUT)               :: idx_max, atom_max
      REAL(KIND=dp), INTENT(INOUT)                       :: j_sq_sum, j_pair_max
      INTEGER, DIMENSION(3), INTENT(INOUT)               :: j_idx_max, j_atom_max
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: atommax_e3c, atommax_e3crel, &
                                                            atommax_j
      INTEGER, DIMENSION(:, :), INTENT(INOUT)            :: atommax_e3c_idx, atommax_e3crel_idx, &
                                                            atommax_j_idx
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: R_full
      INTEGER, INTENT(IN)                                :: n_ao_total
      TYPE(ri_rs_topk_type), INTENT(INOUT)               :: fit_topk
      LOGICAL, INTENT(IN)                                :: do_j_full

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ri_rs_e3c_accumulate'
      INTEGER, PARAMETER                                 :: grid_chunk = 1024

      INTEGER :: atom_j, atom_k, c, handle, ix_hi, ix_lo, ix_R, ix_S, iy_hi, iy_lo, iy_R, iy_S, &
         iz_hi, iz_lo, iz_R, iz_S, j, jk_idx, jsize, jstart, k, ksize, kstart, l, l0, n_jk, &
         natom, P_start, ri
      INTEGER :: egfp, egfq, ikind_ri, info_solve, info_solve_trunc, info_vpp, info_vpp_trunc, &
         j_rec, jk, k_rec, maxl_ri, ncop, ncoq, nset_ri, nsgf_ri, pset, qset, sgfp, sgfq
      INTEGER                                            :: jk_global
      INTEGER, DIMENSION(3)                              :: cell_R_vec, cell_S_vec
      INTEGER, DIMENSION(:), POINTER                     :: lmax_ri, lmin_ri, npgf_ri, nsgf_set_ri
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_ri
      LOGICAL                                            :: any_kept, have_metric, j_ok, screened, trunc_ok
      REAL(KIND=dp)                                      :: err, exact, j_p, rel, resid
      REAL(KIND=dp), DIMENSION(3)                        :: pos_P_2c
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgf_ri, sphi_ri, zet_ri
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: approx_2d, R_blk, rho_chunk, V_PP, &
                                                            V_PP_chol, V_PP_trunc_chol, work2c, &
                                                            Y_blk, Y_trunc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: int_3c_prv, int_3c_sum
      TYPE(cp_libint_t)                                  :: lib2c
      TYPE(gw_3c_ws_type)                                :: ws
      TYPE(libint_potential_type)                        :: v_coulomb, v_trunc

      CALL timeset(routineN, handle)

      natom = SIZE(bs_env%i_ao_start_from_atom)
      P_start = bs_env%i_RI_start_from_atom(atom_P)

      zmax_loc = MAX(zmax_loc, MAXVAL(ABS(Z_lP(1:n_grid, 1:n_loc_ri))))

      IF (cell%perd(1) == 1) THEN; ix_lo = -1; ix_hi = 1; ELSE; ix_lo = 0; ix_hi = 0; END IF
      IF (cell%perd(2) == 1) THEN; iy_lo = -1; iy_hi = 1; ELSE; iy_lo = 0; iy_hi = 0; END IF
      IF (cell%perd(3) == 1) THEN; iz_lo = -1; iz_hi = 1; ELSE; iz_lo = 0; iz_hi = 0; END IF

      CALL gw_3c_ws_create(ws, ctx)
      ALLOCATE (int_3c_prv(max_ao_size, max_ao_size, n_loc_ri))
      ALLOCATE (int_3c_sum(max_ao_size, max_ao_size, n_loc_ri))
      ALLOCATE (approx_2d(max_ao_size*max_ao_size, n_loc_ri))
      ALLOCATE (rho_chunk(grid_chunk, max_ao_size*max_ao_size))

      ! ------------------------------------------------------------------------------------------
      ! J (Coulomb-norm) diagnostic setup: build the on-atom 2c Coulomb block V_PP = (P|Q) over
      ! atom_P's RI aux and Cholesky-factor it. Per pair (mu nu), J_pair = R^T V_PP^-1 R with
      ! R(P) = approx(mu,nu,P) - exact(mu,nu,P): the Coulomb self-energy of the density-fit error,
      ! matching what the offline rirs_common.coulomb_J computes (its E is the coefficient-space
      ! residual, so its trace E^T V E == R^T V^-1 R in integral-residual terms). BARE Coulomb is
      ! used (offline uses int2c2e); for an on-atom block the truncated production metric differs
      ! only by tails beyond the cutoff radius, negligible vs the per-atom aux extent. The fit is
      ! block-diagonal in atom_P: cross-atom-P Coulomb coupling is dropped (the full coupled J is
      ! incompatible with the memory-lean per-atom-P flow). Exact for single-atom references.
      ! ------------------------------------------------------------------------------------------
      ikind_ri = ctx%kind_of(atom_P)
      CALL get_gto_basis_set(ctx%basis_i(ikind_ri)%gto_basis_set, lmax=lmax_ri, lmin=lmin_ri, &
                             npgf=npgf_ri, nset=nset_ri, zet=zet_ri, first_sgf=first_sgf_ri, &
                             nsgf_set=nsgf_set_ri, pgf_radius=rpgf_ri, sphi=sphi_ri, &
                             maxl=maxl_ri, nsgf=nsgf_ri)

      have_metric = (nsgf_ri == n_loc_ri)
      j_ok = have_metric
      trunc_ok = have_metric
      IF (have_metric) THEN
         pos_P_2c(1:3) = ctx%particle_set(atom_P)%r(1:3)
         CALL init_md_ftable(2*maxl_ri)
         v_coulomb%potential_type = do_potential_coulomb
         CALL cp_libint_init_2eri(lib2c, maxl_ri)
         CALL cp_libint_set_contrdepth(lib2c, 1)

         ALLOCATE (V_PP(n_loc_ri, n_loc_ri))
         V_PP(:, :) = 0.0_dp
         DO pset = 1, nset_ri
            ncop = npgf_ri(pset)*ncoset(lmax_ri(pset))
            sgfp = first_sgf_ri(1, pset)
            egfp = sgfp + nsgf_set_ri(pset) - 1
            DO qset = 1, nset_ri
               ncoq = npgf_ri(qset)*ncoset(lmax_ri(qset))
               sgfq = first_sgf_ri(1, qset)
               egfq = sgfq + nsgf_set_ri(qset) - 1
               ALLOCATE (work2c(ncop, ncoq))
               work2c(:, :) = 0.0_dp
               CALL eri_2center(work2c, lmin_ri(pset), lmax_ri(pset), npgf_ri(pset), &
                                zet_ri(:, pset), rpgf_ri(:, pset), pos_P_2c, &
                                lmin_ri(qset), lmax_ri(qset), npgf_ri(qset), &
                                zet_ri(:, qset), rpgf_ri(:, qset), pos_P_2c, &
                                0.0_dp, lib2c, v_coulomb)
               CALL ab_contract(V_PP(sgfp:egfp, sgfq:egfq), work2c, &
                                sphi_ri(:, sgfp:), sphi_ri(:, sgfq:), &
                                ncop, ncoq, nsgf_set_ri(pset), nsgf_set_ri(qset))
               DEALLOCATE (work2c)
            END DO
         END DO
         CALL cp_libint_cleanup_2eri(lib2c)

         ALLOCATE (V_PP_chol(n_loc_ri, n_loc_ri))
         V_PP_chol(:, :) = V_PP(:, :)
         CALL dpotrf('L', n_loc_ri, V_PP_chol, n_loc_ri, info_vpp)
         IF (info_vpp /= 0) j_ok = .FALSE.
         DEALLOCATE (V_PP)
         IF (j_ok) ALLOCATE (Y_blk(n_loc_ri, max_ao_size*max_ao_size))

         v_trunc = bs_env%ri_metric
         CALL cp_libint_init_2eri(lib2c, maxl_ri)
         CALL cp_libint_set_contrdepth(lib2c, 1)

         ALLOCATE (V_PP(n_loc_ri, n_loc_ri))
         V_PP(:, :) = 0.0_dp
         DO pset = 1, nset_ri
            ncop = npgf_ri(pset)*ncoset(lmax_ri(pset))
            sgfp = first_sgf_ri(1, pset)
            egfp = sgfp + nsgf_set_ri(pset) - 1
            DO qset = 1, nset_ri
               ncoq = npgf_ri(qset)*ncoset(lmax_ri(qset))
               sgfq = first_sgf_ri(1, qset)
               egfq = sgfq + nsgf_set_ri(qset) - 1
               ALLOCATE (work2c(ncop, ncoq))
               work2c(:, :) = 0.0_dp
               CALL eri_2center(work2c, lmin_ri(pset), lmax_ri(pset), npgf_ri(pset), &
                                zet_ri(:, pset), rpgf_ri(:, pset), pos_P_2c, &
                                lmin_ri(qset), lmax_ri(qset), npgf_ri(qset), &
                                zet_ri(:, qset), rpgf_ri(:, qset), pos_P_2c, &
                                0.0_dp, lib2c, v_trunc)
               CALL ab_contract(V_PP(sgfp:egfp, sgfq:egfq), work2c, &
                                sphi_ri(:, sgfp:), sphi_ri(:, sgfq:), &
                                ncop, ncoq, nsgf_set_ri(pset), nsgf_set_ri(qset))
               DEALLOCATE (work2c)
            END DO
         END DO
         CALL cp_libint_cleanup_2eri(lib2c)

         ALLOCATE (V_PP_trunc_chol(n_loc_ri, n_loc_ri))
         V_PP_trunc_chol(:, :) = V_PP(:, :)
         CALL dpotrf('L', n_loc_ri, V_PP_trunc_chol, n_loc_ri, info_vpp_trunc)
         IF (info_vpp_trunc /= 0) trunc_ok = .FALSE.
         DEALLOCATE (V_PP)
         IF (trunc_ok) ALLOCATE (Y_trunc(n_loc_ri, max_ao_size*max_ao_size))
         IF (j_ok .OR. trunc_ok) ALLOCATE (R_blk(n_loc_ri, max_ao_size*max_ao_size))
      END IF

      DO atom_j = atom_j_mepos + 1, natom, atom_j_stride
         DO atom_k = 1, natom
            jstart = bs_env%i_ao_start_from_atom(atom_j)
            jsize = bs_env%i_ao_end_from_atom(atom_j) - jstart + 1
            kstart = bs_env%i_ao_start_from_atom(atom_k)
            ksize = bs_env%i_ao_end_from_atom(atom_k) - kstart + 1
            n_jk = jsize*ksize

            ! Exact (image-summed) 3c block; screened far-pair triples contribute ~0 both to the
            ! exact integral and to the reconstruction (their AO tails vanish on atom_P's sphere),
            ! so skipping them matches compute_d_lp and the offline single-atom reference.
            int_3c_sum(1:jsize, 1:ksize, 1:n_loc_ri) = 0.0_dp
            any_kept = .FALSE.

            DO ix_R = ix_lo, ix_hi
            DO iy_R = iy_lo, iy_hi
            DO iz_R = iz_lo, iz_hi
               cell_R_vec = [ix_R, iy_R, iz_R]
               DO ix_S = ix_lo, ix_hi
               DO iy_S = iy_lo, iy_hi
               DO iz_S = iz_lo, iz_hi
                  cell_S_vec = [ix_S, iy_S, iz_S]
                  int_3c_prv(1:jsize, 1:ksize, 1:n_loc_ri) = 0.0_dp
                  CALL build_3c_integral_block_ctx( &
                     int_3c_prv(1:jsize, 1:ksize, 1:n_loc_ri), ctx, ws, &
                     atom_j=atom_j, atom_k=atom_k, atom_i=atom_P, &
                     cell_j=cell_R_vec, cell_k=cell_S_vec, cell_i=[0, 0, 0], &
                     screened=screened)
                  IF (screened) CYCLE
                  any_kept = .TRUE.
                  int_3c_sum(1:jsize, 1:ksize, 1:n_loc_ri) = &
                     int_3c_sum(1:jsize, 1:ksize, 1:n_loc_ri) + &
                     int_3c_prv(1:jsize, 1:ksize, 1:n_loc_ri)
               END DO
               END DO
               END DO
            END DO
            END DO
            END DO

            IF (.NOT. any_kept) CYCLE

            ! int3c_approx((mu nu), P) = sum_l Phi_mu(r_l) Phi_nu(r_l) Z_lP(l, P)  (rho^T Z)
            approx_2d(1:n_jk, 1:n_loc_ri) = 0.0_dp
            DO l0 = 1, n_grid, grid_chunk
               c = MIN(grid_chunk, n_grid - l0 + 1)
               DO k = 1, ksize
                  DO j = 1, jsize
                     jk_idx = (k - 1)*jsize + j
                     DO l = 1, c
                        rho_chunk(l, jk_idx) = phi_val(l0 + l - 1, jstart + j - 1)* &
                                               phi_val(l0 + l - 1, kstart + k - 1)
                     END DO
                  END DO
               END DO
               CALL dgemm("T", "N", n_jk, n_loc_ri, c, &
                          1.0_dp, rho_chunk, grid_chunk, &
                          Z_lP(l0, 1), n_grid, &
                          1.0_dp, approx_2d, max_ao_size*max_ao_size)
            END DO

            DO ri = 1, n_loc_ri
               DO k = 1, ksize
                  DO j = 1, jsize
                     jk_idx = (k - 1)*jsize + j
                     exact = int_3c_sum(j, k, ri)
                     vmax_loc = MAX(vmax_loc, ABS(exact))
                     resid = approx_2d(jk_idx, ri) - exact
                     err = ABS(resid)
                     rel = err/MAX(ABS(exact), ri_rs_relfloor)
                     IF (err > e_abs_max) THEN
                        e_abs_max = err
                        idx_max(1) = jstart + j - 1
                        idx_max(2) = kstart + k - 1
                        idx_max(3) = P_start + ri - 1
                        atom_max(1) = atom_j
                        atom_max(2) = atom_k
                        atom_max(3) = atom_P
                     END IF
                     CALL ri_rs_atommax_update(atommax_e3c, atommax_e3c_idx, atom_j, err, &
                                               jstart + j - 1, kstart + k - 1, P_start + ri - 1)
                     CALL ri_rs_atommax_update(atommax_e3c, atommax_e3c_idx, atom_k, err, &
                                               jstart + j - 1, kstart + k - 1, P_start + ri - 1)
                     CALL ri_rs_atommax_update(atommax_e3c, atommax_e3c_idx, atom_P, err, &
                                               jstart + j - 1, kstart + k - 1, P_start + ri - 1)
                     CALL ri_rs_atommax_update(atommax_e3crel, atommax_e3crel_idx, atom_j, rel, &
                                               jstart + j - 1, kstart + k - 1, P_start + ri - 1)
                     CALL ri_rs_atommax_update(atommax_e3crel, atommax_e3crel_idx, atom_k, rel, &
                                               jstart + j - 1, kstart + k - 1, P_start + ri - 1)
                     CALL ri_rs_atommax_update(atommax_e3crel, atommax_e3crel_idx, atom_P, rel, &
                                               jstart + j - 1, kstart + k - 1, P_start + ri - 1)
                     CALL ri_rs_topk_insert(fit_topk, ri_rs_topk_e3c, resid, ABS(exact), &
                                            jstart + j - 1, kstart + k - 1, &
                                            P_start + ri - 1, atom_j, atom_k, atom_P)
                     CALL ri_rs_topk_insert(fit_topk, ri_rs_topk_e3crel, rel, ABS(exact), &
                                            jstart + j - 1, kstart + k - 1, &
                                            P_start + ri - 1, atom_j, atom_k, atom_P)
                     IF (do_j_full) THEN
                        jk_global = (kstart + k - 2)*n_ao_total + jstart + j - 1
                        R_full(jk_global, P_start + ri - 1) = resid
                     END IF
                     IF (j_ok .OR. trunc_ok) R_blk(ri, jk_idx) = resid
                  END DO
               END DO
            END DO

            ! J: per pair (mu nu) solve V_PP y = R and accumulate J_pair = R^T y = R^T V_PP^-1 R.
            IF (j_ok) THEN
               Y_blk(1:n_loc_ri, 1:n_jk) = R_blk(1:n_loc_ri, 1:n_jk)
               CALL dpotrs('L', n_loc_ri, n_jk, V_PP_chol, n_loc_ri, Y_blk, n_loc_ri, info_solve)
               DO jk = 1, n_jk
                  k_rec = (jk - 1)/jsize + 1
                  j_rec = jk - (k_rec - 1)*jsize
                  j_p = 0.0_dp
                  DO ri = 1, n_loc_ri
                     j_p = j_p + R_blk(ri, jk)*Y_blk(ri, jk)
                     IF (info_solve == 0) THEN
                        CALL ri_rs_topk_insert(fit_topk, ri_rs_topk_mbare, Y_blk(ri, jk), 0.0_dp, &
                                               jstart + j_rec - 1, kstart + k_rec - 1, &
                                               P_start + ri - 1, atom_j, atom_k, atom_P)
                     END IF
                  END DO
                  j_sq_sum = j_sq_sum + j_p
                  IF (j_p > j_pair_max) THEN
                     j_pair_max = j_p
                     j_idx_max(1) = jstart + j_rec - 1
                     j_idx_max(2) = kstart + k_rec - 1
                     j_idx_max(3) = P_start
                     j_atom_max(1) = atom_j
                     j_atom_max(2) = atom_k
                     j_atom_max(3) = atom_P
                  END IF
                  CALL ri_rs_atommax_update(atommax_j, atommax_j_idx, atom_j, j_p, &
                                            jstart + j_rec - 1, kstart + k_rec - 1, P_start)
                  CALL ri_rs_atommax_update(atommax_j, atommax_j_idx, atom_k, j_p, &
                                            jstart + j_rec - 1, kstart + k_rec - 1, P_start)
               END DO
            END IF
            IF (trunc_ok) THEN
               Y_trunc(1:n_loc_ri, 1:n_jk) = R_blk(1:n_loc_ri, 1:n_jk)
               CALL dpotrs('L', n_loc_ri, n_jk, V_PP_trunc_chol, n_loc_ri, &
                           Y_trunc, n_loc_ri, info_solve_trunc)
               IF (info_solve_trunc == 0) THEN
                  DO jk = 1, n_jk
                     k_rec = (jk - 1)/jsize + 1
                     j_rec = jk - (k_rec - 1)*jsize
                     DO ri = 1, n_loc_ri
                        CALL ri_rs_topk_insert(fit_topk, ri_rs_topk_mtrunc, Y_trunc(ri, jk), 0.0_dp, &
                                               jstart + j_rec - 1, kstart + k_rec - 1, &
                                               P_start + ri - 1, atom_j, atom_k, atom_P)
                     END DO
                  END DO
               END IF
            END IF
         END DO
      END DO

      DEALLOCATE (int_3c_prv, int_3c_sum, approx_2d, rho_chunk)
      IF (ALLOCATED(V_PP_chol)) DEALLOCATE (V_PP_chol)
      IF (ALLOCATED(V_PP_trunc_chol)) DEALLOCATE (V_PP_trunc_chol)
      IF (ALLOCATED(R_blk)) DEALLOCATE (R_blk)
      IF (ALLOCATED(Y_blk)) DEALLOCATE (Y_blk)
      IF (ALLOCATED(Y_trunc)) DEALLOCATE (Y_trunc)
      CALL gw_3c_ws_release(ws)

      CALL timestop(handle)

   END SUBROUTINE ri_rs_e3c_accumulate

! **************************************************************************************************
!> \brief Build the full bare-Coulomb auxiliary metric and evaluate R_full^T V_full^-1 R_full.
!>        Diagnostic only: replicated on rank 0 after R_full has been summed over para_env.
!> \param bs_env ...
!> \param ctx shared 3c-integral context (gw_3c_ctx_create)
!> \param para_env full RI-RS communicator
!> \param particle_set ...
!> \param R_full global residual rows (mu,nu), columns auxiliary P
!> \param n_ao_total total number of AO functions
!> \param j_full_ok true if the full Cholesky/solve succeeded on rank 0
!> \param j_full_info dpotrf/dpotrs info, or negative for a basis-size mismatch
!> \param j_sq_sum rank-local contribution to the full J sum (rank 0 only)
!> \param j_pair_max rank-local max over AO pairs (rank 0 only)
!> \param j_idx_max global (mu,nu,0) at the full-J pair max
!> \param j_atom_max atom indices (atom_mu,atom_nu,0) at the full-J pair max
! **************************************************************************************************
   SUBROUTINE ri_rs_j_full_compute(bs_env, ctx, para_env, particle_set, R_full, n_ao_total, &
                                   j_full_ok, j_full_info, j_sq_sum, j_pair_max, &
                                   j_idx_max, j_atom_max)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gw_3c_ctx_type), INTENT(IN)                   :: ctx
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: R_full
      INTEGER, INTENT(IN)                                :: n_ao_total
      LOGICAL, INTENT(OUT)                               :: j_full_ok
      INTEGER, INTENT(OUT)                               :: j_full_info
      REAL(KIND=dp), INTENT(OUT)                         :: j_sq_sum, j_pair_max
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: j_idx_max, j_atom_max

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ri_rs_j_full_compute'

      INTEGER :: atom_mu, atom_P, atom_Q, col_aux, egfp, egfq, handle, ikind_P, ikind_Q, &
         info_solve, maxl_P, maxl_Q, maxl_pair, mu_idx, natom, naux_total, ncoq, ncop, &
         n_ao_pair_total, nu_idx, p_nri, p_start, pset, q_nri, q_start, qset, row_pair, sgfp, sgfq
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_atom
      INTEGER, DIMENSION(:), POINTER                     :: lmax_P, lmax_Q, lmin_P, lmin_Q, &
                                                            npgf_P, npgf_Q, nsgf_set_P, nsgf_set_Q
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_P, first_sgf_Q
      REAL(KIND=dp)                                      :: j_p
      REAL(KIND=dp), DIMENSION(3)                        :: pos_P, pos_Q
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgf_P, rpgf_Q, sphi_P, sphi_Q, &
                                                            zet_P, zet_Q
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: V_full, work2c, Y_full
      TYPE(cp_libint_t)                                  :: lib2c
      TYPE(libint_potential_type)                        :: v_coulomb

      CALL timeset(routineN, handle)

      j_full_ok = .FALSE.
      j_full_info = 0
      j_sq_sum = 0.0_dp
      j_pair_max = -1.0_dp
      j_idx_max(:) = 0
      j_atom_max(:) = 0

      IF (para_env%mepos /= 0) THEN
         CALL timestop(handle)
         RETURN
      END IF

      natom = SIZE(bs_env%i_RI_start_from_atom)
      naux_total = SIZE(R_full, 2)
      n_ao_pair_total = SIZE(R_full, 1)

      ALLOCATE (ao_atom(n_ao_total))
      DO atom_mu = 1, natom
         ao_atom(bs_env%i_ao_start_from_atom(atom_mu):bs_env%i_ao_end_from_atom(atom_mu)) = atom_mu
      END DO

      ALLOCATE (V_full(naux_total, naux_total))
      V_full(:, :) = 0.0_dp
      v_coulomb%potential_type = do_potential_coulomb

      DO atom_P = 1, natom
         p_start = bs_env%i_RI_start_from_atom(atom_P)
         p_nri = bs_env%i_RI_end_from_atom(atom_P) - p_start + 1
         ikind_P = ctx%kind_of(atom_P)
         CALL get_gto_basis_set(ctx%basis_i(ikind_P)%gto_basis_set, lmax=lmax_P, lmin=lmin_P, &
                                npgf=npgf_P, nset=pset, zet=zet_P, first_sgf=first_sgf_P, &
                                nsgf_set=nsgf_set_P, pgf_radius=rpgf_P, sphi=sphi_P, &
                                maxl=maxl_P, nsgf=col_aux)
         IF (col_aux /= p_nri) THEN
            j_full_info = -1
            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,I0,A,I0)') &
                  'RI-RS fit error | WARNING: J_full skipped; RI basis size mismatch on atom ', &
                  atom_P, ' nsgf=', col_aux
            END IF
            DEALLOCATE (V_full, ao_atom)
            CALL timestop(handle)
            RETURN
         END IF

         pos_P(1:3) = particle_set(atom_P)%r(1:3)

         DO atom_Q = 1, natom
            q_start = bs_env%i_RI_start_from_atom(atom_Q)
            q_nri = bs_env%i_RI_end_from_atom(atom_Q) - q_start + 1
            ikind_Q = ctx%kind_of(atom_Q)
            CALL get_gto_basis_set(ctx%basis_i(ikind_Q)%gto_basis_set, lmax=lmax_Q, lmin=lmin_Q, &
                                   npgf=npgf_Q, nset=qset, zet=zet_Q, first_sgf=first_sgf_Q, &
                                   nsgf_set=nsgf_set_Q, pgf_radius=rpgf_Q, sphi=sphi_Q, &
                                   maxl=maxl_Q, nsgf=col_aux)
            IF (col_aux /= q_nri) THEN
               j_full_info = -1
               IF (bs_env%unit_nr > 0) THEN
                  WRITE (bs_env%unit_nr, '(T2,A,I0,A,I0)') &
                     'RI-RS fit error | WARNING: J_full skipped; RI basis size mismatch on atom ', &
                     atom_Q, ' nsgf=', col_aux
               END IF
               DEALLOCATE (V_full, ao_atom)
               CALL timestop(handle)
               RETURN
            END IF

            pos_Q(1:3) = particle_set(atom_Q)%r(1:3)
            maxl_pair = MAX(maxl_P, maxl_Q)
            CALL init_md_ftable(2*maxl_pair)
            CALL cp_libint_init_2eri(lib2c, maxl_pair)
            CALL cp_libint_set_contrdepth(lib2c, 1)

            DO pset = 1, SIZE(lmax_P)
               ncop = npgf_P(pset)*ncoset(lmax_P(pset))
               sgfp = first_sgf_P(1, pset)
               egfp = sgfp + nsgf_set_P(pset) - 1
               DO qset = 1, SIZE(lmax_Q)
                  ncoq = npgf_Q(qset)*ncoset(lmax_Q(qset))
                  sgfq = first_sgf_Q(1, qset)
                  egfq = sgfq + nsgf_set_Q(qset) - 1
                  ALLOCATE (work2c(ncop, ncoq))
                  work2c(:, :) = 0.0_dp
                  CALL eri_2center(work2c, lmin_P(pset), lmax_P(pset), npgf_P(pset), &
                                   zet_P(:, pset), rpgf_P(:, pset), pos_P, &
                                   lmin_Q(qset), lmax_Q(qset), npgf_Q(qset), &
                                   zet_Q(:, qset), rpgf_Q(:, qset), pos_Q, &
                                   0.0_dp, lib2c, v_coulomb)
                  CALL ab_contract(V_full(p_start + sgfp - 1:p_start + egfp - 1, &
                                          q_start + sgfq - 1:q_start + egfq - 1), work2c, &
                                   sphi_P(:, sgfp:), sphi_Q(:, sgfq:), &
                                   ncop, ncoq, nsgf_set_P(pset), nsgf_set_Q(qset))
                  DEALLOCATE (work2c)
               END DO
            END DO

            CALL cp_libint_cleanup_2eri(lib2c)
         END DO
      END DO

      CALL dpotrf('L', naux_total, V_full, naux_total, j_full_info)
      IF (j_full_info /= 0) THEN
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I0)') &
               'RI-RS fit error | WARNING: J_full skipped; dpotrf info = ', j_full_info
         END IF
         DEALLOCATE (V_full, ao_atom)
         CALL timestop(handle)
         RETURN
      END IF

      ALLOCATE (Y_full(naux_total, n_ao_pair_total))
      DO row_pair = 1, n_ao_pair_total
         Y_full(1:naux_total, row_pair) = R_full(row_pair, 1:naux_total)
      END DO

      CALL dpotrs('L', naux_total, n_ao_pair_total, V_full, naux_total, &
                  Y_full, naux_total, info_solve)
      IF (info_solve /= 0) THEN
         j_full_info = info_solve
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I0)') &
               'RI-RS fit error | WARNING: J_full skipped; dpotrs info = ', info_solve
         END IF
         DEALLOCATE (V_full, Y_full, ao_atom)
         CALL timestop(handle)
         RETURN
      END IF

      DO row_pair = 1, n_ao_pair_total
         j_p = 0.0_dp
         DO col_aux = 1, naux_total
            j_p = j_p + R_full(row_pair, col_aux)*Y_full(col_aux, row_pair)
         END DO
         j_sq_sum = j_sq_sum + j_p
         IF (j_p > j_pair_max) THEN
            j_pair_max = j_p
            mu_idx = MOD(row_pair - 1, n_ao_total) + 1
            nu_idx = (row_pair - 1)/n_ao_total + 1
            j_idx_max(1) = mu_idx
            j_idx_max(2) = nu_idx
            j_idx_max(3) = 0
            j_atom_max(1) = ao_atom(mu_idx)
            j_atom_max(2) = ao_atom(nu_idx)
            j_atom_max(3) = 0
         END IF
      END DO

      j_full_ok = .TRUE.

      DEALLOCATE (V_full, Y_full, ao_atom)

      CALL timestop(handle)

   END SUBROUTINE ri_rs_j_full_compute

! **************************************************************************************************
!> \brief Global MAXLOC reduction + print of the RI-RS fit-error diagnostics accumulated per atom
!>        by ri_rs_e3c_accumulate. Collective on para_env (every rank calls the reductions/bcast;
!>        only bs_env%unit_nr > 0 writes). Shared by the molecular and periodic-Gamma RI-RS paths.
!> \param bs_env ...
!> \param para_env full RI-RS communicator
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param e_abs_max per-rank max abs reconstruction error
!> \param idx_max per-rank global (mu,nu,P) at the max
!> \param atom_max per-rank (atom_j,atom_k,atom_P) at the max
!> \param vmax_loc per-rank max |exact 3c|
!> \param zmax_loc per-rank max |Z_lP|
!> \param j_sq_sum per-rank sum over pairs of J_pair = R^T V_PP^-1 R (Coulomb-norm squared)
!> \param j_pair_max per-rank max over pairs of J_pair (squared)
!> \param j_idx_max per-rank global (mu,nu,P_start) at the J_pair max
!> \param j_atom_max per-rank (atom_j,atom_k,atom_P) at the J_pair max
!> \param atommax_e3c per-rank per-atom max abs 3c residual
!> \param atommax_e3c_idx global (mu,nu,P) at each per-atom 3c max
!> \param atommax_e3crel per-rank per-atom max relative 3c residual
!> \param atommax_e3crel_idx global (mu,nu,P) at each per-atom relative 3c max
!> \param atommax_j per-rank per-atom max J_pair, squared
!> \param atommax_j_idx global (mu,nu,P_start) at each per-atom J max
!> \param j_full_ok true on the rank that successfully computed the full-J diagnostic
!> \param j_full_info dpotrf/dpotrs status, or negative for setup failure
!> \param j_full_sq_sum rank-local full-J sum (rank 0 only)
!> \param j_full_pair_max rank-local full-J max over AO pairs (rank 0 only)
!> \param j_full_idx_max global (mu,nu,0) at the full-J pair max
!> \param j_full_atom_max atom indices (atom_mu,atom_nu,0) at the full-J pair max
! **************************************************************************************************
   SUBROUTINE ri_rs_fit_error_finalize(bs_env, para_env, particle_set, atomic_kind_set, &
                                       e_abs_max, idx_max, atom_max, vmax_loc, zmax_loc, &
                                       j_sq_sum, j_pair_max, j_idx_max, j_atom_max, &
                                       atommax_e3c, atommax_e3c_idx, &
                                       atommax_e3crel, atommax_e3crel_idx, atommax_j, atommax_j_idx, &
                                       j_full_ok, j_full_info, j_full_sq_sum, j_full_pair_max, &
                                       j_full_idx_max, j_full_atom_max, fit_topk)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      REAL(KIND=dp), INTENT(IN)                          :: e_abs_max, vmax_loc, zmax_loc
      INTEGER, DIMENSION(3), INTENT(IN)                  :: idx_max, atom_max
      REAL(KIND=dp), INTENT(IN)                          :: j_sq_sum, j_pair_max
      INTEGER, DIMENSION(3), INTENT(IN)                  :: j_idx_max, j_atom_max
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: atommax_e3c, atommax_e3crel, &
                                                            atommax_j
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: atommax_e3c_idx, atommax_e3crel_idx, &
                                                            atommax_j_idx
      LOGICAL, INTENT(IN)                                :: j_full_ok
      INTEGER, INTENT(IN)                                :: j_full_info
      REAL(KIND=dp), INTENT(IN)                          :: j_full_sq_sum, j_full_pair_max
      INTEGER, DIMENSION(3), INTENT(IN)                  :: j_full_idx_max, j_full_atom_max
      TYPE(ri_rs_topk_type), INTENT(IN)                  :: fit_topk

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ri_rs_fit_error_finalize'

      CHARACTER(LEN=2)                                   :: sym_fmu, sym_fk, sym_jk, sym_jmu, sym_jP, &
                                                            sym_k, sym_mu, sym_P
      INTEGER                                            :: handle, ik, j_full_info_glob, j_full_ok_i, src, src_j, &
                                                            src_j_full
      INTEGER, DIMENSION(4)                              :: loc_buf_j_full
      INTEGER, DIMENSION(6)                              :: loc_buf, loc_buf_j
      LOGICAL                                            :: j_full_ok_glob
      REAL(KIND=dp)                                      :: e_abs_glob, e_rel_glob, j_pair_glob, &
                                                            j_full_pair_glob, j_full_sq_glob, &
                                                            sqrt_j_full_glob, sqrt_j_full_pair_glob, &
                                                            j_sq_glob, sqrt_j_glob, &
                                                            sqrt_jpair_glob, v_glob, z_glob

      CALL timeset(routineN, handle)

      e_abs_glob = e_abs_max
      CALL para_env%max(e_abs_glob)
      v_glob = vmax_loc
      CALL para_env%max(v_glob)
      z_glob = zmax_loc
      CALL para_env%max(z_glob)

      ! Locate the rank holding the global max abs error (ties -> highest rank), broadcast its
      ! (mu,nu,P) and (atom_j,atom_k,atom_P). MAX of doubles is exact, so == is safe here.
      src = -1
      IF (e_abs_max == e_abs_glob) src = para_env%mepos
      CALL para_env%max(src)

      loc_buf(1:6) = 0
      IF (para_env%mepos == src) THEN
         loc_buf(1:3) = idx_max(1:3)
         loc_buf(4:6) = atom_max(1:3)
      END IF
      CALL para_env%bcast(loc_buf, src)

      ! J (Coulomb-norm): total = sqrt(sum_pairs J_pair); per-pair worst = sqrt(max_pairs J_pair).
      j_sq_glob = j_sq_sum
      CALL para_env%sum(j_sq_glob)
      sqrt_j_glob = SQRT(MAX(j_sq_glob, 0.0_dp))
      j_pair_glob = j_pair_max
      CALL para_env%max(j_pair_glob)
      sqrt_jpair_glob = SQRT(MAX(j_pair_glob, 0.0_dp))

      src_j = -1
      IF (j_pair_max == j_pair_glob) src_j = para_env%mepos
      CALL para_env%max(src_j)
      loc_buf_j(1:6) = 0
      IF (para_env%mepos == src_j) THEN
         loc_buf_j(1:3) = j_idx_max(1:3)
         loc_buf_j(4:6) = j_atom_max(1:3)
      END IF
      CALL para_env%bcast(loc_buf_j, src_j)

      j_full_info_glob = ABS(j_full_info)
      CALL para_env%max(j_full_info_glob)
      j_full_sq_glob = j_full_sq_sum
      CALL para_env%sum(j_full_sq_glob)
      sqrt_j_full_glob = SQRT(MAX(j_full_sq_glob, 0.0_dp))
      j_full_pair_glob = j_full_pair_max
      CALL para_env%max(j_full_pair_glob)
      sqrt_j_full_pair_glob = SQRT(MAX(j_full_pair_glob, 0.0_dp))
      j_full_ok_i = 0
      IF (j_full_ok) j_full_ok_i = 1
      CALL para_env%max(j_full_ok_i)
      j_full_ok_glob = (j_full_info_glob == 0 .AND. j_full_pair_glob >= 0.0_dp .AND. &
                        j_full_ok_i == 1)

      loc_buf_j_full(1:4) = 0
      IF (j_full_ok_glob) THEN
         src_j_full = -1
         IF (j_full_pair_max == j_full_pair_glob) src_j_full = para_env%mepos
         CALL para_env%max(src_j_full)
         IF (para_env%mepos == src_j_full) THEN
            loc_buf_j_full(1:2) = j_full_idx_max(1:2)
            loc_buf_j_full(3:4) = j_full_atom_max(1:2)
         END IF
         CALL para_env%bcast(loc_buf_j_full, src_j_full)
      END IF

      e_rel_glob = e_abs_glob/MAX(v_glob, 1.0E-300_dp)

      sym_mu = "??"; sym_k = "??"; sym_P = "??"
      IF (loc_buf(4) >= 1) THEN
         ik = particle_set(loc_buf(4))%atomic_kind%kind_number
         sym_mu = atomic_kind_set(ik)%element_symbol
      END IF
      IF (loc_buf(5) >= 1) THEN
         ik = particle_set(loc_buf(5))%atomic_kind%kind_number
         sym_k = atomic_kind_set(ik)%element_symbol
      END IF
      IF (loc_buf(6) >= 1) THEN
         ik = particle_set(loc_buf(6))%atomic_kind%kind_number
         sym_P = atomic_kind_set(ik)%element_symbol
      END IF

      sym_jmu = "??"; sym_jk = "??"; sym_jP = "??"
      IF (loc_buf_j(4) >= 1) THEN
         ik = particle_set(loc_buf_j(4))%atomic_kind%kind_number
         sym_jmu = atomic_kind_set(ik)%element_symbol
      END IF
      IF (loc_buf_j(5) >= 1) THEN
         ik = particle_set(loc_buf_j(5))%atomic_kind%kind_number
         sym_jk = atomic_kind_set(ik)%element_symbol
      END IF
      IF (loc_buf_j(6) >= 1) THEN
         ik = particle_set(loc_buf_j(6))%atomic_kind%kind_number
         sym_jP = atomic_kind_set(ik)%element_symbol
      END IF

      sym_fmu = "??"; sym_fk = "??"
      IF (loc_buf_j_full(3) >= 1) THEN
         ik = particle_set(loc_buf_j_full(3))%atomic_kind%kind_number
         sym_fmu = atomic_kind_set(ik)%element_symbol
      END IF
      IF (loc_buf_j_full(4) >= 1) THEN
         ik = particle_set(loc_buf_j_full(4))%atomic_kind%kind_number
         sym_fk = atomic_kind_set(ik)%element_symbol
      END IF

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(A)') ' '
         WRITE (bs_env%unit_nr, '(T2,A,ES11.3,A,I0,A,I0,A,I0,A)') &
            'RI-RS fit error | e3c Linf abs = ', e_abs_glob, &
            '  at (mu=', loc_buf(1), ' nu=', loc_buf(2), ' P=', loc_buf(3), ')'
         WRITE (bs_env%unit_nr, '(T2,A,ES11.3,A,ES11.3)') &
            'RI-RS fit error | e3c Linf rel = ', e_rel_glob, ' ; max|Z| = ', z_glob
         WRITE (bs_env%unit_nr, '(T2,A,3(2X,I0,1X,A2))') &
            'RI-RS fit error | argmax atoms (mu nu P):', &
            loc_buf(4), sym_mu, loc_buf(5), sym_k, loc_buf(6), sym_P
         WRITE (bs_env%unit_nr, '(T2,A,ES11.3,A,ES11.3,A,I0,A,I0,A)') &
            'RI-RS fit error | J sqrt(tr E^T V E) = ', sqrt_j_glob, &
            ' ; per-pair max = ', sqrt_jpair_glob, &
            ' at (mu=', loc_buf_j(1), ' nu=', loc_buf_j(2), ')'
         WRITE (bs_env%unit_nr, '(T2,A,3(2X,I0,1X,A2))') &
            'RI-RS fit error | J argmax atoms (mu nu P):', &
            loc_buf_j(4), sym_jmu, loc_buf_j(5), sym_jk, loc_buf_j(6), sym_jP
         IF (j_full_ok_glob) THEN
            WRITE (bs_env%unit_nr, '(T2,A,ES11.3,A,ES11.3,A,I0,A,I0,A)') &
               'RI-RS fit error | J_full sqrt(tr R^T V^-1 R) = ', sqrt_j_full_glob, &
               ' ; per-pair max = ', sqrt_j_full_pair_glob, &
               ' at (mu=', loc_buf_j_full(1), ' nu=', loc_buf_j_full(2), ')'
            WRITE (bs_env%unit_nr, '(T2,A,2(2X,I0,1X,A2))') &
               'RI-RS fit error | J_full argmax atoms (mu nu):', &
               loc_buf_j_full(3), sym_fmu, loc_buf_j_full(4), sym_fk
         ELSE
            WRITE (bs_env%unit_nr, '(T2,A,I0)') &
               'RI-RS fit error | J_full skipped; info = ', j_full_info_glob
         END IF
      END IF

      CALL ri_rs_atommax_reduce_print(bs_env, para_env, particle_set, atomic_kind_set, &
                                      atommax_e3c, atommax_e3c_idx, &
                                      atommax_e3crel, atommax_e3crel_idx, &
                                      atommax_j, atommax_j_idx)

      CALL ri_rs_topk_reduce_print(bs_env, para_env, particle_set, atomic_kind_set, fit_topk, &
                                   ri_rs_topk_e3c, &
                                   'RI-RS fit error | top-30 (mu nu|P) abs:', 'P')
      CALL ri_rs_topk_reduce_print(bs_env, para_env, particle_set, atomic_kind_set, fit_topk, &
                                   ri_rs_topk_mbare, &
                                   'RI-RS fit error | top-30 3c.Minv bare (mu nu|Q):', 'Q')
      CALL ri_rs_topk_reduce_print(bs_env, para_env, particle_set, atomic_kind_set, fit_topk, &
                                   ri_rs_topk_mtrunc, &
                                   'RI-RS fit error | top-30 3c.Minv trunc (mu nu|Q):', 'Q')
      CALL ri_rs_topk_reduce_print(bs_env, para_env, particle_set, atomic_kind_set, fit_topk, &
                                   ri_rs_topk_e3crel, &
                                   'RI-RS fit error | top-30 (mu nu|P) rel (floor 1.0E-02):', 'P')

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      CALL timestop(handle)

   END SUBROUTINE ri_rs_fit_error_finalize

! **************************************************************************************************
!> \brief Distributed pdpotrf/pdpotrs solve of D x = b for one atom of the
!>        RI-RS Z_lP build. Called when N_PROCS_PER_ATOM_Z_LP > 1, with a
!>        subgroup of cooperating ranks and an associated BLACS context.
!>        Each rank in the subgroup holds the (replicated) phi_local and the
!>        (replicated) RHS d_lp; it builds its own block-cyclic slice of the
!>        squared+Jacobi-scaled Gram matrix D via tiled DGEMM, factorizes via
!>        cp_fm_cholesky_decompose (UPLO='U'), and solves with
!>        cp_fm_cholesky_solve. The replicated d_lp is updated in place via
!>        cp_fm_get_submatrix.
!> \param phi_local   replicated (n_loc x n_ao) AO values on the sphere grid
!> \param d_vec       replicated Jacobi diagonal = 1 / ||phi_i||^2
!> \param d_lp        replicated RHS in/out (n_loc x n_rhs); on input pre-scaled
!>                    by d_vec, on output also pre-scaled (caller post-scales)
!> \param n_loc       number of sphere-local grid points
!> \param n_ao        number of AO basis functions
!> \param n_rhs       number of RI functions of this atom_P (RHS columns)
!> \param tikhonov    regularisation added to the diagonal
!> \param para_env_sub  sub-communicator of the atom-group
!> \param blacs_env_sub BLACS context on the sub-communicator
!> \param fm_struct_D ...
!> \param fm_struct_b ...
!> \param fm_D ...
!> \param fm_b ...
!> \param info        0 on success, non-zero if pdpotrf or pdpotrs failed
! **************************************************************************************************
   SUBROUTINE solve_D_lp_distributed(phi_local, d_vec, d_lp, n_loc, n_ao, n_rhs, &
                                     tikhonov, para_env_sub, blacs_env_sub, &
                                     fm_struct_D, fm_struct_b, fm_D, fm_b, info)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi_local
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: d_vec
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: d_lp
      INTEGER, INTENT(IN)                                :: n_loc, n_ao, n_rhs
      REAL(KIND=dp), INTENT(IN)                          :: tikhonov
      TYPE(mp_para_env_type), POINTER                    :: para_env_sub
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_D, fm_struct_b
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_D, fm_b
      INTEGER, INTENT(OUT)                               :: info

      CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_D_lp_distributed'

      INTEGER                                            :: handle, i_loc, ig, j_loc, jg, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data

      CALL timeset(routineN, handle)
      info = 0

      NULLIFY (fm_struct_D, fm_struct_b)
      CALL cp_fm_struct_create(fm_struct_D, para_env=para_env_sub, &
                               context=blacs_env_sub, &
                               nrow_global=n_loc, ncol_global=n_loc)
      CALL cp_fm_struct_create(fm_struct_b, para_env=para_env_sub, &
                               context=blacs_env_sub, &
                               nrow_global=n_loc, ncol_global=n_rhs)
      CALL cp_fm_create(fm_D, fm_struct_D)
      CALL cp_fm_create(fm_b, fm_struct_b)

      ! ---- Build the local block-cyclic slice of fm_D --------------------
      ! Tiled DGEMM build: for each (row_tile x col_tile) sub-block we
      ! gather small phi strips, DGEMM into a small gram tile, then write
      ! the squared & d_vec-scaled result directly into fm_D%local_data.
      ! This avoids materialising the full (nrow_local x n_ao),
      ! (n_ao x ncol_local), and (nrow_local x ncol_local) buffers.
      CALL cp_fm_get_info(fm_D, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices, &
                          local_data=local_data)

      IF (nrow_local > 0 .AND. ncol_local > 0) THEN
         BLOCK
            INTEGER, PARAMETER :: ntile = 1024
            INTEGER :: ib, ie, jb, je, mb, kb, ti, tj
            REAL(KIND=dp), ALLOCATABLE :: gram_t(:, :), phi_cols_t(:, :), phi_rows_t(:, :)
            ALLOCATE (phi_rows_t(ntile, n_ao), phi_cols_t(n_ao, ntile), gram_t(ntile, ntile))
            DO ib = 1, nrow_local, ntile
               ie = MIN(ib + ntile - 1, nrow_local)
               mb = ie - ib + 1
               !$OMP PARALLEL DO DEFAULT(NONE) &
               !$OMP SHARED(mb, n_ao, phi_rows_t, phi_local, row_indices, ib) &
               !$OMP PRIVATE(ti, j_loc) SCHEDULE(STATIC)
               DO j_loc = 1, n_ao
                  DO ti = 1, mb
                     phi_rows_t(ti, j_loc) = phi_local(row_indices(ib + ti - 1), j_loc)
                  END DO
               END DO
               !$OMP END PARALLEL DO
               DO jb = 1, ncol_local, ntile
                  je = MIN(jb + ntile - 1, ncol_local)
                  kb = je - jb + 1
                  !$OMP PARALLEL DO DEFAULT(NONE) &
                  !$OMP SHARED(kb, n_ao, phi_cols_t, phi_local, col_indices, jb) &
                  !$OMP PRIVATE(tj, i_loc) SCHEDULE(STATIC)
                  DO tj = 1, kb
                     DO i_loc = 1, n_ao
                        phi_cols_t(i_loc, tj) = phi_local(col_indices(jb + tj - 1), i_loc)
                     END DO
                  END DO
                  !$OMP END PARALLEL DO
                  CALL dgemm('N', 'N', mb, kb, n_ao, &
                             1.0_dp, phi_rows_t, ntile, phi_cols_t, n_ao, &
                             0.0_dp, gram_t, ntile)
                  !$OMP PARALLEL DO DEFAULT(NONE) &
                  !$OMP SHARED(mb, kb, gram_t, d_vec, row_indices, col_indices, ib, jb) &
                  !$OMP SHARED(local_data, tikhonov) &
                  !$OMP PRIVATE(ti, tj, ig, jg) SCHEDULE(STATIC)
                  DO tj = 1, kb
                     jg = col_indices(jb + tj - 1)
                     DO ti = 1, mb
                        ig = row_indices(ib + ti - 1)
                        local_data(ib + ti - 1, jb + tj - 1) = &
                           gram_t(ti, tj)*gram_t(ti, tj)*d_vec(ig)*d_vec(jg)
                        IF (ig == jg) THEN
                           local_data(ib + ti - 1, jb + tj - 1) = &
                              local_data(ib + ti - 1, jb + tj - 1) + tikhonov
                        END IF
                     END DO
                  END DO
                  !$OMP END PARALLEL DO
               END DO
            END DO
            DEALLOCATE (phi_rows_t, phi_cols_t, gram_t)
         END BLOCK
      END IF

      ! ---- Load the replicated d_lp into the block-cyclic fm_b -----------
      CALL cp_fm_set_submatrix(fm_b, d_lp)

      ! ---- pdpotrf (Cholesky factorisation; cp_fm_cholesky_decompose
      !      factors with UPLO='U', so pdpotrs must match)
      CALL cp_fm_cholesky_decompose(fm_D, n=n_loc, info_out=info)
      IF (info /= 0) THEN
         CPABORT("pdpotrf failed in solve_D_lp_distributed")
      END IF

      ! ---- pdpotrs/dpotrs (solve in place on fm_b) -----------------------
      CALL cp_fm_cholesky_solve(fm_D, fm_b, n=n_loc, info_out=info)
      IF (info /= 0) THEN
         CPABORT("pdpotrs failed in solve_D_lp_distributed")
      END IF

      ! ---- Gather distributed solution back into the replicated d_lp ----
      CALL cp_fm_get_submatrix(fm_b, d_lp)

      CALL cp_fm_release(fm_D)
      CALL cp_fm_release(fm_b)
      CALL cp_fm_struct_release(fm_struct_D)
      CALL cp_fm_struct_release(fm_struct_b)

      CALL timestop(handle)

   END SUBROUTINE solve_D_lp_distributed

! **************************************************************************************************
!> \brief Computes the χ(iτ) matrix
!> \param bs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
! **************************************************************************************************

   SUBROUTINE get_mat_chi_Gamma_tau(bs_env, mat_chi_Gamma_tau, mat_phi_mu_l, mat_Z_lP)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l, mat_Z_lP

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mat_chi_Gamma_tau'

      INTEGER                                            :: handle, i, i_t, ispin, npcol
      INTEGER, DIMENSION(:), POINTER                     :: blk_ao, blk_grid, dist_col_ao, &
                                                            dist_col_grid, dist_row_grid
      REAL(KIND=dp)                                      :: t1, tau
      TYPE(dbcsr_distribution_type)                      :: dist_grid_grid, dist_phi
      TYPE(dbcsr_type)                                   :: matrix_chi_grid, matrix_chi_grid_spin, &
                                                            matrix_G_occ_grid, matrix_G_vir_grid

      CALL timeset(routineN, handle)

      ! =========================================================================
      ! 1. SETUP CORE TOPOLOGIES
      ! =========================================================================
      CALL dbcsr_get_info(mat_phi_mu_l, distribution=dist_phi, row_blk_size=blk_grid, col_blk_size=blk_ao)
      CALL dbcsr_distribution_get(dist_phi, row_dist=dist_row_grid, col_dist=dist_col_ao)

      ! Determine number of MPI process columns
      npcol = MAXVAL(dist_col_ao) + 1

      ! Build a perfectly safe column distribution for the Grid dimension
      ALLOCATE (dist_col_grid(SIZE(blk_grid)))
      DO i = 1, SIZE(blk_grid)
         dist_col_grid(i) = MOD(i - 1, npcol)
      END DO

      CALL dbcsr_distribution_new(dist_grid_grid, template=dist_phi, &
                                  row_dist=dist_row_grid, col_dist=dist_col_grid)

      CALL dbcsr_create(matrix_G_occ_grid, "G_occ_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)
      CALL dbcsr_create(matrix_G_vir_grid, "G_vir_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)
      CALL dbcsr_create(matrix_chi_grid, "chi_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)
      CALL dbcsr_create(matrix_chi_grid_spin, "chi_grid_spin", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)

      ! =========================================================================
      ! 2. MAIN IMAGINARY TIME LOOP
      ! =========================================================================
      DO i_t = 1, bs_env%num_time_freq_points
         t1 = m_walltime()

         tau = bs_env%imag_time_points(i_t)
         CALL dbcsr_set(matrix_chi_grid, 0.0_dp)

         ! ----------------------------------------------------------------------
         ! A. SPIN LOOP (Allocations safely encapsulated in wrappers)
         ! ----------------------------------------------------------------------
         DO ispin = 1, bs_env%n_spin

            ! G^occ_µλ(i|τ|) = sum_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
            ! G^occ_ll'(i|τ|) = sum_µν Φ_µ(r_l) G^occ_µν Φ_ν(r_l')
            CALL build_G_grid(bs_env, tau, ispin, .TRUE., .FALSE., mat_phi_mu_l, &
                              matrix_G_occ_grid, bs_env%eps_filter)

            ! G^vir_µλ(i|τ|) = sum_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
            ! G^vir_ll'(i|τ|) = sum_µν Φ_µ(r_l) G^vir_µν Φ_ν(r_l')
            CALL build_G_grid(bs_env, tau, ispin, .FALSE., .TRUE., mat_phi_mu_l, &
                              matrix_G_vir_grid, bs_env%eps_filter)

            ! -------------------------------------------------------------------
            ! B. ELEMENT-WISE HADAMARD PRODUCT
            ! -------------------------------------------------------------------
            ! χ_ll'(iτ) = G^occ_ll'(i|τ|) * G^vir_ll'(i|τ|)
            CALL hadamard_product(matrix_G_occ_grid, matrix_G_vir_grid, matrix_chi_grid_spin, bs_env%spin_degeneracy)

            ! Accumulate spin contributions
            CALL dbcsr_add(matrix_chi_grid, matrix_chi_grid_spin, 1.0_dp, 1.0_dp)

         END DO ! ispin

         ! ----------------------------------------------------------------------
         ! C. TRANSFORM TO AUXILIARY BASIS & EXPORT DIRECTLY
         ! χ_aux        = Z^T * χ_grid * Z
         ! χ_PQ(iτ) = sum_ll' Z_lP χ_ll'(iτ) Z_l'Q
         ! Result is dumped into the final array mat_chi_Gamma_tau!
         ! ----------------------------------------------------------------------
         CALL contract_A_B_A("T", "N", mat_Z_lP, matrix_chi_grid, &
                             mat_chi_Gamma_tau(i_t)%matrix, bs_env%eps_filter)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I13,A,I3,A,F7.1,A)') &
               'Computed χ(iτ,k=0) for time point', i_t, ' /', bs_env%num_time_freq_points, &
               ',    Execution time', m_walltime() - t1, ' s'
         END IF

      END DO ! i_t

      ! =========================================================================
      ! 3. FINAL CLEANUP
      ! =========================================================================
      CALL dbcsr_release(matrix_G_occ_grid)
      CALL dbcsr_release(matrix_G_vir_grid)
      CALL dbcsr_release(matrix_chi_grid)
      CALL dbcsr_release(matrix_chi_grid_spin)
      CALL dbcsr_distribution_release(dist_grid_grid)
      DEALLOCATE (dist_col_grid)

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE get_mat_chi_Gamma_tau

! **************************************************************************************************
!> \brief Computes Green's Function in grid basis
!> \param bs_env ...
!> \param tau ...
!> \param ispin ...
!> \param occ ...
!> \param vir ...
!> \param mat_phi_mu_l ...
!> \param matrix_G_grid ...
!> \param eps_filter ...
! **************************************************************************************************

   SUBROUTINE build_G_grid(bs_env, tau, ispin, occ, vir, mat_phi_mu_l, matrix_G_grid, eps_filter)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: tau
      INTEGER, INTENT(IN)                                :: ispin
      LOGICAL, INTENT(IN)                                :: occ, vir
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l, matrix_G_grid
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_G_grid'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: blk_ao, dist_row_ao
      TYPE(cp_fm_type), POINTER                          :: fm_G
      TYPE(dbcsr_distribution_type)                      :: dist_ao_ao
      TYPE(dbcsr_type)                                   :: matrix_G_ao

      CALL timeset(routineN, handle)

      ! 1. Select the FM matrix based on occ/vir flags
      IF (occ) THEN
         fm_G => bs_env%fm_Gocc
      ELSE
         fm_G => bs_env%fm_Gvir
      END IF

      ! 2. Compute Dense FM Green's Function
      ! G^occ/vir_µλ(i|τ|) = sum_G^occ/vir_µλn^occ/vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
      CALL G_occ_vir(bs_env, tau, fm_G, ispin, occ=occ, vir=vir)

      ! 3. Setup AO DBCSR Topology and Create Matrix dynamically
      CALL setup_square_topology(mat_phi_mu_l, 'COL', dist_ao_ao, blk_ao, dist_row_ao)

      CALL dbcsr_create(matrix_G_ao, name="G_ao", dist=dist_ao_ao, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=blk_ao, col_blk_size=blk_ao)

      ! 4. Convert FM to Sparse DBCSR
      CALL copy_fm_to_dbcsr(fm_G, matrix_G_ao, keep_sparsity=.FALSE.)

      ! 5. Transform to Grid Basis: G_grid = phi * G_ao * phi^T
      ! G^occ/vir_ll'(i|τ|) = sum_µν Φ_µ(r_l) G^occ/vir_µν Φ_ν(r_l')
      CALL contract_A_B_A("N", "T", mat_phi_mu_l, matrix_G_ao, matrix_G_grid, eps_filter)

      ! 6. Release AO matrix and topology
      CALL release_dbcsr_topology_and_matrices(dist=dist_ao_ao, mapped_dist=dist_row_ao, m1=matrix_G_ao)

      CALL timestop(handle)

   END SUBROUTINE build_G_grid

! **************************************************************************************************
!> \brief Generalized routine to compute OUT = A * B * A^T  OR  OUT = A^T * B * A using DBCSR
!> \param transA_left ...
!> \param transA_right ...
!> \param matrix_A ...
!> \param matrix_B ...
!> \param matrix_out ...
!> \param eps_filter ...
! **************************************************************************************************

   SUBROUTINE contract_A_B_A(transA_left, transA_right, matrix_A, matrix_B, matrix_out, eps_filter)

      CHARACTER(LEN=1), INTENT(IN)                       :: transA_left, transA_right
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_A, matrix_B, matrix_out
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract_A_B_A'

      INTEGER                                            :: handle
      TYPE(dbcsr_type)                                   :: matrix_tmp

      CALL timeset(routineN, handle)

      CALL dbcsr_create(matrix_tmp, template=matrix_A)

      IF (transA_left == "N" .AND. transA_right == "T") THEN
         ! Path 1: Out = A * B * A^T
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_A, matrix_B, &
                             0.0_dp, matrix_tmp, filter_eps=eps_filter)
         CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp, matrix_A, &
                             0.0_dp, matrix_out, filter_eps=eps_filter)

      ELSE IF (transA_left == "T" .AND. transA_right == "N") THEN
         ! Path 2: Out = A^T * B * A
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_B, matrix_A, &
                             0.0_dp, matrix_tmp, filter_eps=eps_filter)
         CALL dbcsr_multiply("T", "N", 1.0_dp, matrix_A, matrix_tmp, &
                             0.0_dp, matrix_out, filter_eps=eps_filter)
      ELSE
         CPABORT("Unsupported transposition pair in contract_A_B_A")
      END IF

      CALL dbcsr_release(matrix_tmp)

      CALL timestop(handle)

   END SUBROUTINE contract_A_B_A

! **************************************************************************************************
!> \brief Computes C = A ◦ B (Element-wise Hadamard product) for sparse DBCSR matrices.
!> \param matrix_A ...
!> \param matrix_B ...
!> \param matrix_C ...
!> \param fac (Scaling factor applied to the product)
! **************************************************************************************************

   SUBROUTINE hadamard_product(matrix_A, matrix_B, matrix_C, fac)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_A, matrix_B, matrix_C
      REAL(KIND=dp), INTENT(IN)                          :: fac

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'hadamard_product'

      INTEGER                                            :: col, handle, row
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_B, blk_C
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_copy(matrix_C, matrix_A)

      CALL dbcsr_iterator_start(iter, matrix_C)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, blk_C)

         CALL dbcsr_get_block_p(matrix_B, row, col, blk_B, found)

         IF (found) THEN
            blk_C(:, :) = fac*blk_C(:, :)*blk_B(:, :)
         ELSE
            ! If B is sparse here, the product is zero
            blk_C(:, :) = 0.0_dp
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE hadamard_product

! **************************************************************************************************
!> \brief Compute screened Coulomb interaction matrix
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_time ...
! **************************************************************************************************

   SUBROUTINE compute_W(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_W'

      INTEGER                                            :: handle, i_t, j_w
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_fm_type)                                   :: fm_M_inv_V_sqrt, fm_V, fm_V_sqrt

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      CALL create_fm_W_MIC_time(bs_env, fm_W_time)

      ! 1. Allocate V and M matrices
      CALL cp_fm_create(fm_V, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(fm_V_sqrt, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(fm_M_inv_V_sqrt, bs_env%fm_RI_RI%matrix_struct)

      ! Compute V and M^-1 * V^0.5
      CALL compute_V_MinvVsqrt(bs_env, qs_env, fm_V, fm_V_sqrt, fm_M_inv_V_sqrt)

      ! 2. Loop over frequencies
      DO j_w = 1, bs_env%num_time_freq_points
         ! Fourier transformation of χ_PQ(iτ) to χ_PQ(iω_j)
         CALL compute_fm_chi_Gamma_freq(bs_env, bs_env%fm_chi_Gamma_freq, j_w, mat_chi_Gamma_tau)

         ! ε(iω_j) = Id - V^0.5*M^-1*χ(iω_j)*M^-1*V^0.5
         ! W(iω_j) = V^0.5*(ε^-1(iω_j)-Id)*V^0.5
         CALL compute_fm_W_freq(bs_env, bs_env%fm_chi_Gamma_freq, fm_V_sqrt, &
                                fm_M_inv_V_sqrt, bs_env%fm_W_MIC_freq)

         ! Fourier transform from W_PQ^MIC(iω_j) to W_PQ^MIC(iτ)
         CALL Fourier_transform_w_to_t(bs_env, fm_W_time, bs_env%fm_W_MIC_freq, j_w)
      END DO

      ! M^-1*W^MIC(iτ)*M^-1
      CALL multiply_fm_W_MIC_time_with_Minv_Gamma(bs_env, qs_env, fm_W_time)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T55,A,F10.1,A)') &
            'Computed W(iτ),', ' Execution time', m_walltime() - t1, ' s'
      END IF

      CALL dbcsr_deallocate_matrix_set(mat_chi_Gamma_tau)

      ! Cleanup
      CALL cp_fm_release(fm_V)
      CALL cp_fm_release(fm_V_sqrt)
      CALL cp_fm_release(fm_M_inv_V_sqrt)

      ! Marek : Fourier transform W^MIC(itau) back to get it at a specific im.frequency point - iomega = 0
      IF (bs_env%rtp_method == rtp_method_bse) THEN
         t1 = m_walltime()
         CALL cp_fm_create(bs_env%fm_W_MIC_freq_zero, bs_env%fm_W_MIC_freq%matrix_struct)
         ! Set to zero
         CALL cp_fm_set_all(bs_env%fm_W_MIC_freq_zero, 0.0_dp)
         ! Sum over all times
         DO i_t = 1, bs_env%num_time_freq_points
            ! Add the relevant structure with correct weight
            CALL cp_fm_scale_and_add(1.0_dp, bs_env%fm_W_MIC_freq_zero, &
                                     bs_env%imag_time_weights_freq_zero(i_t), fm_W_time(i_t))
         END DO
         ! Done, save to file
         CALL fm_write(bs_env%fm_W_MIC_freq_zero, 0, "W_freq_rtp", qs_env)
         ! Report calculation
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T55,A,F10.1,A)') &
               'Computed W(0),', ' Execution time', m_walltime() - t1, ' s'
         END IF
      END IF

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE compute_W

! **************************************************************************************************
!> \brief Computes V, V^0.5, and M^-1 * V^0.5
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_V ...
!> \param fm_V_sqrt ...
!> \param fm_Minv_Vsqrt ...
! **************************************************************************************************

   SUBROUTINE compute_V_MinvVsqrt(bs_env, qs_env, fm_V, fm_V_sqrt, fm_Minv_Vsqrt)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_V, fm_V_sqrt, fm_Minv_Vsqrt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_V_MinvVsqrt'

      INTEGER                                            :: handle, info, n_RI, ndep
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type)                                   :: fm_work
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_M
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_V_kp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI
      CALL cp_fm_create(fm_work, fm_V%matrix_struct)

      ! -----------------------------------------------------------------------
      ! 1. Build Coulomb Matrix V(k=0) using the kp-routine but only for ikp=1
      ! -----------------------------------------------------------------------
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell, &
                      qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set)

      ALLOCATE (mat_V_kp(1:1, 1:2))
      NULLIFY (mat_V_kp(1, 1)%matrix, mat_V_kp(1, 2)%matrix)
      ALLOCATE (mat_V_kp(1, 1)%matrix, mat_V_kp(1, 2)%matrix)

      CALL dbcsr_create(mat_V_kp(1, 1)%matrix, template=bs_env%mat_RI_RI%matrix)
      CALL dbcsr_reserve_all_blocks(mat_V_kp(1, 1)%matrix)
      CALL dbcsr_set(mat_V_kp(1, 1)%matrix, 0.0_dp)

      ! Dummy imaginary part just to satisfy the routine
      CALL dbcsr_create(mat_V_kp(1, 2)%matrix, template=bs_env%mat_RI_RI%matrix)
      CALL dbcsr_reserve_all_blocks(mat_V_kp(1, 2)%matrix)
      CALL dbcsr_set(mat_V_kp(1, 2)%matrix, 0.0_dp)

      bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig

      CALL build_2c_coulomb_matrix_kp(mat_V_kp, bs_env%kpoints_chi_eps_W, "RI_AUX", cell, &
                                      particle_set, qs_kind_set, atomic_kind_set, &
                                      bs_env%size_lattice_sum_V, operator_coulomb, 1, 1)

      ! Copy real part to fm_V
      CALL copy_dbcsr_to_fm(mat_V_kp(1, 1)%matrix, fm_V)

      CALL dbcsr_deallocate_matrix(mat_V_kp(1, 1)%matrix)
      CALL dbcsr_deallocate_matrix(mat_V_kp(1, 2)%matrix)
      DEALLOCATE (mat_V_kp)

      ! -----------------------------------------------------------------------
      ! 2. Get RI-Metric Matrix M(k=0)
      ! -----------------------------------------------------------------------
      CALL RI_2c_integral_mat(qs_env, fm_M, fm_V, n_RI, bs_env%ri_metric, &
                              do_kpoints=.FALSE., regularization_RI=bs_env%regularization_RI)

      ! -----------------------------------------------------------------------
      ! 3. M -> M^-1
      ! -----------------------------------------------------------------------
      CALL cp_fm_cholesky_decompose(fm_M(1, 1), info_out=info)
      IF (info == 0) THEN
         CALL cp_fm_cholesky_invert(fm_M(1, 1))
         CALL cp_fm_uplo_to_full(fm_M(1, 1), fm_work)
      ELSE
         ! Fallback if Cholesky fails due to conditioning
         CALL cp_fm_power(fm_M(1, 1), fm_work, -1.0_dp, bs_env%eps_eigval_mat_RI, ndep)
         CALL cp_fm_to_fm(fm_work, fm_M(1, 1))
      END IF

      ! -----------------------------------------------------------------------
      ! 4. V -> V^0.5
      ! -----------------------------------------------------------------------
      CALL cp_fm_to_fm(fm_V, fm_V_sqrt)
      CALL cp_fm_cholesky_decompose(fm_V_sqrt, info_out=info)
      IF (info == 0) THEN
         CALL clean_lower_part(fm_V_sqrt)
      ELSE
         CALL cp_fm_power(fm_V, fm_V_sqrt, 0.5_dp, bs_env%eps_eigval_mat_RI, ndep)
      END IF

      ! -----------------------------------------------------------------------
      ! 5. M^-1 * V^0.5
      ! -----------------------------------------------------------------------
      CALL parallel_gemm("N", "T", n_RI, n_RI, n_RI, 1.0_dp, fm_M(1, 1), fm_V_sqrt, &
                         0.0_dp, fm_Minv_Vsqrt)

      CALL cp_fm_release(fm_M)
      CALL cp_fm_release(fm_work)

      CALL timestop(handle)

   END SUBROUTINE compute_V_MinvVsqrt

! **************************************************************************************************
!> \brief Computes W(iω) from χ_PQ(iω_j)
!> \param bs_env ...
!> \param fm_chi_freq_j ...
!> \param fm_V_sqrt ...
!> \param fm_Minv_Vsqrt ...
!> \param fm_W_freq_j ...
! **************************************************************************************************

   SUBROUTINE compute_fm_W_freq(bs_env, fm_chi_freq_j, fm_V_sqrt, fm_Minv_Vsqrt, fm_W_freq_j)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_chi_freq_j, fm_V_sqrt, fm_Minv_Vsqrt
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_W_freq_j

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_fm_W_freq'

      INTEGER                                            :: handle, info, n_RI, ndep
      TYPE(cp_fm_type)                                   :: fm_eps_freq_j, fm_work

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI

      CALL cp_fm_create(fm_eps_freq_j, fm_chi_freq_j%matrix_struct)
      CALL cp_fm_create(fm_work, fm_chi_freq_j%matrix_struct)

      ! -----------------------------------------------------------------------
      ! 1. ε(iω_j) = Id - (M^-1 * V^0.5)^T * χ(iω_j) * (M^-1 * V^0.5)
      ! -----------------------------------------------------------------------
      ! work = χ(iω_j) * (M^-1 * V^0.5)
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, &
                         fm_chi_freq_j, fm_Minv_Vsqrt, 0.0_dp, fm_work)

      ! eps_work = (M^-1 * V^0.5)^T * work
      CALL parallel_gemm('T', 'N', n_RI, n_RI, n_RI, 1.0_dp, &
                         fm_Minv_Vsqrt, fm_work, 0.0_dp, fm_eps_freq_j)

      ! ε(iω_j) = Id - eps_work  -->  -eps_work + Id
      CALL fm_add_on_diag(fm_eps_freq_j, 1.0_dp)

      ! Force perfect symmetry before Cholesky to avoid info != 0 due to GEMM noise
      CALL cp_fm_uplo_to_full(fm_eps_freq_j, fm_work)

      ! -----------------------------------------------------------------------
      ! 2. W(iω_j) = V^0.5^T * (ε^-1(iω_j) - Id) * V^0.5
      ! -----------------------------------------------------------------------

      ! a) Cholesky decomposition of ε(iω_j)
      CALL cp_fm_cholesky_decompose(fm_eps_freq_j, info_out=info)

      ! b) Inversion
      IF (info == 0) THEN
         CALL cp_fm_cholesky_invert(fm_eps_freq_j)
         CALL cp_fm_uplo_to_full(fm_eps_freq_j, fm_work)
      ELSE
         ! Fallback to expensive diagonalization if Cholesky fails
         CALL cp_fm_power(fm_eps_freq_j, fm_work, -1.0_dp, bs_env%eps_eigval_mat_RI, ndep)
         CALL cp_fm_to_fm(fm_work, fm_eps_freq_j)
      END IF

      ! c) ε^-1(iω_j) - Id
      CALL fm_add_on_diag(fm_eps_freq_j, -1.0_dp)

      ! d) work = (ε^-1(iω_j) - Id) * V^0.5
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, fm_eps_freq_j, fm_V_sqrt, &
                         0.0_dp, fm_work)

      ! e) W(iw) = V^0.5^T * work
      CALL parallel_gemm('T', 'N', n_RI, n_RI, n_RI, 1.0_dp, fm_V_sqrt, fm_work, &
                         0.0_dp, fm_W_freq_j)

      ! Cleanup
      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(fm_eps_freq_j)

      CALL timestop(handle)

   END SUBROUTINE compute_fm_W_freq

! **************************************************************************************************
!> \brief Adds a real scalar value to the diagonal of a real full matrix (fm)
!> \param fm ...
!> \param alpha ...
! **************************************************************************************************

   SUBROUTINE fm_add_on_diag(fm, alpha)
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm
      REAL(KIND=dp), INTENT(IN)                          :: alpha

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fm_add_on_diag'

      INTEGER                                            :: handle, i_global, i_row, j_col, &
                                                            j_global, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=fm, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      DO j_col = 1, ncol_local
         j_global = col_indices(j_col)
         DO i_row = 1, nrow_local
            i_global = row_indices(i_row)
            IF (j_global == i_global) THEN
               fm%local_data(i_row, j_col) = fm%local_data(i_row, j_col) + alpha
            END IF
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE fm_add_on_diag

! **************************************************************************************************
!> \brief Zeroes out the strictly lower triangular part of a real matrix
!> \param fm_mat ...
! **************************************************************************************************
   SUBROUTINE clean_lower_part(fm_mat)
      TYPE(cp_fm_type)                                   :: fm_mat

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'clean_lower_part'

      INTEGER                                            :: handle, i_row, j_col, j_global, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=fm_mat, &
                          nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)

      DO j_col = 1, ncol_local
         j_global = col_indices(j_col)
         DO i_row = 1, nrow_local
            IF (j_global < row_indices(i_row)) fm_mat%local_data(i_row, j_col) = 0.0_dp
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE clean_lower_part

! **************************************************************************************************
!> \brief Computes the exact exchange part of the GW self-energy
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
!> \param fm_Sigma_x_Gamma ...
! **************************************************************************************************

   SUBROUTINE compute_Sigma_x(bs_env, qs_env, mat_phi_mu_l, mat_Z_lP, fm_Sigma_x_Gamma)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l, mat_Z_lP
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_x'

      INTEGER                                            :: handle, ispin
      INTEGER, DIMENSION(:), POINTER                     :: blk_aux, blk_grid, dist_col_grid, &
                                                            dist_row_aux
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Vtr_Gamma
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux, dist_grid_grid
      TYPE(dbcsr_type)                                   :: mat_Sigma_x_Gamma, matrix_D_grid, &
                                                            matrix_Sigma_x_grid, matrix_V_aux, &
                                                            matrix_V_grid

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      ALLOCATE (fm_Sigma_x_Gamma(bs_env%n_spin))
      DO ispin = 1, bs_env%n_spin
         CALL cp_fm_create(fm_Sigma_x_Gamma(ispin), bs_env%fm_s_Gamma%matrix_struct)
      END DO

      CALL dbcsr_create(mat_Sigma_x_Gamma, template=bs_env%mat_ao_ao%matrix)

      ! =========================================================================
      ! 1. SETUP CORE TOPOLOGIES
      ! =========================================================================
      CALL setup_square_topology(mat_phi_mu_l, 'ROW', dist_grid_grid, blk_grid, dist_col_grid)
      CALL setup_square_topology(mat_Z_lP, 'COL', dist_aux_aux, blk_aux, dist_row_aux)

      ! =========================================================================
      ! 2. COMPUTE V^tr_ll'
      ! =========================================================================
      CALL RI_2c_integral_mat(qs_env, fm_Vtr_Gamma, bs_env%fm_RI_RI, bs_env%n_RI, &
                              bs_env%trunc_coulomb, do_kpoints=.FALSE.)

      ! M^-1 * V^tr * M^-1 directly modifies fm_Vtr_Gamma(:, 1)
      CALL multiply_fm_W_MIC_time_with_Minv_Gamma(bs_env, qs_env, fm_Vtr_Gamma(:, 1))

      CALL dbcsr_create(matrix_V_aux, "V_aux", dist_aux_aux, dbcsr_type_no_symmetry, blk_aux, blk_aux)
      CALL dbcsr_create(matrix_V_grid, "V_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)

      CALL copy_fm_to_dbcsr(fm_Vtr_Gamma(1, 1), matrix_V_aux, keep_sparsity=.FALSE.)

      ! V^tr_ll' = sum_PQ Z_lP V^trunc_PQ Z_l'Q
      CALL contract_A_B_A("N", "T", mat_Z_lP, matrix_V_aux, matrix_V_grid, bs_env%eps_filter)
      CALL dbcsr_release(matrix_V_aux)

      ! =========================================================================
      ! 3. SPIN LOOP FOR EXACT EXCHANGE
      ! =========================================================================
      DO ispin = 1, bs_env%n_spin

         ! Density matrix on grid is essentially G_occ at tau = 0.0
         ! D_µν  = sum_n^occ C_µn C_νn
         ! D_ll' = sum_µν Φ_µ(r_l) D_µν Φ_ν(r_l')
         CALL dbcsr_create(matrix_D_grid, "D_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)
         CALL build_G_grid(bs_env, 0.0_dp, ispin, .TRUE., .FALSE., mat_phi_mu_l, matrix_D_grid, bs_env%eps_filter)

         ! Element-wise Hadamard product: Σ^x_grid = D_grid ◦ V_grid
         ! Σ^x_ll' = D_ll' * V^tr_ll'
         CALL dbcsr_create(matrix_Sigma_x_grid, template=matrix_V_grid)
         CALL hadamard_product(matrix_D_grid, matrix_V_grid, matrix_Sigma_x_grid, 1.0_dp)

         CALL dbcsr_release(matrix_D_grid)

         ! Transform back to AO basis: Σ^x_ao = -1.0 * phi^T * Σ^x_grid * phi
         ! Σ^x_λσ   = -sum_ll' Φ_λ(r_l) Σ^x_ll' Φ_σ(r_l')
         CALL contract_A_B_A("T", "N", mat_phi_mu_l, matrix_Sigma_x_grid, mat_Sigma_x_Gamma, bs_env%eps_filter)
         CALL dbcsr_scale(mat_Sigma_x_Gamma, -1.0_dp)

         CALL dbcsr_release(matrix_Sigma_x_grid)

         ! Data I/O and Export to CP2K Full Matrices
         CALL copy_dbcsr_to_fm(mat_Sigma_x_Gamma, fm_Sigma_x_Gamma(ispin))

      END DO ! ispin

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Computed Σ^x(k=0),', ' Execution time', m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      ! =========================================================================
      ! 4. CLEANUP
      ! =========================================================================
      CALL release_dbcsr_topology_and_matrices(dist=dist_grid_grid, mapped_dist=dist_col_grid, &
                                               m1=mat_Sigma_x_Gamma, m2=matrix_V_grid)
      CALL release_dbcsr_topology_and_matrices(dist=dist_aux_aux, mapped_dist=dist_row_aux)

      CALL cp_fm_release(fm_Vtr_Gamma)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_x

! **************************************************************************************************
!> \brief Computes the correlation part of the GW self-energy
!> \param bs_env ...
!> \param fm_W_time ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
!> \param fm_Sigma_c_Gamma_time ...
! **************************************************************************************************

   SUBROUTINE compute_Sigma_c(bs_env, fm_W_time, mat_phi_mu_l, mat_Z_lP, fm_Sigma_c_Gamma_time)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_time
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l, mat_Z_lP
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_c'

      INTEGER                                            :: handle, i_t, ispin
      INTEGER, DIMENSION(:), POINTER                     :: blk_aux, blk_grid, dist_col_grid, &
                                                            dist_row_aux
      REAL(KIND=dp)                                      :: t1, tau
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux, dist_grid_grid
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_Sigma_neg_tau, mat_Sigma_pos_tau
      TYPE(dbcsr_type) :: matrix_G_occ_grid, matrix_G_vir_grid, matrix_Sigma_neg_grid, &
         matrix_Sigma_pos_grid, matrix_W_aux, matrix_W_grid

      CALL timeset(routineN, handle)

      ! =========================================================================
      ! 1. SETUP CORE TOPOLOGIES AND PRE-ALLOCATE OUTPUT ARRAYS
      ! =========================================================================
      CALL setup_square_topology(mat_phi_mu_l, 'ROW', dist_grid_grid, blk_grid, dist_col_grid)
      CALL setup_square_topology(mat_Z_lP, 'COL', dist_aux_aux, blk_aux, dist_row_aux)

      ! Pre-allocate local DBCSR matrices to act as targets for final output
      NULLIFY (mat_Sigma_neg_tau, mat_Sigma_pos_tau)
      ALLOCATE (mat_Sigma_neg_tau(bs_env%num_time_freq_points, bs_env%n_spin))
      ALLOCATE (mat_Sigma_pos_tau(bs_env%num_time_freq_points, bs_env%n_spin))

      DO i_t = 1, bs_env%num_time_freq_points
         DO ispin = 1, bs_env%n_spin
            ALLOCATE (mat_Sigma_neg_tau(i_t, ispin)%matrix)
            ALLOCATE (mat_Sigma_pos_tau(i_t, ispin)%matrix)
            CALL dbcsr_create(mat_Sigma_neg_tau(i_t, ispin)%matrix, template=bs_env%mat_ao_ao%matrix)
            CALL dbcsr_create(mat_Sigma_pos_tau(i_t, ispin)%matrix, template=bs_env%mat_ao_ao%matrix)
         END DO
      END DO

      ! =========================================================================
      ! 2. MAIN IMAGINARY TIME LOOP
      ! =========================================================================
      DO i_t = 1, bs_env%num_time_freq_points
         tau = bs_env%imag_time_points(i_t)

         ! -------------------------------------------------------------------
         ! Compute W_grid = Z * W_aux * Z^T
         ! -------------------------------------------------------------------
         CALL dbcsr_create(matrix_W_aux, "W_aux", dist_aux_aux, dbcsr_type_no_symmetry, blk_aux, blk_aux)
         CALL dbcsr_create(matrix_W_grid, "W_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)

         CALL copy_fm_to_dbcsr(fm_W_time(i_t), matrix_W_aux, keep_sparsity=.FALSE.)

         ! W^MIC_ll'(iτ,k=0) = sum_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
         CALL contract_A_B_A("N", "T", mat_Z_lP, matrix_W_aux, matrix_W_grid, bs_env%eps_filter)

         CALL dbcsr_release(matrix_W_aux) ! Clean up aux basis immediately

         DO ispin = 1, bs_env%n_spin
            t1 = m_walltime()

            ! -------------------------------------------------------------------
            ! A. Transform Green's Functions to the Grid
            ! -------------------------------------------------------------------
            CALL dbcsr_create(matrix_G_occ_grid, "G_occ_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)
            CALL dbcsr_create(matrix_G_vir_grid, "G_vir_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)

            ! G^occ_µλ(i|τ|) = sum_G^occ_µλn^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
            ! G^occ_ll'(i|τ|) = sum_µν Φ_µ(r_l) G^occ_µν Φ_ν(r_l')
            CALL build_G_grid(bs_env, tau, ispin, .TRUE., .FALSE., mat_phi_mu_l, matrix_G_occ_grid, bs_env%eps_filter)

            ! G^vir_µλ(i|τ|) = sum_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
            ! G^vir_ll'(i|τ|) = sum_µν Φ_µ(r_l) G^vir_µν Φ_ν(r_l')
            CALL build_G_grid(bs_env, tau, ispin, .FALSE., .TRUE., mat_phi_mu_l, matrix_G_vir_grid, bs_env%eps_filter)

            ! -------------------------------------------------------------------
            ! B. Element-wise Hadamard Products for Sigma_c on Grid
            ! Σ_neg_grid = G_occ_grid ◦ W_grid
            ! Σ_pos_grid = G_vir_grid ◦ W_grid
            ! -------------------------------------------------------------------
            CALL dbcsr_create(matrix_Sigma_neg_grid, template=matrix_W_grid)
            CALL dbcsr_create(matrix_Sigma_pos_grid, template=matrix_W_grid)

            ! Σ^c_ll'(iτ)   = -G^occ_ll'(i|τ|) * W^MIC_ll'(iτ), for τ < 0
            CALL hadamard_product(matrix_G_occ_grid, matrix_W_grid, matrix_Sigma_neg_grid, 1.0_dp)

            ! Σ^c_ll'(iτ)   =  G^vir_ll'(i|τ|) * W^MIC_ll'(iτ), for τ > 0
            CALL hadamard_product(matrix_G_vir_grid, matrix_W_grid, matrix_Sigma_pos_grid, 1.0_dp)

            ! Instantly purge massive G_grid arrays to save memory
            CALL dbcsr_release(matrix_G_occ_grid)
            CALL dbcsr_release(matrix_G_vir_grid)

            ! -------------------------------------------------------------------
            ! C. Transform Sigma back to AO Basis
            ! Σ_AO = phi^T * Σ_grid * phi
            ! -------------------------------------------------------------------

            ! Σ^c_λσ(iτ)    = sum_ll' Φ_λ(r_l) Σ^c_ll'(iτ) Φ_σ(r_l'), for τ < 0
            CALL contract_A_B_A("T", "N", mat_phi_mu_l, matrix_Sigma_neg_grid, &
                                mat_Sigma_neg_tau(i_t, ispin)%matrix, bs_env%eps_filter)
            CALL dbcsr_scale(mat_Sigma_neg_tau(i_t, ispin)%matrix, -1.0_dp)

            ! Σ^c_λσ(iτ)    = sum_ll' Φ_λ(r_l) Σ^c_ll'(iτ) Φ_σ(r_l'), for τ > 0
            CALL contract_A_B_A("T", "N", mat_phi_mu_l, matrix_Sigma_pos_grid, &
                                mat_Sigma_pos_tau(i_t, ispin)%matrix, bs_env%eps_filter)

            ! Purge Grid Sigma arrays
            CALL dbcsr_release(matrix_Sigma_neg_grid)
            CALL dbcsr_release(matrix_Sigma_pos_grid)

            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,I10,A,I3,A,F7.1,A)') &
                  'Computed Σ^c(iτ) for time point ', i_t, ' /', bs_env%num_time_freq_points, &
                  ',    Execution time', m_walltime() - t1, ' s'
            END IF

         END DO ! ispin

         ! Release the W_grid for this time point
         CALL dbcsr_release(matrix_W_grid)

      END DO ! i_t

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      ! -------------------------------------------------------------------------
      ! 3. FINALIZE AND CLEANUP
      ! -------------------------------------------------------------------------
      CALL fill_fm_Sigma_c_Gamma_time(fm_Sigma_c_Gamma_time, bs_env, &
                                      mat_Sigma_pos_tau, mat_Sigma_neg_tau)

      CALL cp_fm_release(fm_W_time)

      CALL dbcsr_deallocate_matrix_set(mat_Sigma_neg_tau)
      CALL dbcsr_deallocate_matrix_set(mat_Sigma_pos_tau)

      CALL release_dbcsr_topology_and_matrices(dist=dist_grid_grid, mapped_dist=dist_col_grid)
      CALL release_dbcsr_topology_and_matrices(dist=dist_aux_aux, mapped_dist=dist_row_aux)

      CALL delete_unnecessary_files(bs_env)
      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c

! **************************************************************************************************
!> \brief DBCSR Topology Generation
!> \param matrix_template ...
!> \param dim_type ...
!> \param square_dist ...
!> \param blk_sizes ...
!> \param mapped_dist ...
! **************************************************************************************************

   SUBROUTINE setup_square_topology(matrix_template, dim_type, square_dist, blk_sizes, mapped_dist)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_template
      CHARACTER(LEN=*), INTENT(IN)                       :: dim_type
      TYPE(dbcsr_distribution_type), INTENT(OUT)         :: square_dist
      INTEGER, DIMENSION(:), INTENT(OUT), POINTER        :: blk_sizes, mapped_dist

      INTEGER                                            :: i, np
      INTEGER, DIMENSION(:), POINTER                     :: col_blk, col_dist, row_blk, row_dist
      TYPE(dbcsr_distribution_type)                      :: dist_template

      CALL dbcsr_get_info(matrix_template, distribution=dist_template, &
                          row_blk_size=row_blk, col_blk_size=col_blk)
      CALL dbcsr_distribution_get(dist_template, row_dist=row_dist, col_dist=col_dist)

      IF (TRIM(dim_type) == 'ROW') THEN
         ! Creates ROW x ROW (e.g., Grid x Grid from mat_phi_mu_l)
         blk_sizes => row_blk
         np = MAXVAL(col_dist) + 1  ! npcol
         ALLOCATE (mapped_dist(SIZE(blk_sizes)))
         DO i = 1, SIZE(blk_sizes)
            mapped_dist(i) = MOD(i - 1, np)
         END DO
         CALL dbcsr_distribution_new(square_dist, template=dist_template, &
                                     row_dist=row_dist, col_dist=mapped_dist)

      ELSE IF (TRIM(dim_type) == 'COL') THEN
         ! Creates COL x COL (e.g., Aux x Aux from mat_Z_lP)
         blk_sizes => col_blk
         np = MAXVAL(row_dist) + 1  ! nprow
         ALLOCATE (mapped_dist(SIZE(blk_sizes)))
         DO i = 1, SIZE(blk_sizes)
            mapped_dist(i) = MOD(i - 1, np)
         END DO
         CALL dbcsr_distribution_new(square_dist, template=dist_template, &
                                     row_dist=mapped_dist, col_dist=col_dist)
      END IF

   END SUBROUTINE setup_square_topology

! **************************************************************************************************
!> \brief DBCSR matrices deallocation
!> \param dist ...
!> \param mapped_dist    ...
!> \param m1 ...
!> \param m2 ...
!> \param m3 ...
!> \param m4 ...
! **************************************************************************************************

   SUBROUTINE release_dbcsr_topology_and_matrices(dist, mapped_dist, m1, m2, m3, m4)

      TYPE(dbcsr_distribution_type), INTENT(INOUT), &
         OPTIONAL                                        :: dist
      INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL, &
         POINTER                                         :: mapped_dist
      TYPE(dbcsr_type), INTENT(INOUT), OPTIONAL          :: m1, m2, m3, m4

      IF (PRESENT(dist)) CALL dbcsr_distribution_release(dist)
      IF (PRESENT(mapped_dist)) THEN
         IF (ASSOCIATED(mapped_dist)) THEN
            DEALLOCATE (mapped_dist)
            NULLIFY (mapped_dist)
         END IF
      END IF
      IF (PRESENT(m1)) CALL dbcsr_release(m1)
      IF (PRESENT(m2)) CALL dbcsr_release(m2)
      IF (PRESENT(m3)) CALL dbcsr_release(m3)
      IF (PRESENT(m4)) CALL dbcsr_release(m4)

   END SUBROUTINE release_dbcsr_topology_and_matrices

END MODULE gw_non_periodic_ri_rs
