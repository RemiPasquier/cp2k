!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Work with symmetry
!> \par History
!> \author jgh
! **************************************************************************************************
MODULE cp_symmetry
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              real_to_scaled
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cryssym,                         ONLY: crys_sym_gen,&
                                              csym_type,&
                                              print_crys_symmetry,&
                                              release_csym_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE molsym,                          ONLY: molecular_symmetry,&
                                              molsym_type,&
                                              print_symmetry,&
                                              release_molsym
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: massunit
   USE string_utilities,                ONLY: uppercase
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters (in this module)

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_symmetry'

   PUBLIC :: write_symmetry

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Write symmetry information to output
!> \param particle_set  Atom coordinates and types
!> \param cell          Cell information
!> \param input_section Input
!> \par History
!> \author jgh
! **************************************************************************************************
   SUBROUTINE write_symmetry(particle_set, cell, input_section)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(section_vals_type), POINTER                   :: input_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'write_symmetry'

      CHARACTER(LEN=2), ALLOCATABLE, DIMENSION(:)        :: element
      CHARACTER(LEN=8)                                   :: csymm, esymm
      INTEGER                                            :: handle, i, iw, natom, plevel
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atype, z
      LOGICAL                                            :: check, molecular, pall, pcoor, pinertia, &
                                                            prmat, psymmele
      REAL(KIND=dp)                                      :: eps_geo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: coord, scoord
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(csym_type)                                    :: crys_sym
      TYPE(molsym_type), POINTER                         :: mol_sym
      TYPE(section_vals_type), POINTER                   :: section

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      NULLIFY (section)

      logger => cp_get_default_logger()
      iw = cp_print_key_unit_nr(logger=logger, &
                                basis_section=input_section, &
                                print_key_path="PRINT%SYMMETRY", &
                                extension=".symLog")

      IF (iw > 0) THEN
         section => section_vals_get_subs_vals(section_vals=input_section, &
                                               subsection_name="PRINT%SYMMETRY")
         CALL section_vals_val_get(section_vals=section, &
                                   keyword_name="MOLECULE", l_val=molecular)
         CALL section_vals_val_get(section_vals=section, &
                                   keyword_name="EPS_GEO", r_val=eps_geo)
         IF (molecular) THEN

            NULLIFY (mol_sym)

            natom = SIZE(particle_set)
            ALLOCATE (coord(3, natom), z(natom), weight(natom), atype(natom), element(natom))

            DO i = 1, natom
               CALL get_atomic_kind(particle_set(i)%atomic_kind, z=z(i))
               CALL get_atomic_kind(atomic_kind=particle_set(i)%atomic_kind, &
                                    kind_number=atype(i), element_symbol=element(i), mass=weight(i))
               coord(1:3, i) = particle_set(i)%r(1:3)
            END DO
            weight(:) = weight(:)/massunit

            CALL molecular_symmetry(mol_sym, eps_geo, coord, atype, weight)

            CALL section_vals_val_get(section_vals=section, &
                                      keyword_name="STANDARD_ORIENTATION", l_val=pcoor)
            CALL section_vals_val_get(section_vals=section, &
                                      keyword_name="INERTIA", l_val=pinertia)
            CALL section_vals_val_get(section_vals=section, &
                                      keyword_name="SYMMETRY_ELEMENTS", l_val=psymmele)
            CALL section_vals_val_get(section_vals=section, &
                                      keyword_name="ALL", l_val=pall)
            plevel = 0
            IF (pcoor) plevel = plevel + 1
            IF (pinertia) plevel = plevel + 10
            IF (psymmele) plevel = plevel + 100
            IF (pall) plevel = 1111111111

            CALL print_symmetry(mol_sym, coord, atype, element, z, weight, iw, plevel)

            CALL section_vals_val_get(section_vals=section, &
                                      keyword_name="CHECK_SYMMETRY", c_val=esymm)
            CALL uppercase(esymm)
            IF (TRIM(esymm) /= "NONE") THEN
               csymm = mol_sym%point_group_symbol
               CALL uppercase(csymm)
               check = TRIM(ADJUSTL(csymm)) == TRIM(ADJUSTL(esymm))
               IF (.NOT. check) THEN
                  CALL cp_warn(__LOCATION__, "Symmetry check failed: "// &
                               "Expected symmetry:"//TRIM(ADJUSTL(esymm))// &
                               "Calculated symmetry:"//TRIM(ADJUSTL(csymm)))
               END IF
               CPASSERT(check)
            END IF

            DEALLOCATE (coord, z, weight, atype, element)

            CALL release_molsym(mol_sym)

         ELSE
            ! Crystal symmetry

            natom = SIZE(particle_set)
            ALLOCATE (scoord(3, natom), atype(natom))

            DO i = 1, natom
               CALL get_atomic_kind(atomic_kind=particle_set(i)%atomic_kind, kind_number=atype(i))
               CALL real_to_scaled(scoord(1:3, i), particle_set(i)%r(1:3), cell)
            END DO

            CALL crys_sym_gen(crys_sym, scoord, atype, cell%hmat, delta=eps_geo, iounit=iw)

            CALL section_vals_val_get(section_vals=section, &
                                      keyword_name="ROTATION_MATRICES", l_val=prmat)
            CALL section_vals_val_get(section_vals=section, &
                                      keyword_name="ALL", l_val=pall)
            plevel = 0
            IF (prmat) plevel = plevel + 1
            IF (pall) plevel = 1111111111
            crys_sym%plevel = plevel

            CALL print_crys_symmetry(crys_sym)

            DEALLOCATE (scoord, atype)

            CALL release_csym_type(crys_sym)

         END IF

      END IF
      CALL cp_print_key_finished_output(iw, logger, input_section, "PRINT%SYMMETRY")

      CALL timestop(handle)

   END SUBROUTINE write_symmetry

! **************************************************************************************************

END MODULE cp_symmetry
