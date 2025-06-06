!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the kpoints section of the input
!> \par History
!>      init [07.2014]
!> \author JGH
! **************************************************************************************************
MODULE input_cp2k_kpoints
   USE bibliography,                    ONLY: MacDonald1978,&
                                              Monkhorst1976
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_create,&
                                              section_type
   USE input_val_types,                 ONLY: char_t,&
                                              real_t
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE string_utilities,                ONLY: newline,&
                                              s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_kpoints'

   INTEGER, PARAMETER                      :: use_real_wfn = 101, &
                                              use_complex_wfn = 100

   PUBLIC :: create_kpoints_section
   PUBLIC :: create_kpoint_set_section
   PUBLIC :: use_real_wfn, use_complex_wfn

CONTAINS

! **************************************************************************************************
!> \brief Creates the Kpoints section
!>  SECTION: &kpoint... &end
!>      SCHEME [None, Gamma, Monkhorst-Pack, MacDonald, General]
!>        { nx ny nz }
!>        { nx ny nz  sx sy sz }
!>      KPOINT           x1 y1 z1 w1
!>      SYMMETRY         [on, off]
!>      WAVEFUNCTION     [real, complex]
!>      FULL_GRID        [on, off]
!>      VERBOSE          [on, off]
!>      EPS_GEO          value
!>      PARALLEL_GROUP_SIZE   [-1,0,n]
!>
!> \param section the section to create
!> \author JGH
! **************************************************************************************************
   SUBROUTINE create_kpoints_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="KPOINTS", &
                          description="Sets up the kpoints.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SCHEME", &
                          description="Kpoint scheme to be used. Available options are:"//newline// &
                          "- `NONE`"//newline// &
                          "- `GAMMA`"//newline// &
                          "- `MONKHORST-PACK`"//newline// &
                          "- `MACDONALD`"//newline// &
                          "- `GENERAL`"//newline// &
                          newline// &
                          "For `MONKHORST-PACK` and `MACDONALD` the number of k points in all "// &
                          "3 dimensions has to be supplied along with the keyword. "// &
                          "E.g. `MONKHORST-PACK 12 12 8`", &
                          usage="SCHEME {KPMETHOD} {integer} {integer} ..", &
                          citations=(/Monkhorst1976, MacDonald1978/), &
                          n_var=-1, type_of_var=char_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="KPOINT", &
                          description="Specify kpoint coordinates and weight. ", &
                          usage="KPOINT x  y  z  w", repeats=.TRUE., &
                          n_var=4, type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="UNITS", &
                          description="Special k-points are defined either in units"// &
                          " of reciprocal lattice vectors or in Cartesian coordinates in uints of 2Pi/len."// &
                          " B_VECTOR: in multiples of the reciprocal lattice vectors (b)."// &
                          " CART_ANGSTROM: In units of 2*Pi/Angstrom."// &
                          " CART_BOHR: In units of 2*Pi/Bohr.", &
                          usage="UNITS <value>", type_of_var=char_t, default_c_val="B_VECTOR")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SYMMETRY", &
                          description="Use symmetry to reduce the number of kpoints.", &
                          usage="SYMMETRY <LOGICAL>", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FULL_GRID", &
                          description="Use full non-reduced kpoint grid.", &
                          usage="FULL_GRID <LOGICAL>", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VERBOSE", &
                          description="Verbose output information.", &
                          usage="VERBOSE <LOGICAL>", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_GEO", &
                          description="Accuracy in symmetry determination.", &
                          usage="EPS_GEO <real>", &
                          default_r_val=1.0e-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PARALLEL_GROUP_SIZE", &
                          description="Number of processors to be used for a single kpoint."// &
                          " This number must divide the total number of processes."// &
                          " The number of groups must divide the total number of kpoints."// &
                          " Value=-1 (smallest possible number of processes per group, satisfying the constraints)."// &
                          " Value=0 (all processes)."// &
                          " Value=n (exactly n processes).", &
                          usage="PARALLEL_GROUP_SIZE <integer>", &
                          default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WAVEFUNCTIONS", &
                          description="Use real/complex wavefunctions if possible.", &
                          usage="WAVEFUNCTIONS REAL", &
                          default_i_val=use_complex_wfn, &
                          enum_c_vals=s2a("REAL", "COMPLEX"), &
                          enum_desc=s2a("Use real wavefunctions (if possible by kpoints specified).", &
                                        "Use complex wavefunctions (default)."), &
                          enum_i_vals=(/use_real_wfn, use_complex_wfn/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_kpoints_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
!> \param section_name ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE create_kpoint_set_section(section, section_name)
      TYPE(section_type), POINTER                        :: section
      CHARACTER(LEN=*), OPTIONAL                         :: section_name

      CHARACTER(len=default_path_length)                 :: my_section_name
      TYPE(keyword_type), POINTER                        :: keyword

      IF (PRESENT(section_name)) THEN
         my_section_name = section_name
      ELSE
         my_section_name = "KPOINT_SET"
      END IF

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name=my_section_name, &
                          description="Specifies a k-point line to be calculated.", &
                          n_keywords=0, n_subsections=0, repeats=.TRUE.)
      ! keywords
      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SPECIAL_POINT", &
                          description="Name and coordinates of a special k-point", &
                          usage="SPECIAL_POINT GAMMA 0.0 0.0 0.0", n_var=-1, type_of_var=char_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      !
      CALL keyword_create(keyword, __LOCATION__, name="NPOINTS", &
                          description="Number of k-points along the line.", &
                          default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      !
      CALL keyword_create(keyword, __LOCATION__, name="UNITS", &
                          description="Special k-points are defined either in units"// &
                          " of reciprocal lattice vectors or in Cartesian coordinates in uints of 2Pi/len."// &
                          " B_VECTOR: in multiples of the reciprocal lattice vectors (b)."// &
                          " CART_ANGSTROM: In units of 2*Pi/Angstrom."// &
                          " CART_BOHR: In units of 2*Pi/Bohr.", &
                          usage="UNITS <value>", type_of_var=char_t, default_c_val="B_VECTOR")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_kpoint_set_section

END MODULE input_cp2k_kpoints
