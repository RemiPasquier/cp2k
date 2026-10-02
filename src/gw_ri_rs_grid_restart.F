!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Automatic restart I/O for optimized RI-RS grids.
!>
!> The versioned formatted file stores atom-relative points, positions, the four grid/metric
!> parameters, MAX_ITER of the writing run and AO/RI basis names. Only a matching file is reused;
!> any difference or corrupt record leaves the freshly initialized grid for cold optimization.
!> Basis definitions are not serialized: changes under the same name are not detected.
! **************************************************************************************************
MODULE gw_ri_rs_grid_restart
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE machine,                         ONLY: m_mov
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type,&
                                              rirs_grid_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: read_ri_rs_grid_restart, write_ri_rs_grid_restart
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_grid_restart'
   CHARACTER(LEN=*), PARAMETER, PRIVATE :: file_header = '# CP2K RI-RS grid restart v1'
   REAL(KIND=dp), PARAMETER, PRIVATE    :: position_tolerance = 1.0E-8_dp
   REAL(KIND=dp), PARAMETER, PRIVATE    :: relative_tolerance = 1.0E-12_dp
CONTAINS
! **************************************************************************************************
!> \brief Try to load the automatic grid restart on the I/O rank and broadcast the selected grid.
!> \param bs_env Bandstructure state containing the freshly initialized grid.
!> \param exact_reuse True only when all recorded metadata match and the file parses strictly.
! **************************************************************************************************
   SUBROUTINE read_ri_rs_grid_restart(bs_env, exact_reuse)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      LOGICAL, INTENT(OUT)                               :: exact_reuse

      CHARACTER(LEN=default_path_length)                 :: filename, reason
      INTEGER                                            :: handle, iatom, status
      LOGICAL                                            :: basis_diff, kind_diff, max_iter_diff, &
                                                            parameter_diff, position_diff
      TYPE(rirs_grid_type), ALLOCATABLE                  :: loaded_grids(:)

      CALL timeset('read_ri_rs_grid_restart', handle)
      filename = TRIM(bs_env%prefix)//'RI_RS_grid.dat'
      status = 0
      reason = ''
      IF (bs_env%para_env%is_source()) THEN
         CALL read_grid_on_source(bs_env, filename, loaded_grids, reason, basis_diff, kind_diff, &
                                  parameter_diff, position_diff, max_iter_diff)
         IF (LEN_TRIM(reason) == 0) THEN
            IF (position_diff) reason = 'positions, '
            IF (parameter_diff) reason = TRIM(reason)//' grid/metric parameters, '
            IF (kind_diff) reason = TRIM(reason)//' kind metadata, '
            IF (basis_diff) reason = TRIM(reason)//' basis names, '
            IF (max_iter_diff) reason = TRIM(reason)//' MAX_ITER, '
            IF (LEN_TRIM(reason) > 0) THEN
               reason = TRIM(ADJUSTL(reason(:LEN_TRIM(reason) - 1)))//' differ'
            ELSE
               status = 1
               DO iatom = 1, bs_env%n_atom
                  bs_env%ri_rs%atomic_grids(iatom)%raw_points(:, :) = loaded_grids(iatom)%raw_points
               END DO
            END IF
         END IF
      END IF
      CALL bs_env%para_env%bcast(status)
      CALL bs_env%para_env%bcast(reason)
      DO iatom = 1, bs_env%n_atom
         CALL bs_env%para_env%bcast(bs_env%ri_rs%atomic_grids(iatom)%raw_points)
      END DO
      exact_reuse = status == 1
      IF (bs_env%unit_nr > 0) THEN
         IF (exact_reuse) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I0,A)') 'RI-RS grid reused without optimization: ', &
               SUM(bs_env%ri_rs%atomic_grids(:)%npts), ' points'
         ELSE
            WRITE (bs_env%unit_nr, '(A,A)') 'WARNING: RI-RS grid restart not used: ', TRIM(reason)
         END IF
      END IF
      CALL timestop(handle)
   END SUBROUTINE read_ri_rs_grid_restart
! **************************************************************************************************
!> \brief Validate one grid restart on the I/O rank without changing the initialized grids.
!> \param bs_env Bandstructure state containing the freshly initialized grid.
!> \param filename Restart filename.
!> \param grids Temporary atom-relative grids, installed only on exact reuse.
!> \param reason Nonempty reason when parsing or hard identity fails.
!> \param basis_diff Whether AO or RI basis names differ.
!> \param kind_diff Whether the number or elements of kinds differ.
!> \param parameter_diff Whether a recorded numerical parameter differs.
!> \param position_diff Whether an atomic position differs by more than 1e-8 bohr.
!> \param max_iter_diff Whether MAX_ITER differs from the writing run.
! **************************************************************************************************
   SUBROUTINE read_grid_on_source(bs_env, filename, grids, reason, basis_diff, kind_diff, parameter_diff, &
                                  position_diff, max_iter_diff)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      TYPE(rirs_grid_type), ALLOCATABLE, INTENT(OUT)     :: grids(:)
      CHARACTER(LEN=*), INTENT(OUT)                      :: reason
      LOGICAL, INTENT(OUT)                               :: basis_diff, kind_diff, parameter_diff, &
                                                            position_diff, max_iter_diff

      CHARACTER(LEN=4*default_path_length)               :: line
      CHARACTER(LEN=default_path_length)                 :: fields(6)
      INTEGER                                            :: current_nkind, iatom, ikind, ios, ipt, &
                                                            jkind, kind_file, max_iter_file, &
                                                            natom_file, nkind_file, npts, u
      LOGICAL                                            :: valid
      REAL(KIND=dp)                                      :: parameters_file(4), point(3), &
                                                            position_file(3)

      reason = ''
      basis_diff = .FALSE.
      kind_diff = .FALSE.
      parameter_diff = .FALSE.
      position_diff = .FALSE.
      max_iter_diff = .FALSE.
      natom_file = 0
      nkind_file = 0
      ! The v1 ATOM records store the element, not the kind: if two kinds of one element have different
      ! basis names, a changed per-atom assignment cannot be detected, so such a file is never reused.
      current_nkind = SIZE(bs_env%ri_rs%atomic_kind_set)
      DO ikind = 1, current_nkind
         DO jkind = 1, ikind - 1
            IF (bs_env%ri_rs%atomic_kind_set(ikind)%element_symbol /= &
                bs_env%ri_rs%atomic_kind_set(jkind)%element_symbol) CYCLE
            IF (bs_env%basis_set_AO(ikind)%gto_basis_set%name /= &
                bs_env%basis_set_AO(jkind)%gto_basis_set%name .OR. &
                bs_env%basis_set_RI(ikind)%gto_basis_set%name /= &
                bs_env%basis_set_RI(jkind)%gto_basis_set%name) THEN
               reason = 'per-atom basis assignment is unavailable in v1'
               RETURN
            END IF
         END DO
      END DO
      OPEN (NEWUNIT=u, FILE=filename, STATUS='OLD', ACTION='READ', FORM='FORMATTED', IOSTAT=ios)
      IF (ios /= 0) THEN
         reason = 'file is unreadable'
         RETURN
      END IF
      CALL read_grid_line(u, line, ios)
      IF (ios /= 0 .OR. TRIM(line) /= file_header) reason = 'corrupt or missing format header'
      IF (LEN_TRIM(reason) == 0) THEN
         CALL read_grid_record(u, fields(:2), valid)
         IF (valid) CALL parse_grid_integer(fields(2), natom_file, valid)
         IF (.NOT. valid .OR. fields(1) /= 'NATOM') reason = 'corrupt NATOM record'
      END IF
      IF (LEN_TRIM(reason) == 0) THEN
         IF (natom_file /= bs_env%n_atom) reason = 'atom count differs'
      END IF
      IF (LEN_TRIM(reason) == 0) THEN
         CALL read_grid_record(u, fields(:5), valid)
         IF (valid) CALL parse_grid_reals(fields(2:5), parameters_file, valid)
         IF (.NOT. valid .OR. fields(1) /= 'PARAM') reason = 'corrupt PARAM record'
      END IF
      IF (LEN_TRIM(reason) == 0) THEN
         parameter_diff = .NOT. ALL(ABS(parameters_file - grid_parameters(bs_env)) <= &
                                    relative_tolerance*MAX(1.0_dp, ABS(grid_parameters(bs_env))))
         CALL read_grid_record(u, fields(:2), valid)
         IF (valid) CALL parse_grid_integer(fields(2), max_iter_file, valid)
         IF (.NOT. valid .OR. fields(1) /= 'MAX_ITER') THEN
            reason = 'corrupt MAX_ITER record'
         ELSE IF (max_iter_file < 1) THEN
            reason = 'corrupt MAX_ITER record'
         END IF
      END IF
      IF (LEN_TRIM(reason) == 0) THEN
         max_iter_diff = max_iter_file /= bs_env%ri_rs%grid_opt%max_iter
         CALL read_grid_record(u, fields(:2), valid)
         IF (valid) CALL parse_grid_integer(fields(2), nkind_file, valid)
         IF (.NOT. valid .OR. fields(1) /= 'NKIND' .OR. nkind_file < 1 .OR. nkind_file > bs_env%n_atom) THEN
            reason = 'corrupt NKIND record'
         END IF
      END IF
      IF (LEN_TRIM(reason) == 0) kind_diff = nkind_file /= current_nkind
      DO ikind = 1, nkind_file
         IF (LEN_TRIM(reason) > 0) EXIT
         CALL read_grid_record(u, fields(:5), valid)
         IF (valid) CALL parse_grid_integer(fields(2), kind_file, valid)
         IF (.NOT. valid .OR. fields(1) /= 'KIND') THEN
            reason = 'corrupt KIND record'
         ELSE IF (kind_file /= ikind) THEN
            reason = 'corrupt KIND record'
         ELSE IF (ikind <= current_nkind) THEN
            kind_diff = kind_diff .OR. fields(3) /= bs_env%ri_rs%atomic_kind_set(ikind)%element_symbol
            basis_diff = basis_diff .OR. fields(4) /= bs_env%basis_set_AO(ikind)%gto_basis_set%name .OR. &
                         fields(5) /= bs_env%basis_set_RI(ikind)%gto_basis_set%name
         END IF
      END DO
      IF (LEN_TRIM(reason) == 0) ALLOCATE (grids(bs_env%n_atom))
      DO iatom = 1, bs_env%n_atom
         IF (LEN_TRIM(reason) > 0) EXIT
         CALL read_grid_record(u, fields, valid)
         IF (valid) CALL parse_grid_integer(fields(2), npts, valid)
         IF (valid) CALL parse_grid_reals(fields(3:5), position_file, valid)
         IF (.NOT. valid .OR. fields(1) /= 'ATOM') THEN
            WRITE (reason, '(A,I0)') 'corrupt ATOM record ', iatom
         ELSE IF (fields(6) /= bs_env%ri_rs%atomic_kind_set( &
                  bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number)%element_symbol) THEN
            WRITE (reason, '(A,I0,A)') 'element of atom ', iatom, ' differs'
         ELSE IF (npts /= bs_env%ri_rs%atomic_grids(iatom)%npts) THEN
            WRITE (reason, '(A,I0,A)') 'point count of atom ', iatom, ' differs'
         ELSE
            position_diff = position_diff .OR. &
                            ANY(ABS(position_file - bs_env%ri_rs%particle_set(iatom)%r) > position_tolerance)
            grids(iatom)%npts = npts
            ALLOCATE (grids(iatom)%raw_points(3, npts))
            DO ipt = 1, npts
               CALL read_grid_record(u, fields(:3), valid)
               IF (valid) CALL parse_grid_reals(fields(:3), point, valid)
               IF (.NOT. valid) THEN
                  WRITE (reason, '(A,I0)') 'corrupt point record for atom ', iatom
                  EXIT
               END IF
               grids(iatom)%raw_points(:, ipt) = point
            END DO
         END IF
      END DO
      IF (LEN_TRIM(reason) == 0) THEN
         CALL read_grid_record(u, fields(:1), valid)
         IF (.NOT. valid .OR. fields(1) /= 'END') reason = 'corrupt or missing END record'
      END IF
      IF (LEN_TRIM(reason) == 0) THEN
         CALL read_grid_line(u, line, ios)
         IF (.NOT. IS_IOSTAT_END(ios)) reason = 'corrupt trailing data after END record'
      END IF
      CLOSE (u, IOSTAT=ios)
      IF (ios /= 0 .AND. LEN_TRIM(reason) == 0) reason = 'file close failed'
   END SUBROUTINE read_grid_on_source
! **************************************************************************************************
!> \brief Read exactly one complete formatted record, rejecting buffer overflow.
!> \param u File unit.
!> \param line Record buffer, larger than every record produced by the v1 writer.
!> \param ios Zero on success; EOF, read errors or overlong records otherwise.
! **************************************************************************************************
   SUBROUTINE read_grid_line(u, line, ios)
      INTEGER, INTENT(IN)                                :: u
      CHARACTER(LEN=*), INTENT(OUT)                      :: line
      INTEGER, INTENT(OUT)                               :: ios

      ! Nonadvancing input distinguishes a whole record from a silently truncated line.
      line = ''
      READ (u, '(A)', ADVANCE='NO', IOSTAT=ios) line
      IF (IS_IOSTAT_EOR(ios)) THEN
         ios = 0
      ELSE IF (ios == 0) THEN
         ios = 1
      END IF
   END SUBROUTINE read_grid_line
! **************************************************************************************************
!> \brief Split one record into exactly the required number of whitespace-separated fields.
!> \param u File unit.
!> \param fields Expected fields; neither null fields nor continuation records are accepted.
!> \param valid True when the whole record contains exactly SIZE(fields) bounded tokens.
! **************************************************************************************************
   SUBROUTINE read_grid_record(u, fields, valid)
      INTEGER, INTENT(IN)                                :: u
      CHARACTER(LEN=*), INTENT(OUT)                      :: fields(:)
      LOGICAL, INTENT(OUT)                               :: valid

      CHARACTER(LEN=4*default_path_length)               :: line
      INTEGER                                            :: first, i, ios, last, nfield

      fields = ''
      valid = .FALSE.
      CALL read_grid_line(u, line, ios)
      IF (ios /= 0) RETURN
      nfield = 0
      last = LEN_TRIM(line)
      i = 1
      DO WHILE (i <= last)
         IF (line(i:i) == ' ' .OR. line(i:i) == ACHAR(9)) THEN
            i = i + 1
            CYCLE
         END IF
         first = i
         DO WHILE (i <= last)
            IF (line(i:i) == ' ' .OR. line(i:i) == ACHAR(9)) EXIT
            i = i + 1
         END DO
         nfield = nfield + 1
         IF (nfield > SIZE(fields) .OR. i - first > LEN(fields)) RETURN
         fields(nfield) = line(first:i - 1)
         IF (fields(nfield) == '/' .OR. fields(nfield) == ',' .OR. &
             fields(nfield) == "''" .OR. fields(nfield) == '""') RETURN
      END DO
      valid = nfield == SIZE(fields)
   END SUBROUTINE read_grid_record
! **************************************************************************************************
!> \brief Convert a single integer token, excluding list-directed nulls and control syntax.
!> \param field Token.
!> \param value Parsed integer.
!> \param valid True when the token converts without an I/O error.
! **************************************************************************************************
   SUBROUTINE parse_grid_integer(field, value, valid)
      CHARACTER(LEN=*), INTENT(IN)                       :: field
      INTEGER, INTENT(OUT)                               :: value
      LOGICAL, INTENT(OUT)                               :: valid

      INTEGER                                            :: ios

      value = 0
      valid = .FALSE.
      IF (LEN_TRIM(field) == 0) RETURN
      IF (VERIFY(TRIM(field), '+-0123456789') /= 0) RETURN
      READ (field, *, IOSTAT=ios) value
      IF (ios /= 0) THEN
         value = 0
         RETURN
      END IF
      valid = .TRUE.
   END SUBROUTINE parse_grid_integer
! **************************************************************************************************
!> \brief Convert individual real tokens; reject nulls, repeat counts, comments and non-finite values.
!> \param fields Tokens.
!> \param values Parsed real values.
!> \param valid True only when every token converts to a finite real value.
! **************************************************************************************************
   SUBROUTINE parse_grid_reals(fields, values, valid)
      CHARACTER(LEN=*), INTENT(IN)                       :: fields(:)
      REAL(KIND=dp), INTENT(OUT)                         :: values(:)
      LOGICAL, INTENT(OUT)                               :: valid

      INTEGER                                            :: i, ios

      values = 0.0_dp
      valid = .FALSE.
      DO i = 1, SIZE(fields)
         IF (LEN_TRIM(fields(i)) == 0) RETURN
         IF (VERIFY(TRIM(fields(i)), '+-0123456789.eEdD') /= 0) RETURN
         READ (fields(i), *, IOSTAT=ios) values(i)
         IF (ios /= 0) RETURN
         IF (.NOT. ieee_is_finite(values(i))) RETURN
      END DO
      valid = .TRUE.
   END SUBROUTINE parse_grid_reals
! **************************************************************************************************
!> \brief Write the automatic formatted grid restart from the I/O rank.
!> \param bs_env Bandstructure state containing the final atom-relative grids.
! **************************************************************************************************
   SUBROUTINE write_ri_rs_grid_restart(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=default_path_length)                 :: filename
      CHARACTER(LEN=default_path_length+4)               :: temporary_filename
      INTEGER                                            :: handle, iatom, ikind, ios, ipt, u
      LOGICAL                                            :: target_exists, temporary_exists

      CALL timeset('write_ri_rs_grid_restart', handle)
      filename = TRIM(bs_env%prefix)//'RI_RS_grid.dat'
      temporary_filename = TRIM(filename)//'.tmp'
      IF (bs_env%para_env%is_source()) THEN
         OPEN (NEWUNIT=u, FILE=temporary_filename, STATUS='REPLACE', ACTION='WRITE', FORM='FORMATTED', IOSTAT=ios)
         CALL check_write(ios)
         WRITE (u, '(A)', IOSTAT=ios) file_header
         CALL check_write(ios)
         WRITE (u, '(A,1X,I0)', IOSTAT=ios) 'NATOM', bs_env%n_atom
         CALL check_write(ios)
         WRITE (u, '(A,4(1X,ES24.16))', IOSTAT=ios) 'PARAM', grid_parameters(bs_env)
         CALL check_write(ios)
         WRITE (u, '(A,1X,I0)', IOSTAT=ios) 'MAX_ITER', bs_env%ri_rs%grid_opt%max_iter
         CALL check_write(ios)
         WRITE (u, '(A,1X,I0)', IOSTAT=ios) 'NKIND', SIZE(bs_env%ri_rs%atomic_kind_set)
         CALL check_write(ios)
         DO ikind = 1, SIZE(bs_env%ri_rs%atomic_kind_set)
            WRITE (u, '(A,1X,I0,3(1X,A))', IOSTAT=ios) 'KIND', ikind, &
               TRIM(bs_env%ri_rs%atomic_kind_set(ikind)%element_symbol), &
               TRIM(bs_env%basis_set_AO(ikind)%gto_basis_set%name), &
               TRIM(bs_env%basis_set_RI(ikind)%gto_basis_set%name)
            CALL check_write(ios)
         END DO
         DO iatom = 1, bs_env%n_atom
            ikind = bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number
            WRITE (u, '(A,1X,I0,3(1X,ES24.16),1X,A)', IOSTAT=ios) 'ATOM', &
               bs_env%ri_rs%atomic_grids(iatom)%npts, bs_env%ri_rs%particle_set(iatom)%r, &
               TRIM(bs_env%ri_rs%atomic_kind_set(ikind)%element_symbol)
            CALL check_write(ios)
            DO ipt = 1, bs_env%ri_rs%atomic_grids(iatom)%npts
               WRITE (u, '(3(1X,ES24.16))', IOSTAT=ios) bs_env%ri_rs%atomic_grids(iatom)%raw_points(:, ipt)
               CALL check_write(ios)
            END DO
         END DO
         WRITE (u, '(A)', IOSTAT=ios) 'END'
         CALL check_write(ios)
         CLOSE (u, IOSTAT=ios)
         CALL check_write(ios)
         CALL m_mov(TRIM(temporary_filename), TRIM(filename))
         INQUIRE (FILE=TRIM(filename), EXIST=target_exists)
         INQUIRE (FILE=TRIM(temporary_filename), EXIST=temporary_exists)
         IF (.NOT. target_exists .OR. temporary_exists) THEN
            CPABORT('RI-RS grid restart: moving the temporary file into place failed')
         END IF
      END IF
      CALL bs_env%para_env%sync()
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(T2,A,I0,A)') 'RI-RS grid restart written: ', &
         SUM(bs_env%ri_rs%atomic_grids(:)%npts), ' points'
      CALL timestop(handle)
   END SUBROUTINE write_ri_rs_grid_restart
! **************************************************************************************************
!> \brief Return the four numerical parameters recorded with a grid.
!> \param bs_env Bandstructure state.
!> \return RS_AO_RATIO, cluster cutoff, RI metric cutoff and TIKHONOV.
! **************************************************************************************************
   FUNCTION grid_parameters(bs_env) RESULT(parameters)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp)                                      :: parameters(4)

      parameters = [bs_env%ri_rs%grid_opt%rs_ao_ratio, bs_env%ri_rs%grid_opt%cutoff_atomic_cluster, &
                    bs_env%ri_metric%cutoff_radius, bs_env%ri_rs%tikhonov]
   END FUNCTION grid_parameters
! **************************************************************************************************
!> \brief Abort after a failed grid-restart write operation.
!> \param ios I/O status.
! **************************************************************************************************
   SUBROUTINE check_write(ios)
      INTEGER, INTENT(IN)                                :: ios

      IF (ios /= 0) THEN
         CPABORT('RI-RS grid restart write failed')
      END IF
   END SUBROUTINE check_write
END MODULE gw_ri_rs_grid_restart
