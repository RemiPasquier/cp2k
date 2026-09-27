!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Writing and reading RI-RS grid restart files.
!>
!>        The file stores the atom-relative RI-RS grid of every atom together with the data the grid
!>        depends on. On reading, CP2K aborts unless these match the current run: atom count, order and
!>        elements; for every kind the ordered AO and RI basis definition (angular momenta, exponents and
!>        contraction coefficients); cell matrix and periodicity; RS_AO_RATIO, CUTOFF_ATOMIC_CLUSTER, the
!>        RI metric cutoff radius and TIKHONOV. Numerical basis, cell and parameter values are compared
!>        using abs(file - run) <= 1e-12 max(1, abs(run)); exact reuse additionally requires every position
!>        component to agree within 1e-8 bohr. The file records the writing invocation's mode, MAX_ITER
!>        and EPS_FILTER, which are printed, not checked; re-saving a read grid does not preserve its
!>        original optimization settings or achieved quality. Every record is read completely and must
!>        consist of exactly the expected fields; numeric fields are unquoted and are converted one value
!>        at a time (see get_record). Wrong field counts, invalid or non-finite numbers and trailing data
!>        abort the run. Names containing a double quote cannot be saved.
!> \par History
!>      09.2026 created
! **************************************************************************************************
MODULE gw_ri_rs_grid_restart
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type,&
                                              rirs_grid_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: write_ri_rs_grid_restart, read_ri_rs_grid_restart

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_grid_restart'
   CHARACTER(LEN=*), PARAMETER, PRIVATE :: file_header = "# CP2K RI-RS grid restart v1"
   INTEGER, PARAMETER, PRIVATE          :: max_fields = 10, max_line = 1024, &
                                           max_npts_atom = 10000000
   CHARACTER(LEN=*), PARAMETER, PRIVATE :: field_chars = &
                                           "0123456789+-._ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
   REAL(KIND=dp), PARAMETER, PRIVATE    :: rel_tol = 1.0E-12_dp, pos_tol = 1.0E-8_dp

CONTAINS

! **************************************************************************************************
!> \brief Write the atom-relative RI-RS grids (source rank only).
!> \param bs_env ...
!> \param filename ...
!> \param mode how the grid was obtained (recorded, not checked)
! **************************************************************************************************
   SUBROUTINE write_ri_rs_grid_restart(bs_env, filename, mode)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      CHARACTER(LEN=*), INTENT(IN)                       :: filename, mode

      CHARACTER(len=*), PARAMETER :: routineN = 'write_ri_rs_grid_restart'

      INTEGER                                            :: handle, i, iatom, ikind, ios, ipt, &
                                                            nkind, u
      REAL(KIND=dp), ALLOCATABLE                         :: ao(:), ri(:)

      CALL timeset(routineN, handle)
      CPASSERT(ALLOCATED(bs_env%ri_rs%atomic_grids))
      IF (bs_env%para_env%is_source()) THEN
         nkind = SIZE(bs_env%ri_rs%atomic_kind_set)
         ! Names are written quoted; check all of them before an existing file is replaced.
         DO ikind = 1, nkind
            IF (INDEX(bs_env%ri_rs%atomic_kind_set(ikind)%element_symbol, '"') > 0 .OR. &
                INDEX(bs_env%basis_set_AO(ikind)%gto_basis_set%name, '"') > 0 .OR. &
                INDEX(bs_env%basis_set_RI(ikind)%gto_basis_set%name, '"') > 0) THEN
               CPABORT("RI-RS grid restart: cannot save a name containing a double quote")
            END IF
         END DO
         CALL open_file(file_name=filename, unit_number=u, file_action="WRITE", file_status="REPLACE")
         WRITE (u, IOSTAT=ios, FMT="(A)") file_header
         CALL check_write(ios)
         WRITE (u, IOSTAT=ios, FMT="(A,1X,I0)") "NATOM", bs_env%n_atom
         CALL check_write(ios)
         WRITE (u, IOSTAT=ios, FMT="(A,1X,I0)") "NKIND", nkind
         CALL check_write(ios)
         WRITE (u, IOSTAT=ios, FMT="(A,3(1X,I0))") "PERIODIC", bs_env%ri_rs%cell%perd
         CALL check_write(ios)
         WRITE (u, IOSTAT=ios, FMT="(A,9(1X,ES24.16))") "CELL", bs_env%ri_rs%cell%hmat
         CALL check_write(ios)
         WRITE (u, IOSTAT=ios, FMT="(A,4(1X,ES24.16))") "PARAM", grid_params(bs_env)
         CALL check_write(ios)
         WRITE (u, IOSTAT=ios, FMT="(A,1X,A,1X,I0,1X,ES24.16)") "PROVENANCE", '"'//TRIM(mode)//'"', &
            bs_env%ri_rs%grid_opt%max_iter, bs_env%eps_filter
         CALL check_write(ios)
         DO ikind = 1, nkind
            CALL serialize_basis(bs_env%basis_set_AO(ikind)%gto_basis_set, ao)
            CALL serialize_basis(bs_env%basis_set_RI(ikind)%gto_basis_set, ri)
            WRITE (u, IOSTAT=ios, FMT="(A,1X,I0,3(1X,A),2(1X,I0))") "KIND", ikind, &
               '"'//TRIM(bs_env%ri_rs%atomic_kind_set(ikind)%element_symbol)//'"', &
               '"'//TRIM(bs_env%basis_set_AO(ikind)%gto_basis_set%name)//'"', &
               '"'//TRIM(bs_env%basis_set_RI(ikind)%gto_basis_set%name)//'"', SIZE(ao), SIZE(ri)
            CALL check_write(ios)
            DO i = 1, SIZE(ao)
               WRITE (u, IOSTAT=ios, FMT="(ES24.16)") ao(i)
               CALL check_write(ios)
            END DO
            DO i = 1, SIZE(ri)
               WRITE (u, IOSTAT=ios, FMT="(ES24.16)") ri(i)
               CALL check_write(ios)
            END DO
            DEALLOCATE (ao, ri)
         END DO
         DO iatom = 1, bs_env%n_atom
            ikind = bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number
            WRITE (u, IOSTAT=ios, FMT="(A,2(1X,I0),3(1X,ES24.16),1X,I0)") "ATOM", iatom, ikind, &
               bs_env%ri_rs%particle_set(iatom)%r(:), bs_env%ri_rs%atomic_grids(iatom)%npts
            CALL check_write(ios)
            DO ipt = 1, bs_env%ri_rs%atomic_grids(iatom)%npts
               WRITE (u, IOSTAT=ios, FMT="(3(1X,ES24.16))") bs_env%ri_rs%atomic_grids(iatom)%raw_points(:, ipt)
               CALL check_write(ios)
            END DO
         END DO
         CALL close_file(unit_number=u)
      END IF
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, "(T2,A,I0,A,A)") "RI-RS grid points written: ", &
         SUM(bs_env%ri_rs%atomic_grids(:)%npts), " to ", filename
      CALL timestop(handle)

   END SUBROUTINE write_ri_rs_grid_restart

! **************************************************************************************************
!> \brief Read an RI-RS grid restart file and validate it against the current run.
!> \param bs_env ...
!> \param filename ...
!> \param exact .TRUE.: atom positions must match (the grid is used unchanged);
!>              .FALSE.: positions may differ (warm start of the grid optimization)
!> \param grids atom-relative grids read from the file
! **************************************************************************************************
   SUBROUTINE read_ri_rs_grid_restart(bs_env, filename, exact, grids)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      LOGICAL, INTENT(IN)                                :: exact
      TYPE(rirs_grid_type), ALLOCATABLE, INTENT(OUT)     :: grids(:)

      CHARACTER(len=*), PARAMETER :: routineN = 'read_ri_rs_grid_restart'

      CHARACTER(LEN=default_string_length)               :: ao_name, element, mode, ri_name
      CHARACTER(LEN=max_line)                            :: line
      INTEGER                                            :: handle, i, iatom, ikind, ios, ipt, &
                                                            max_iter_file, n, nao_file, nkind, &
                                                            nri_file, perd_file(3), u
      INTEGER(KIND=int_8)                                :: npts_total
      INTEGER, DIMENSION(2, max_fields)                  :: fb
      REAL(KIND=dp)                                      :: cell_file(9), eps_file, max_shift, &
                                                            params_file(4), r_file(3), x(3)
      REAL(KIND=dp), ALLOCATABLE                         :: ao(:), ri(:), values(:)

      CALL timeset(routineN, handle)
      nkind = SIZE(bs_env%ri_rs%atomic_kind_set)
      ALLOCATE (grids(bs_env%n_atom))
      CALL open_file(file_name=filename, unit_number=u, file_action="READ", file_status="OLD")

      CALL get_record(u, line, 0, fb)
      IF (TRIM(line) /= file_header) THEN
         CPABORT("Unsupported RI-RS grid restart file")
      END IF

      CALL get_record(u, line, 2, fb, "NATOM")
      IF (field_int(line, fb(:, 2), "NATOM") /= bs_env%n_atom) THEN
         CPABORT("RI-RS grid restart: atom count differs")
      END IF
      CALL get_record(u, line, 2, fb, "NKIND")
      IF (field_int(line, fb(:, 2), "NKIND") /= nkind) THEN
         CPABORT("RI-RS grid restart: kind count differs")
      END IF
      CALL get_record(u, line, 4, fb, "PERIODIC")
      DO i = 1, 3
         perd_file(i) = field_int(line, fb(:, i + 1), "PERIODIC")
      END DO
      IF (ANY(perd_file /= bs_env%ri_rs%cell%perd)) THEN
         CPABORT("RI-RS grid restart: periodicity differs")
      END IF
      CALL get_record(u, line, 10, fb, "CELL")
      DO i = 1, 9
         cell_file(i) = field_real(line, fb(:, i + 1), "CELL")
      END DO
      IF (.NOT. matches(cell_file, RESHAPE(bs_env%ri_rs%cell%hmat, [9]))) THEN
         CPABORT("RI-RS grid restart: cell differs")
      END IF
      CALL get_record(u, line, 5, fb, "PARAM")
      DO i = 1, 4
         params_file(i) = field_real(line, fb(:, i + 1), "PARAM")
      END DO
      IF (.NOT. matches(params_file, grid_params(bs_env))) THEN
         CPABORT("RI-RS grid restart: grid or RI metric parameters differ")
      END IF
      CALL get_record(u, line, 4, fb, "PROVENANCE")
      CALL field_str(line, fb(:, 2), "PROVENANCE", mode)
      max_iter_file = field_int(line, fb(:, 3), "PROVENANCE")
      eps_file = field_real(line, fb(:, 4), "PROVENANCE")
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, "(T2,A,A,A,I0,A,ES10.3)") &
         "RI-RS grid restart file written by: ", TRIM(mode), ", MAX_ITER ", max_iter_file, &
         ", EPS_FILTER ", eps_file

      DO ikind = 1, nkind
         CALL get_record(u, line, 7, fb, "KIND")
         IF (field_int(line, fb(:, 2), "KIND") /= ikind) THEN
            CPABORT("RI-RS grid restart: bad KIND record")
         END IF
         CALL field_str(line, fb(:, 3), "KIND", element)
         CALL field_str(line, fb(:, 4), "KIND", ao_name)
         CALL field_str(line, fb(:, 5), "KIND", ri_name)
         nao_file = field_int(line, fb(:, 6), "KIND")
         nri_file = field_int(line, fb(:, 7), "KIND")
         IF (TRIM(element) /= TRIM(bs_env%ri_rs%atomic_kind_set(ikind)%element_symbol) .OR. &
             TRIM(ao_name) /= TRIM(bs_env%basis_set_AO(ikind)%gto_basis_set%name) .OR. &
             TRIM(ri_name) /= TRIM(bs_env%basis_set_RI(ikind)%gto_basis_set%name)) THEN
            CPABORT("RI-RS grid restart: element or basis set name differs")
         END IF
         CALL serialize_basis(bs_env%basis_set_AO(ikind)%gto_basis_set, ao)
         CALL serialize_basis(bs_env%basis_set_RI(ikind)%gto_basis_set, ri)
         IF (nao_file /= SIZE(ao) .OR. nri_file /= SIZE(ri)) THEN
            CPABORT("RI-RS grid restart: basis set definition differs")
         END IF
         ALLOCATE (values(nao_file + nri_file))
         DO i = 1, nao_file + nri_file
            CALL get_record(u, line, 1, fb)
            values(i) = field_real(line, fb(:, 1), "basis value")
         END DO
         IF (.NOT. matches(values, [ao, ri])) THEN
            CPABORT("RI-RS grid restart: basis set definition differs")
         END IF
         DEALLOCATE (ao, ri, values)
      END DO

      npts_total = 0
      max_shift = 0.0_dp
      DO iatom = 1, bs_env%n_atom
         CALL get_record(u, line, 7, fb, "ATOM")
         i = field_int(line, fb(:, 2), "ATOM")
         n = field_int(line, fb(:, 3), "ATOM")
         IF (i /= iatom .OR. n /= bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number) THEN
            CPABORT("RI-RS grid restart: bad ATOM record")
         END IF
         DO i = 1, 3
            r_file(i) = field_real(line, fb(:, i + 3), "ATOM")
         END DO
         n = field_int(line, fb(:, 7), "ATOM")
         max_shift = MAX(max_shift, MAXVAL(ABS(r_file - bs_env%ri_rs%particle_set(iatom)%r(:))))
         IF (exact .AND. max_shift > pos_tol) THEN
            CPABORT("RI-RS grid restart: atom positions differ")
         END IF
         IF (n < 1 .OR. n > max_npts_atom) THEN
            CPABORT("RI-RS grid restart: point count out of range")
         END IF
         npts_total = npts_total + n
         IF (npts_total >= HUGE(1)) THEN
            CPABORT("RI-RS grid restart: too many grid points")
         END IF
         grids(iatom)%npts = n
         ALLOCATE (grids(iatom)%raw_points(3, n))
         DO ipt = 1, n
            CALL get_record(u, line, 3, fb)
            DO i = 1, 3
               x(i) = field_real(line, fb(:, i), "grid point")
            END DO
            grids(iatom)%raw_points(:, ipt) = x
         END DO
      END DO

      ! Require EOF after the last grid-point record; reject any additional record.
      READ (u, "(A)", ADVANCE="NO", SIZE=n, IOSTAT=ios) line
      IF (.NOT. IS_IOSTAT_END(ios)) THEN
         CPABORT("RI-RS grid restart: data after the last atom")
      END IF
      CALL close_file(unit_number=u)
      ! All ranks read the file; wait for every reader before a SAVE_GRID_FILE_NAME with the same name
      ! may replace it.
      CALL bs_env%para_env%sync()

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, "(T2,A,A,A,ES10.3,A)") &
         "RI-RS grid read from: ", filename, " (largest atom displacement ", max_shift, " bohr)"
      CALL timestop(handle)

   END SUBROUTINE read_ri_rs_grid_restart

! **************************************************************************************************
!> \brief Read one complete restart record with exactly nfields explicit values.
!>
!>        Fields are separated by blanks or tabs. A field is either a double-quoted character field
!>        (possibly empty, no quote inside), or a nonempty unquoted run of letters, digits and the
!>        characters + - . _ ; every other character aborts. field_int/field_real reject quoted fields,
!>        so list-directed separators, slashes, repeat counts and null values never reach a numeric
!>        conversion and each numeric field supplies exactly one value.
!> \param u unit
!> \param line the record
!> \param nfields expected number of fields (0: the caller compares the whole record)
!> \param fb first and last character of every field in line (quotes excluded)
!> \param key expected first field (optional)
! **************************************************************************************************
   SUBROUTINE get_record(u, line, nfields, fb, key)
      INTEGER, INTENT(IN)                                :: u
      CHARACTER(LEN=*), INTENT(OUT)                      :: line
      INTEGER, INTENT(IN)                                :: nfields
      INTEGER, DIMENSION(2, max_fields), INTENT(OUT)     :: fb
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: key

      INTEGER                                            :: ifield, ios, nread, pos, q

      line = ""
      fb = 0
      READ (u, "(A)", ADVANCE="NO", SIZE=nread, IOSTAT=ios) line
      ! Require EOR; ios == 0 means the buffer filled before EOR was observed.
      IF (.NOT. IS_IOSTAT_EOR(ios)) THEN
         CPABORT("RI-RS grid restart: truncated, overlong or unreadable record")
      END IF
      IF (nread >= LEN(line)) THEN
         CPABORT("RI-RS grid restart: overlong record")
      END IF
      IF (nfields == 0) RETURN
      CPASSERT(nfields <= max_fields)
      ifield = 0
      pos = 1
      DO WHILE (pos <= nread)
         IF (is_blank(line(pos:pos))) THEN
            pos = pos + 1
            CYCLE
         END IF
         ifield = ifield + 1
         IF (ifield > nfields) THEN
            CPABORT("RI-RS grid restart: wrong number of fields in a record")
         END IF
         IF (line(pos:pos) == '"') THEN
            ! q: position of the closing quote relative to the opening one; q == 0: no closing quote,
            ! q == 1: empty character field
            q = INDEX(line(pos + 1:nread), '"')
            IF (q == 0) THEN
               CPABORT("RI-RS grid restart: bad quoted field")
            END IF
            fb(:, ifield) = [pos + 1, pos + q - 1]
            pos = pos + q + 1
            IF (pos <= nread) THEN
               IF (.NOT. is_blank(line(pos:pos))) THEN
                  CPABORT("RI-RS grid restart: bad quoted field")
               END IF
            END IF
         ELSE
            fb(1, ifield) = pos
            DO WHILE (pos <= nread)
               IF (is_blank(line(pos:pos))) EXIT
               IF (VERIFY(line(pos:pos), field_chars) /= 0) THEN
                  CPABORT("RI-RS grid restart: invalid character in a record")
               END IF
               pos = pos + 1
            END DO
            fb(2, ifield) = pos - 1
         END IF
      END DO
      IF (ifield /= nfields) THEN
         CPABORT("RI-RS grid restart: wrong number of fields in a record")
      END IF
      IF (PRESENT(key)) THEN
         IF (line(fb(1, 1):fb(2, 1)) /= key) THEN
            CPABORT("RI-RS grid restart: unexpected record "//line(fb(1, 1):fb(2, 1)))
         END IF
      END IF

   END SUBROUTINE get_record

! **************************************************************************************************
!> \brief Blank or tab.
!> \param c ...
!> \return ...
! **************************************************************************************************
   LOGICAL FUNCTION is_blank(c)
      CHARACTER, INTENT(IN)                              :: c

      is_blank = (c == " " .OR. c == CHAR(9))

   END FUNCTION is_blank

! **************************************************************************************************
!> \brief Convert one unquoted field of a record from get_record to an integer; abort otherwise.
!> \param line the record
!> \param b first and last character of the field
!> \param what record name for the error message
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION field_int(line, b, what)
      CHARACTER(LEN=*), INTENT(IN)                       :: line
      INTEGER, DIMENSION(2), INTENT(IN)                  :: b
      CHARACTER(LEN=*), INTENT(IN)                       :: what

      INTEGER                                            :: ios

      IF (b(1) > 1) THEN
         IF (line(b(1) - 1:b(1) - 1) == '"') THEN
            CPABORT("RI-RS grid restart: quoted numeric field in a "//what//" record")
         END IF
      END IF
      READ (line(b(1):b(2)), *, IOSTAT=ios) field_int
      IF (ios /= 0) THEN
         CPABORT("RI-RS grid restart: bad "//what//" record")
      END IF

   END FUNCTION field_int

! **************************************************************************************************
!> \brief Convert one unquoted field of a record from get_record to a finite real; abort otherwise.
!> \param line the record
!> \param b first and last character of the field
!> \param what record name for the error message
!> \return ...
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION field_real(line, b, what)
      CHARACTER(LEN=*), INTENT(IN)                       :: line
      INTEGER, DIMENSION(2), INTENT(IN)                  :: b
      CHARACTER(LEN=*), INTENT(IN)                       :: what

      INTEGER                                            :: ios

      IF (b(1) > 1) THEN
         IF (line(b(1) - 1:b(1) - 1) == '"') THEN
            CPABORT("RI-RS grid restart: quoted numeric field in a "//what//" record")
         END IF
      END IF
      READ (line(b(1):b(2)), *, IOSTAT=ios) field_real
      IF (ios /= 0) THEN
         CPABORT("RI-RS grid restart: bad "//what//" record")
      END IF
      IF (.NOT. ieee_is_finite(field_real)) THEN
         CPABORT("RI-RS grid restart: non-finite value in a "//what//" record")
      END IF

   END FUNCTION field_real

! **************************************************************************************************
!> \brief Copy one field of a record from get_record to a string; abort if it does not fit.
!> \param line the record
!> \param b first and last character of the field
!> \param what record name for the error message
!> \param s the field
! **************************************************************************************************
   SUBROUTINE field_str(line, b, what, s)
      CHARACTER(LEN=*), INTENT(IN)                       :: line
      INTEGER, DIMENSION(2), INTENT(IN)                  :: b
      CHARACTER(LEN=*), INTENT(IN)                       :: what
      CHARACTER(LEN=*), INTENT(OUT)                      :: s

      IF (b(2) - b(1) + 1 > LEN(s)) THEN
         CPABORT("RI-RS grid restart: overlong field in a "//what//" record")
      END IF
      s = line(b(1):b(2))

   END SUBROUTINE field_str

! **************************************************************************************************
!> \brief Abort through CP2K if a restart-file WRITE failed.
!> \param ios IOSTAT of the WRITE
! **************************************************************************************************
   SUBROUTINE check_write(ios)
      INTEGER, INTENT(IN)                                :: ios

      IF (ios /= 0) THEN
         CPABORT("RI-RS grid restart: write failed")
      END IF

   END SUBROUTINE check_write

! **************************************************************************************************
!> \brief Complete, ordered numerical definition of a Gaussian basis set.
!> \param basis ...
!> \param values nset, then per set: npgf, lmin, lmax, nshell, exponents, and per shell l and the
!>               contraction coefficients
! **************************************************************************************************
   SUBROUTINE serialize_basis(basis, values)
      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: values(:)

      INTEGER                                            :: iset, ishell, n, npgf

      n = 1
      DO iset = 1, basis%nset
         n = n + 4 + basis%npgf(iset) + basis%nshell(iset)*(1 + basis%npgf(iset))
      END DO
      ALLOCATE (values(n))
      n = 1
      values(n) = REAL(basis%nset, dp)
      DO iset = 1, basis%nset
         npgf = basis%npgf(iset)
         values(n + 1:n + 4) = REAL([npgf, basis%lmin(iset), basis%lmax(iset), basis%nshell(iset)], dp)
         values(n + 5:n + 4 + npgf) = basis%zet(1:npgf, iset)
         n = n + 4 + npgf
         DO ishell = 1, basis%nshell(iset)
            values(n + 1) = REAL(basis%l(ishell, iset), dp)
            values(n + 2:n + 1 + npgf) = basis%gcc(1:npgf, ishell, iset)
            n = n + 1 + npgf
         END DO
      END DO

   END SUBROUTINE serialize_basis

! **************************************************************************************************
!> \brief Parameters that determine the optimized grid.
!> \param bs_env ...
!> \return RS_AO_RATIO, CUTOFF_ATOMIC_CLUSTER, RI metric cutoff radius, TIKHONOV
! **************************************************************************************************
   FUNCTION grid_params(bs_env) RESULT(p)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp)                                      :: p(4)

      p = [bs_env%ri_rs%grid_opt%rs_ao_ratio, bs_env%ri_rs%grid_opt%cutoff_atomic_cluster, &
           bs_env%ri_metric%cutoff_radius, bs_env%ri_rs%tikhonov]

   END FUNCTION grid_params

! **************************************************************************************************
!> \brief All values finite and |a - b| <= rel_tol max(1, |b|).
!> \param a ...
!> \param b ...
!> \return ...
! **************************************************************************************************
   LOGICAL FUNCTION matches(a, b)
      REAL(KIND=dp), INTENT(IN)                          :: a(:), b(:)

      matches = .FALSE.
      IF (SIZE(a) /= SIZE(b)) RETURN
      IF (.NOT. ALL(ieee_is_finite(a))) RETURN
      IF (.NOT. ALL(ieee_is_finite(b))) RETURN
      matches = ALL(ABS(a - b) <= rel_tol*MAX(1.0_dp, ABS(b)))

   END FUNCTION matches

END MODULE gw_ri_rs_grid_restart
