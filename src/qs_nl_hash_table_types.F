!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief A simple hash table of integer keys, using hash function:
!>          H(k) = (k*p) mod n + 1
!>        where:
!>          k = key
!>          p = a prime number >= n
!>          n = size of the hash table
!>         And collision resolvation is done by open addressing with linear
!>         probing.
!>
!>         The table consists of an array of (key,val) pairs, and
!>         there are no intermediate buckets. For every new entry (k,v):
!>         We first look up slot H(k), and if it already contains an entry,
!>         then move to the next empty slot using a predefined linear probing
!>         sequence (e.g. iterate from slots H(k) to n, and then 1 to H(k)-1).
!>         When we look up, we use the same probing sequence.
!>
!>         Derived from  qs_fb_hash_table_types.F  (Mark Tucker, Jun 2016)
! **************************************************************************************************
MODULE qs_nl_hash_table_types

   USE kinds,                           ONLY: int_8
   USE qs_hash_table_functions,         ONLY: hash_table_matching_prime
   USE qs_neighbor_list_types,          ONLY: neighbor_list_task_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! public types
   PUBLIC :: nl_hash_table_obj

! public methods
   PUBLIC :: nl_hash_table_create, & !create new table
             nl_hash_table_release, & !destroy existing table
             nl_hash_table_add, & !add a new entry to the table
             nl_hash_table_get_from_index, & !return the value from the specified index of the table
             nl_hash_table_is_null, &
             nl_hash_table_status

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_nl_hash_table_types'

! key value indicating an empty slot
   INTEGER(KIND=int_8), PARAMETER, PRIVATE :: EMPTY_KEY = -1_int_8
! Parameters related to automatic resizing of the hash_table:
! Resize by EXPAND_FACTOR if total no. slots / no. of filled slots < ENLARGE_RATIO
   INTEGER, PARAMETER, PRIVATE :: ENLARGE_RATIO = 1
   INTEGER, PARAMETER, PRIVATE :: REDUCE_RATIO = 3
   INTEGER, PARAMETER, PRIVATE :: EXPAND_FACTOR = 2
   INTEGER, PARAMETER, PRIVATE :: SHRINK_FACTOR = 2

! **************************************************************************************************
!> \brief hash table entry data type
!> \param key       : key of the entry
!> \param val       : value of the entry
! **************************************************************************************************
   TYPE nl_hash_table_element
      INTEGER(KIND=int_8) :: key = -1_int_8
      TYPE(neighbor_list_task_type), POINTER :: val => NULL()
   END TYPE nl_hash_table_element

! **************************************************************************************************
!> \brief data defining a hash table using open addressing for collision
!>        resolvation. Uses simple entry structure to be memory efficient
!>        as well as small overhead
!> \param table     : hash table data area
!> \param nelements : number of non-empty slots in table
!> \param nmax      : max number of slots in table
!> \param prime     : prime number used in the hash function
! **************************************************************************************************
   TYPE nl_hash_table_data
      TYPE(nl_hash_table_element), DIMENSION(:), POINTER :: table => NULL()
      INTEGER :: nelements = -1
      INTEGER :: nmax = -1
      INTEGER :: prime = -1
   END TYPE nl_hash_table_data

! **************************************************************************************************
!> \brief the object container which allows for the creation of an array
!>        of pointers to nl_hash_table objects
!> \param obj : pointer to the nl_hash_table object
! **************************************************************************************************
   TYPE nl_hash_table_obj
      TYPE(nl_hash_table_data), POINTER, PRIVATE :: obj => NULL()
   END TYPE nl_hash_table_obj

CONTAINS

! **************************************************************************************************
!> \brief Add element to a hash table, auto resize if necessary
!> \param hash_table : the nl_hash_table object
!> \param key        : key of the element
!> \param val        : value of the element
! **************************************************************************************************
   RECURSIVE SUBROUTINE nl_hash_table_add(hash_table, key, val)
      TYPE(nl_hash_table_obj), INTENT(INOUT)             :: hash_table
      INTEGER(KIND=int_8), INTENT(IN)                    :: key
      TYPE(neighbor_list_task_type), INTENT(IN), POINTER :: val

      INTEGER                                            :: islot
      LOGICAL                                            :: check_ok

      check_ok = nl_hash_table_has_data(hash_table)
      CPASSERT(check_ok)

      ! check hash table size, if too small rehash in a larger table
      IF (hash_table%obj%nelements*ENLARGE_RATIO .GE. hash_table%obj%nmax) THEN
         CALL nl_hash_table_rehash(hash_table=hash_table, nmax=hash_table%obj%nmax*EXPAND_FACTOR)
      END IF

      ! find the right slot for the given key
      islot = nl_hash_table_linear_probe(hash_table, key)
      CPASSERT(islot > 0)

      ! add a new task to the list of tasks with that key
      IF (hash_table%obj%table(islot)%key == EMPTY_KEY) THEN
         hash_table%obj%nelements = hash_table%obj%nelements + 1
         hash_table%obj%table(islot)%key = key
      END IF

      ! If a task exists, we make our new task point to that i.e. adding it to the beginning of the list
      IF (ASSOCIATED(hash_table%obj%table(islot)%val)) THEN
         val%next => hash_table%obj%table(islot)%val
      END IF

      ! store the (maybe new) first item in the list in the hash table
      hash_table%obj%table(islot)%val => val
   END SUBROUTINE nl_hash_table_add

! **************************************************************************************************
!> \brief Creates and initialises an empty nl_hash_table object
!> \param hash_table : the nl_hash_table object, its content must be NULL and cannot be UNDEFINED
!> \param nmax       : total size of the table, optional. If absent default size is 1.
! **************************************************************************************************
   SUBROUTINE nl_hash_table_create(hash_table, nmax)
      TYPE(nl_hash_table_obj), INTENT(INOUT)             :: hash_table
      INTEGER, INTENT(IN), OPTIONAL                      :: nmax

      INTEGER                                            :: my_nmax
      LOGICAL                                            :: check_ok

      check_ok = .NOT. nl_hash_table_has_data(hash_table)
      CPASSERT(check_ok)
      ALLOCATE (hash_table%obj)
      NULLIFY (hash_table%obj%table)
      hash_table%obj%nmax = 0
      hash_table%obj%nelements = 0
      hash_table%obj%prime = 2
      my_nmax = 1
      IF (PRESENT(nmax)) my_nmax = nmax
      CALL nl_hash_table_init(hash_table=hash_table, nmax=my_nmax)

   END SUBROUTINE nl_hash_table_create

! **************************************************************************************************
!> \brief Retrieve value from a hash table given a specified index
!> \param hash_table : the nl_hash_table object
!> \param idx        : the index to retrieve the data for
!> \param val        : output value, might be unassociated if there is no data with that index
! **************************************************************************************************
   SUBROUTINE nl_hash_table_get_from_index(hash_table, idx, val)
      TYPE(nl_hash_table_obj), INTENT(IN)                :: hash_table
      INTEGER, INTENT(IN)                                :: idx
      TYPE(neighbor_list_task_type), INTENT(OUT), &
         POINTER                                         :: val

      LOGICAL                                            :: check_ok

      CPASSERT((idx .GT. 0) .AND. (idx .LE. hash_table%obj%nmax))

      check_ok = nl_hash_table_has_data(hash_table)
      CPASSERT(check_ok)

      val => hash_table%obj%table(idx)%val

   END SUBROUTINE nl_hash_table_get_from_index

! **************************************************************************************************
!> \brief check if the object has data associated to it
!> \param hash_table : the nl_hash_table object in question
!> \return : true if hash_table%obj is associated, false otherwise
! **************************************************************************************************
   PURE FUNCTION nl_hash_table_has_data(hash_table) RESULT(res)
      TYPE(nl_hash_table_obj), INTENT(IN)                :: hash_table
      LOGICAL                                            :: res

      res = ASSOCIATED(hash_table%obj)
   END FUNCTION nl_hash_table_has_data

! **************************************************************************************************
!> \brief Initialises a nl_hash_table object
!> \param hash_table : the nl_hash_table object, its content must be NULL and cannot be UNDEFINED
!> \param nmax       : new size of the table, optional. If absent use the old size
! **************************************************************************************************
   SUBROUTINE nl_hash_table_init(hash_table, nmax)
      TYPE(nl_hash_table_obj), INTENT(INOUT)             :: hash_table
      INTEGER, INTENT(IN), OPTIONAL                      :: nmax

      INTEGER                                            :: ii, my_nmax, two_to_power
      LOGICAL                                            :: check_ok

      check_ok = nl_hash_table_has_data(hash_table)
      CPASSERT(check_ok)
      my_nmax = hash_table%obj%nmax
      IF (PRESENT(nmax)) my_nmax = nmax

      ! table length should always be power of 2. Find the least
      ! power that is greater or equal to my_nmax
      two_to_power = 1 ! = 2**0
      DO WHILE (two_to_power .LT. my_nmax)
         two_to_power = 2*two_to_power
      END DO
      my_nmax = two_to_power

      IF (ASSOCIATED(hash_table%obj%table)) THEN
         IF (SIZE(hash_table%obj%table) .NE. my_nmax) THEN
            DEALLOCATE (hash_table%obj%table)
            ALLOCATE (hash_table%obj%table(my_nmax))
         END IF
      ELSE
         ALLOCATE (hash_table%obj%table(my_nmax))
      END IF
      hash_table%obj%nmax = my_nmax
      hash_table%obj%prime = hash_table_matching_prime(my_nmax)

      ! initiate element to be "empty"
      DO ii = 1, hash_table%obj%nmax
         hash_table%obj%table(ii)%key = EMPTY_KEY
         NULLIFY (hash_table%obj%table(ii)%val)
      END DO
      hash_table%obj%nelements = 0
   END SUBROUTINE nl_hash_table_init

! **************************************************************************************************
!> \brief Initialises a nl_hash_table object
!> \param hash_table : the nl_hash_table object, its content must be NULL and cannot be UNDEFINED
!> \param key ...
!> \param is_null ...
! **************************************************************************************************
   SUBROUTINE nl_hash_table_is_null(hash_table, key, is_null)
      TYPE(nl_hash_table_obj), INTENT(IN)                :: hash_table
      INTEGER, INTENT(IN)                                :: key
      LOGICAL, INTENT(OUT)                               :: is_null

      LOGICAL                                            :: check_ok

      check_ok = nl_hash_table_has_data(hash_table)
      CPASSERT(check_ok)
      check_ok = (key .LE. hash_table%obj%nmax)
      CPASSERT(check_ok)

      is_null = .FALSE.
      IF (EMPTY_KEY == hash_table%obj%table(key)%key) THEN !.OR.
         !NULLIFY(hash_table%obj%table(key)%val)
         is_null = .TRUE.
      END IF
   END SUBROUTINE nl_hash_table_is_null

! **************************************************************************************************
!> \brief Rehash table. If nmax is present, then also change the table size
!>        to MAX(nmax, number_of_non_empty_elements).
!> \param hash_table      : the nl_hash_table object
!> \param nmax [OPTIONAL] : maximum size of the rehashed table
! **************************************************************************************************
   RECURSIVE SUBROUTINE nl_hash_table_rehash(hash_table, nmax)
      TYPE(nl_hash_table_obj), INTENT(INOUT)             :: hash_table
      INTEGER, INTENT(IN), OPTIONAL                      :: nmax

      INTEGER                                            :: ii, my_nmax
      TYPE(nl_hash_table_element), ALLOCATABLE, &
         DIMENSION(:)                                    :: tmp_table

      IF (.NOT. nl_hash_table_has_data(hash_table)) THEN
         CALL nl_hash_table_create(hash_table, nmax)
         RETURN
      END IF
      IF (PRESENT(nmax)) THEN
         my_nmax = MAX(nmax, hash_table%obj%nelements)
      ELSE
         my_nmax = hash_table%obj%nmax
      END IF
      ALLOCATE (tmp_table(hash_table%obj%nmax))
      tmp_table(:) = hash_table%obj%table(:)
      CALL nl_hash_table_release(hash_table)
      CALL nl_hash_table_create(hash_table=hash_table, nmax=my_nmax)
      DO ii = 1, SIZE(tmp_table)
         IF (tmp_table(ii)%key .NE. EMPTY_KEY) THEN
            CALL nl_hash_table_add(hash_table=hash_table, &
                                   key=tmp_table(ii)%key, &
                                   val=tmp_table(ii)%val)
         END IF
      END DO
      DEALLOCATE (tmp_table)
   END SUBROUTINE nl_hash_table_rehash

! **************************************************************************************************
!> \brief releases the hash table.  Note that deallocating tasks stored in the table
!>        is the responsibility of the caller
!> \param hash_table : the nl_hash_table object in question
! **************************************************************************************************
   SUBROUTINE nl_hash_table_release(hash_table)
      TYPE(nl_hash_table_obj), INTENT(INOUT)             :: hash_table

      IF (ASSOCIATED(hash_table%obj)) THEN
         IF (ASSOCIATED(hash_table%obj%table)) THEN
            DEALLOCATE (hash_table%obj%table)
         END IF
         DEALLOCATE (hash_table%obj)
      ELSE
         NULLIFY (hash_table%obj)
      END IF
   END SUBROUTINE nl_hash_table_release

! **************************************************************************************************
!> \brief outputs the current information about the table
!> \param hash_table : the nl_hash_table object in question
!> \param nelements  : number of non-empty slots in the table
!> \param nmax       : maximum number of slots in the table
!> \param prime      : the prime used in the hash function
! **************************************************************************************************
   SUBROUTINE nl_hash_table_status(hash_table, nelements, nmax, prime)
      TYPE(nl_hash_table_obj), INTENT(INOUT)             :: hash_table
      INTEGER, INTENT(OUT), OPTIONAL                     :: nelements, nmax, prime

      LOGICAL                                            :: check_ok

      check_ok = nl_hash_table_has_data(hash_table)
      CPASSERT(check_ok)
      IF (PRESENT(nelements)) nelements = hash_table%obj%nelements
      IF (PRESENT(nmax)) nmax = hash_table%obj%nmax
      IF (PRESENT(prime)) prime = hash_table%obj%prime
   END SUBROUTINE nl_hash_table_status

! **************************************************************************************************
!> \brief Linear probing algorithm for the hash table
!> \param hash_table : the nl_hash_table object
!> \param key        : key to locate
!> \return : slot location in the table correspond to key, 0 if key not found
! **************************************************************************************************
   PURE FUNCTION nl_hash_table_linear_probe(hash_table, key) RESULT(islot)
      TYPE(nl_hash_table_obj), INTENT(IN)                :: hash_table
      INTEGER(KIND=int_8), INTENT(IN)                    :: key
      INTEGER                                            :: islot

      INTEGER                                            :: guess

      ! first guess is mapped by the hash_function
      guess = nl_hash_table_hash_function(hash_table, key)

      ! then search for key and stop at first empty slot from guess to
      ! nmax.  using the same linear probe for adding and retrieving
      ! makes all non-empty keys being put before the first empty slot.
      DO islot = guess, hash_table%obj%nmax
         IF ((hash_table%obj%table(islot)%key == key) .OR. &
             (hash_table%obj%table(islot)%key == EMPTY_KEY)) RETURN
      END DO

      ! if unsuccessful, search from 1 to guess
      DO islot = 1, guess - 1
         IF ((hash_table%obj%table(islot)%key == key) .OR. &
             (hash_table%obj%table(islot)%key == EMPTY_KEY)) RETURN
      END DO

      ! if not found and table is full set islot to 0
      islot = 0
   END FUNCTION nl_hash_table_linear_probe

! **************************************************************************************************
!> \brief Hash function
!> \param hash_table : the nl_hash_table object
!> \param key        : key to locate
!> \return : slot location in the table correspond to key, 0 if key not found
! **************************************************************************************************
   PURE FUNCTION nl_hash_table_hash_function(hash_table, key) RESULT(hash)
      TYPE(nl_hash_table_obj), INTENT(IN)                :: hash_table
      INTEGER(KIND=int_8), INTENT(IN)                    :: key
      INTEGER                                            :: hash

      INTEGER(KIND=int_8)                                :: hash_8, nmax_8, prime_8

      nmax_8 = INT(hash_table%obj%nmax, int_8)
      prime_8 = INT(hash_table%obj%prime, int_8)

      ! IAND with nmax-1 is equivalent to MOD nmax if nmax is alway a power of 2.
      hash_8 = IAND(key*prime_8, nmax_8 - 1) + 1_int_8
      hash = INT(hash_8)
   END FUNCTION nl_hash_table_hash_function

END MODULE qs_nl_hash_table_types

