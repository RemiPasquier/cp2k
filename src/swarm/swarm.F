!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Swarm-framwork, provides a convenient master/worker architecture.
!> \author Ole Schuett
! **************************************************************************************************
MODULE swarm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_unit_nr
   USE global_types,                    ONLY: global_environment_type
   USE input_section_types,             ONLY: section_type,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length
   USE message_passing,                 ONLY: mp_para_env_type
   USE swarm_master,                    ONLY: swarm_master_finalize,&
                                              swarm_master_init,&
                                              swarm_master_steer,&
                                              swarm_master_type
   USE swarm_message,                   ONLY: swarm_message_add,&
                                              swarm_message_free,&
                                              swarm_message_get,&
                                              swarm_message_type
   USE swarm_mpi,                       ONLY: swarm_mpi_finalize,&
                                              swarm_mpi_init,&
                                              swarm_mpi_recv_command,&
                                              swarm_mpi_recv_report,&
                                              swarm_mpi_send_command,&
                                              swarm_mpi_send_report,&
                                              swarm_mpi_type
   USE swarm_worker,                    ONLY: swarm_worker_execute,&
                                              swarm_worker_finalize,&
                                              swarm_worker_init,&
                                              swarm_worker_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'swarm'

   PUBLIC :: run_swarm

CONTAINS

! **************************************************************************************************
!> \brief Central driver routine of the swarm framework, called by cp2k_runs.F
!> \param input_declaration ...
!> \param root_section ...
!> \param para_env ...
!> \param globenv ...
!> \param input_path ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE run_swarm(input_declaration, root_section, para_env, globenv, input_path)
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(section_vals_type), POINTER                   :: root_section
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_environment_type), POINTER             :: globenv
      CHARACTER(LEN=*), INTENT(IN)                       :: input_path

      CHARACTER(len=*), PARAMETER                        :: routineN = 'run_swarm'

      INTEGER                                            :: handle, iw, n_workers
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iw = cp_print_key_unit_nr(logger, root_section, &
                                "SWARM%PRINT%MASTER_RUN_INFO", extension=".masterLog")

      IF (iw > 0) WRITE (iw, "(A)") " SWARM| Ready to roll :-)"

      CALL section_vals_val_get(root_section, "SWARM%NUMBER_OF_WORKERS", &
                                i_val=n_workers)

      IF (n_workers == 1 .AND. para_env%num_pe == 1) THEN
         IF (iw > 0) WRITE (iw, "(A)") " SWARM| Running in single worker mode."
         CALL swarm_serial_driver(input_declaration, root_section, input_path, para_env, globenv)
      ELSE
         IF (iw > 0) WRITE (iw, "(A)") " SWARM| Running in master / workers mode."
         !printkey iw passed on for output from swarm_mpi_init()
         CALL swarm_parallel_driver(n_workers, input_declaration, root_section, input_path, para_env, globenv, iw)
      END IF

      CALL timestop(handle)
   END SUBROUTINE run_swarm

! **************************************************************************************************
!> \brief Special driver for using only a single worker.
!> \param input_declaration ...
!> \param root_section ...
!> \param input_path ...
!> \param para_env ...
!> \param globenv ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_serial_driver(input_declaration, root_section, input_path, para_env, globenv)
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(section_vals_type), POINTER                   :: root_section
      CHARACTER(LEN=*), INTENT(IN)                       :: input_path
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_environment_type), POINTER             :: globenv

      INTEGER                                            :: handle
      LOGICAL                                            :: should_stop
      TYPE(swarm_master_type)                            :: master
      TYPE(swarm_message_type)                           :: cmd, report
      TYPE(swarm_worker_type)                            :: worker

      CALL swarm_master_init(master, para_env, globenv, root_section, n_workers=1)
      CALL swarm_worker_init(worker, para_env, input_declaration, root_section, &
                             input_path, worker_id=1)

      CALL swarm_message_add(report, "worker_id", 1)
      CALL swarm_message_add(report, "status", "initial_hello")

      should_stop = .FALSE.
      DO WHILE (.NOT. should_stop)
         CALL timeset("swarm_worker_await_reply", handle)
         CALL swarm_master_steer(master, report, cmd)
         CALL timestop(handle)
         CALL swarm_message_free(report)
         CALL swarm_worker_execute(worker, cmd, report, should_stop)
         CALL swarm_message_free(cmd)
      END DO

      CALL swarm_message_free(report)
      CALL swarm_worker_finalize(worker)
      CALL swarm_master_finalize(master)

   END SUBROUTINE swarm_serial_driver

! **************************************************************************************************
!> \brief Normal driver routine for parallelized runs.
!> \param n_workers ...
!> \param input_declaration ...
!> \param root_section ...
!> \param input_path ...
!> \param para_env ...
!> \param globenv ...
!> \param iw ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_parallel_driver(n_workers, input_declaration, root_section, input_path, para_env, globenv, iw)
      INTEGER, INTENT(IN)                                :: n_workers
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(section_vals_type), POINTER                   :: root_section
      CHARACTER(LEN=*), INTENT(IN)                       :: input_path
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_environment_type), POINTER             :: globenv
      INTEGER, INTENT(IN)                                :: iw

      INTEGER                                            :: worker_id
      TYPE(swarm_mpi_type)                               :: swarm_mpi

      CALL swarm_mpi_init(swarm_mpi, para_env, root_section, n_workers, worker_id, iw)

      IF (ASSOCIATED(swarm_mpi%worker)) THEN
         CALL swarm_parallel_worker_driver(swarm_mpi, input_declaration, worker_id, root_section, input_path)
      ELSE
         CALL swarm_parallel_master_driver(swarm_mpi, n_workers, root_section, globenv)
      END IF

      CALL swarm_mpi_finalize(swarm_mpi, root_section)

   END SUBROUTINE swarm_parallel_driver

! **************************************************************************************************
!> \brief Worker's driver routine for parallelized runs.
!> \param swarm_mpi ...
!> \param input_declaration ...
!> \param worker_id ...
!> \param root_section ...
!> \param input_path ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_parallel_worker_driver(swarm_mpi, input_declaration, worker_id, root_section, input_path)
      TYPE(swarm_mpi_type), INTENT(IN)                   :: swarm_mpi
      TYPE(section_type), POINTER                        :: input_declaration
      INTEGER, INTENT(IN)                                :: worker_id
      TYPE(section_vals_type), POINTER                   :: root_section
      CHARACTER(LEN=*), INTENT(IN)                       :: input_path

      INTEGER                                            :: handle
      LOGICAL                                            :: should_stop
      TYPE(swarm_message_type)                           :: cmd, report
      TYPE(swarm_worker_type)                            :: worker

      CALL swarm_worker_init(worker, swarm_mpi%worker, input_declaration, &
                             root_section, input_path, worker_id=worker_id)

      CALL swarm_message_add(report, "worker_id", worker_id)
      CALL swarm_message_add(report, "status", "initial_hello")

      should_stop = .FALSE.
      DO WHILE (.NOT. should_stop)
         CALL timeset("swarm_worker_await_reply", handle)
         CALL swarm_mpi_send_report(swarm_mpi, report)
         CALL swarm_message_free(report)
         CALL swarm_mpi_recv_command(swarm_mpi, cmd)
         CALL timestop(handle)
         CALL swarm_worker_execute(worker, cmd, report, should_stop)
         CALL swarm_message_free(cmd)
      END DO

      CALL swarm_message_free(report)
      CALL swarm_worker_finalize(worker)

   END SUBROUTINE swarm_parallel_worker_driver

! **************************************************************************************************
!> \brief Master's driver routine for parallelized runs.
!> \param swarm_mpi ...
!> \param n_workers ...
!> \param root_section ...
!> \param globenv ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_parallel_master_driver(swarm_mpi, n_workers, root_section, globenv)
      TYPE(swarm_mpi_type), INTENT(IN)                   :: swarm_mpi
      INTEGER, INTENT(IN)                                :: n_workers
      TYPE(section_vals_type), POINTER                   :: root_section
      TYPE(global_environment_type), POINTER             :: globenv

      CHARACTER(len=default_string_length)               :: command
      INTEGER                                            :: i_shutdowns, j, wid
      LOGICAL, DIMENSION(n_workers)                      :: is_waiting
      TYPE(swarm_master_type)                            :: master
      TYPE(swarm_message_type)                           :: cmd, report

      is_waiting(:) = .FALSE.

      CALL swarm_master_init(master, swarm_mpi%master, globenv, root_section, n_workers)

      i_shutdowns = 0
      j = 0

      DO WHILE (i_shutdowns < n_workers)
         ! Each iteration if the loop does s.th. different depending on j.
         ! First (j==0) it receives one report with (blocking) MPI,
         ! then it searches through the list is_waiting.
         j = MOD(j + 1, n_workers + 1)
         IF (j == 0) THEN
            CALL swarm_mpi_recv_report(swarm_mpi, report)
         ELSE IF (is_waiting(j)) THEN
            is_waiting(j) = .FALSE.
            CALL swarm_message_add(report, "worker_id", j)
            CALL swarm_message_add(report, "status", "wait_done")
         ELSE
            CYCLE
         END IF

         CALL swarm_master_steer(master, report, cmd)
         CALL swarm_message_free(report)

         CALL swarm_message_get(cmd, "command", command)
         IF (TRIM(command) == "wait") THEN
            CALL swarm_message_get(cmd, "worker_id", wid)
            is_waiting(wid) = .TRUE.
         ELSE
            CALL swarm_mpi_send_command(swarm_mpi, cmd)
            IF (TRIM(command) == "shutdown") i_shutdowns = i_shutdowns + 1
         END IF
         CALL swarm_message_free(cmd)
      END DO

      CALL swarm_master_finalize(master)

   END SUBROUTINE swarm_parallel_master_driver

END MODULE swarm

