PROGRAM bc_model

! Breast Cancer Model solved with DDE_SOLVER.
!
! The DDE is defined in the module define_DDEs. The problem
! is solved here with DDE_SOLVER and its output written to
! a file.

  USE bc_define_DDEs
  USE DDE_SOLVER_M

  IMPLICIT NONE

  ! The quantities
  !
  !   NEQN  = number of equations
  !   NLAGS = number of delays
  !
  ! are defined in the module define_DDEs as PARAMETERs so
  ! they can be used for dimensioning arrays here. They are
  ! passed to the solver in the array NVAR.
  INTEGER, DIMENSION(2) :: NVAR = (/NEQN,NLAGS/)

  ! Constant delay and history:
  DOUBLE PRECISION, DIMENSION(NLAGS) :: DELAYS = (/ 1.0D0 /)
  DOUBLE PRECISION, DIMENSION(NEQN)  :: HISTORY = (/ 1.0D0, 1.0D0, 1.0D0 /)

  INTEGER, PARAMETER :: NOUT = 501
  DOUBLE PRECISION, DIMENSION(NOUT) :: TSPAN
  INTEGER :: I, J

  TYPE(DDE_SOL)  :: SOL
  TYPE(DDE_OPTS) :: OPTS

  tau = 1.0D0

  ! Build the output time array
  DO I = 1, NOUT
    TSPAN(I) = (I-1) * 0.02D0
  END DO

  SOL = DDE_SOLVER(NVAR,DDES,DELAYS,HISTORY,TSPAN=TSPAN)

  ! Was the solver successful?
  IF (SOL%FLAG == 0) THEN

    ! Write the solution to a file for subsequent plotting.
    OPEN(UNIT=8, FILE='data/csv_files/bc_model_dde_solver.csv')
    WRITE(UNIT=8, FMT='(A)') 'Time,u1,u2,u3'
    DO I = 1, SOL%NPTS
      WRITE(UNIT=8, FMT='(F12.8, 3(",", F12.8))') SOL%T(I), (SOL%Y(I,J), J=1,NEQN)
    END DO

    PRINT *,' Normal return from DDE_SOLVER with results'
    PRINT *," written to the file 'data/csv_files/bc_model_dde_solver.csv'."

  ELSE
    PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
    SOL%FLAG
  END IF

  ! Integration statistics:
  CALL PRINT_STATS(SOL)
  CALL RELEASE_ARRAYS(SOL,OPTS)

END PROGRAM bc_model