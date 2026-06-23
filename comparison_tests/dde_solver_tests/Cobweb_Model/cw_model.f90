PROGRAM cw_model

! Cobweb Model solved with DDE_SOLVER.
!
! The DDE is defined in the module define_DDEs. The problem
! is solved here with DDE_SOLVER and its output written to
! a file.

  USE cw_define_DDEs
  USE DDE_SOLVER_M

  IMPLICIT NONE

  INTEGER, DIMENSION(2) :: NVAR = (/NEQN,NLAGS/)

  ! Constant delay and history:
  DOUBLE PRECISION, DIMENSION(NLAGS) :: DELAYS  = (/ 1.5D0 /)
  DOUBLE PRECISION, DIMENSION(NEQN)  :: HISTORY = (/ 0.4D0, 0.4D0 /)

  INTEGER, PARAMETER :: NOUT = 10001
  DOUBLE PRECISION, DIMENSION(NOUT) :: TSPAN
  INTEGER :: I, J

  TYPE(DDE_SOL)  :: SOL
  TYPE(DDE_OPTS) :: OPTS

  tau = 1.5D0
  
  ! Build the output time array
  DO I = 1, NOUT
    TSPAN(I) = (I-1) * 1.0D0
  END DO

  SOL = DDE_SOLVER(NVAR,DDES,DELAYS,HISTORY,TSPAN=TSPAN)

  IF (SOL%FLAG == 0) THEN

    OPEN(UNIT=8, FILE='data/csv_files/cw_model_dde_solver.csv')
    WRITE(UNIT=8, FMT='(A)') 'Time,price,expected_price'
    DO I = 1, SOL%NPTS
      WRITE(UNIT=8, FMT='(F14.10, 2(",", F14.10))') SOL%T(I), (SOL%Y(I,J), J=1,NEQN)
    END DO

    PRINT *,' Normal return from DDE_SOLVER with results'
    PRINT *," written to the file 'data/csv_files/cw_model_dde_solver.csv'."

  ELSE
    PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
    SOL%FLAG
  END IF

  CALL PRINT_STATS(SOL)
  CALL RELEASE_ARRAYS(SOL,OPTS)

END PROGRAM cw_model