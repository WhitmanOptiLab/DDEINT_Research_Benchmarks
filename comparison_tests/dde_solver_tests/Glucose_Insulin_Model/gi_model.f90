PROGRAM gi_model

! Glucose-Insulin Regulation Model solved with DDE_SOLVER.
!
! The DDE is defined in the module define_DDEs. The problem
! is solved here with DDE_SOLVER and its output written to
! a file.

  USE gi_define_DDEs
  USE DDE_SOLVER_M

  IMPLICIT NONE

  INTEGER, DIMENSION(2) :: NVAR = (/NEQN,NLAGS/)

  ! Constant delay and history:
  DOUBLE PRECISION, DIMENSION(NLAGS) :: DELAYS  = (/ 5.0D0 /)
  DOUBLE PRECISION, DIMENSION(NEQN)  :: HISTORY = (/ 1.0D0, 1.0D0 /)

  TYPE(DDE_SOL)  :: SOL
  TYPE(DDE_OPTS) :: OPTS

  INTEGER :: I,J

  tau = 5.0D0

  SOL = DDE_SOLVER(NVAR,DDES,DELAYS,HISTORY,TSPAN=(/ 0.0D0, 200.0D0 /))

  IF (SOL%FLAG == 0) THEN

    OPEN(UNIT=8, FILE='data/csv_files/gi_model_dde_solver.csv')
    WRITE(UNIT=8, FMT='(A)') 'Time,Glucose,Insulin'
    DO I = 1, SOL%NPTS
      WRITE(UNIT=8, FMT='(F14.10, 2(",", F14.10))') SOL%T(I), (SOL%Y(I,J), J=1,NEQN)
    END DO

    PRINT *,' Normal return from DDE_SOLVER with results'
    PRINT *," written to the file 'data/csv_files/gi_model_dde_solver.csv'."

  ELSE
    PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
    SOL%FLAG
  END IF

  CALL PRINT_STATS(SOL)
  CALL RELEASE_ARRAYS(SOL,OPTS)

END PROGRAM gi_model