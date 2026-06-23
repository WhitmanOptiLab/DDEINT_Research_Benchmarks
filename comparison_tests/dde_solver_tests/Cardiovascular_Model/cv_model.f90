PROGRAM cv_model

! Cardiovascular Model solved with DDE_SOLVER.
!
! The DDE is defined in the module define_DDEs. The problem
! is solved here with DDE_SOLVER and its output written to
! a file.

  USE cv_define_DDEs
  USE DDE_SOLVER_M

  IMPLICIT NONE

  INTEGER, DIMENSION(2) :: NVAR = (/NEQN,NLAGS/)

  ! Constant delay and history:
  DOUBLE PRECISION, DIMENSION(NLAGS) :: DELAYS = (/ 4.0D0 /)
  DOUBLE PRECISION, DIMENSION(NEQN) :: HISTORY = (/ 93.0D0, &
    (1.0D0/(1.0D0 + 1.05D0/0.068D0))*93.0D0, &
    (1.0D0/(1.05D0*67.9D0))*(1.0D0/(1.0D0 + 0.068D0/1.05D0))*93.0D0 /)

  INTEGER, PARAMETER :: NOUT = 50001
  DOUBLE PRECISION, DIMENSION(NOUT) :: TSPAN
  INTEGER :: I, J

  TYPE(DDE_SOL)  :: SOL
  TYPE(DDE_OPTS) :: OPTS

  tau = 4.0D0

  ! Build the output time array
  DO I = 1, NOUT
    TSPAN(I) = (I-1) * 0.2D0
  END DO

  SOL = DDE_SOLVER(NVAR,DDES,DELAYS,HISTORY,TSPAN=TSPAN)

  IF (SOL%FLAG == 0) THEN

    OPEN(UNIT=8, FILE='data/csv_files/cv_model_dde_solver.csv')
    WRITE(UNIT=8, FMT='(A)') 'Time,Pa,Pv,HR'
    DO I = 1, SOL%NPTS
      WRITE(UNIT=8, FMT='(F14.10, 3(",", F14.10))') SOL%T(I), (SOL%Y(I,J), J=1,NEQN)
    END DO

    PRINT *,' Normal return from DDE_SOLVER with results'
    PRINT *," written to the file 'data/csv_files/cv_model_dde_solver.csv'."

  ELSE
    PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
    SOL%FLAG
  END IF

  CALL PRINT_STATS(SOL)
  CALL RELEASE_ARRAYS(SOL,OPTS)

END PROGRAM cv_model