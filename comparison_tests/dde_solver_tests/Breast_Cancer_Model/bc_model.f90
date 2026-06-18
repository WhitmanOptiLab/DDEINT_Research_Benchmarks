PROGRAM bc_model

! Breast Cancer Model solved with DDE_SOLVER.
!
! The DDE is defined in the module define_DDEs. The problem
! is solved here with DDE_SOLVER and its output written to
! a file.

  USE define_DDEs
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

  TYPE(DDE_SOL)  :: SOL
  TYPE(DDE_OPTS) :: OPTS

  !   SOL%NPTS         -- NPTS, number of output points.
  !
  !   SOL%T(NPTS)      -- values of independent variable, T.
  !
  !   SOL%Y(NPTS,NEQN) -- values of dependent variable, Y,
  !                       corresponding to values of SOL%T.

  INTEGER :: I,J

  ! tau is a global variable in define_DDEs.
  ! Assign its value here.
  tau = 1.0D0

  SOL = DDE_SOLVER(NVAR,DDES,DELAYS,HISTORY,TSPAN=(/ 0.0D0, 10.0D0 /))

  ! Was the solver successful?
  IF (SOL%FLAG == 0) THEN

    ! Write the solution to a file for subsequent plotting.
    OPEN(UNIT=8, FILE='bc_model.dat')
    DO I = 1,SOL%NPTS
      WRITE(UNIT=8,FMT='(4D12.4)') SOL%T(I),(SOL%Y(I,J),J=1,NEQN)
    ENDDO
    CLOSE(UNIT=8)

    PRINT *,' Normal return from DDE_SOLVER with results'
    PRINT *," written to the file 'data/bc_model.dat'."

  ELSE
    PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
    SOL%FLAG
  END IF

  ! Integration statistics:
  CALL PRINT_STATS(SOL)
  CALL RELEASE_ARRAYS(SOL,OPTS)

END PROGRAM bc_model