MODULE cv_define_DDEs

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NEQN=3,NLAGS=1

  !! Physical parameters assigned in the main program
  DOUBLE PRECISION :: tau


CONTAINS

    SUBROUTINE DDES(T,Y,Z,DY)
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(:) :: Y,DY
    DOUBLE PRECISION, DIMENSION(:,:) :: Z
    INTENT (IN) :: T, Y, Z
    INTENT (OUT) :: DY
    ! Physical parameters
    DOUBLE PRECISION :: ca=1.55D0, cv=519.0D0, r1=0.068D0, &
                        Vstr=67.9D0, alphas=93.0D0, alphap=93.0D0, &
                        alphaH=0.84D0, betas=7.0D0, betap=7.0D0, &
                        betaH=1.17D0, gammaH=0.0D0
    ! Time-varying R
    DOUBLE PRECISION :: R_t
    ! Local variables
    DOUBLE PRECISION :: Paoft, Pvoft, Hoft, Patau, Ts, Tp

        IF (T <= 600.0D0) THEN
            R_t = 1.05D0
        ELSE
            R_t = 0.21D0 * EXP(600.0D0 - T) + 0.84D0
        END IF

        Paoft = Y(1)
        Pvoft = Y(2)
        Hoft  = Y(3)
        Patau = Z(1,1)

        DY(1) = -(1.0D0/(ca*R_t))*Paoft + (1.0D0/(ca*R_t))*Pvoft + (1.0D0/ca)*Vstr*Hoft
        DY(2) =  (1.0D0/(cv*R_t))*Paoft - (1.0D0/(cv*R_t) + 1.0D0/(cv*r1))*Pvoft
        Ts    = 1.0D0 / (1.0D0 + (Patau/alphas)**betas)
        Tp    = 1.0D0 / (1.0D0 + (alphap/Paoft)**betap)
        DY(3) = (alphaH*Ts) / (1.0D0 + gammaH*Tp) - betaH*Tp
    RETURN
    END SUBROUTINE DDES

END MODULE cv_define_DDEs