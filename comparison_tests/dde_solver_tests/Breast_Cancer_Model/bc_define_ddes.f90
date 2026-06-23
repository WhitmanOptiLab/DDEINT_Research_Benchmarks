MODULE bc_define_DDEs

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
    DOUBLE PRECISION :: p0=0.2D0, q0=0.3D0, v0=1.0D0, d0=5.0D0, &
                        p1=0.2D0, q1=0.3D0, v1=1.0D0, d1=1.0D0, &
                        d2=1.0D0, beta0=1.0D0, beta1=1.0D0

    ! Local variables
    DOUBLE PRECISION :: u1, u2, u3, hist3, denom0, denom1, &
                        du1dt, du2dt, du3dt

    ! Local solution variables and delayed 
    u1 = Y(1)
    u2 = Y(2)
    u3 = Y(3)
    hist3 = Z(3,1)

    ! Denominantors
    denom0 = 1.0D0 + beta0 * hist3**2
    denom1 = 1.0D0 + beta1 * hist3**2

    ! Local derivatives:
    du1dt = (v0/(1.0D0 + beta0 * hist3**2)) * (p0 - q0)*u1 - d0*u1
    du2dt = (v0/(1.0D0 + beta0 * hist3**2)) * (1.0D0 - p0 + q0)*u1 + &
            (v1/(1.0D0 + beta1 * hist3**2)) * (p1 - q1)*u2 - d1*u2
    du3dt = (v1/(1.0D0 + beta1 * hist3**2)) * (1.0D0 - p1 + q1)*u2 - &
            d2*u3

    ! Derivatives:
    DY = (/ du1dt, du2dt, du3dt/)

    RETURN 
  END SUBROUTINE DDES
END MODULE bc_define_DDEs