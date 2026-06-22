MODULE gi_define_DDEs

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NEQN=2,NLAGS=1

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
    DOUBLE PRECISION :: V_G=10.0D0, V_I=5.0D0, S_I=1.0D0, G_b=90.0D0, &
                        I_b=10.0D0, n1=0.1D0, gamma=0.05D0

    ! Local variables
    DOUBLE PRECISION :: G_Tau, I_Tau

    ! Delayed values
    G_Tau = Z(1,1)
    I_Tau = Z(2,1)

    ! Derivatives:
    DY(1) = (Y(1) - G_b) / V_G - (S_I * I_Tau * Y(1)) / V_G + 10.0D0
    DY(2) = -n1 * (Y(2) - I_b) + (gamma * G_Tau) / V_I

    RETURN 
  END SUBROUTINE DDES
END MODULE gi_define_DDEs