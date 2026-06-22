MODULE cw_define_DDEs

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
    DOUBLE PRECISION :: a=-0.5D0, b=1.5D0, c1=0.8D0, d=0.0D0, &
                        beta=0.6D0, speed=1.2D0

    ! Local variables
    DOUBLE PRECISION :: p_now, p_exp_delayed, demand, supply

    ! Local solution variables and delayed
    p_now         = Y(1)
    p_exp_delayed = Z(2,1)

    ! Derivatives:
    demand = a * p_now + b
    supply = c1 * p_exp_delayed + d

    DY(1) = speed * (demand - supply)
    DY(2) = beta * (p_now - Y(2))

    RETURN 
  END SUBROUTINE DDES
END MODULE cw_define_DDEs