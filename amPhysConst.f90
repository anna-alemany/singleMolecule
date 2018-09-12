! Module with the updated values of various physical constants.
! Written by Alessandro Mossa.

MODULE amPhysConst
  USE nrtype
  IMPLICIT NONE

  REAL(SP), PARAMETER :: k_B = 0.013806504_sp ! pN x nm / K
  REAL(SP), PARAMETER :: T0 = -273.15_sp ! C

CONTAINS

  FUNCTION fric(r,T,mol)
    REAL(SP), INTENT(IN) :: r, T, mol
    REAL(SP) :: fric
    ! compute the attrition coefficient of the bead given the temperature,
    ! the radius of the bead and the NaCl concentration. If r is given in nm, 
    ! and T in K, then the friction coefficient is expressed in (pN*s)/nm
    ! Reference: J. Chem. Eng. Data 1996, 41, 516-520.
    REAL(SP) :: eta
    REAL(SP), PARAMETER :: A = 2.414E-5 ! Pa*s
    REAL(SP), PARAMETER :: B = 247.8    ! K
    REAL(SP), PARAMETER :: C = 140.0    ! K
    REAL(SP), PARAMETER :: AA = 0.0061
    REAL(SP), PARAMETER :: BB = 0.0794
    REAL(SP), PARAMETER :: DD = 0.01142
    REAL(SP), PARAMETER :: EE = 0.000619
    eta = A*10**(B/(T-C))   ! viscosity of pure water, in Pa*s
    eta = eta*(1.0_sp+AA*sqrt(mol)+BB*mol+DD*mol**2+EE*mol**3.5_sp)
    fric = 3.0E-6_sp*TWOPI*eta*r
  END FUNCTION fric

  SUBROUTINE timeconvert(time,days,hours,minutes,seconds)
    REAL(SP), INTENT(IN) :: time
    INTEGER(I4B), INTENT(OUT) :: days, hours, minutes
    REAL(SP), INTENT(OUT) :: seconds
    ! convert the amount of time (given in seconds) into human readable units
    REAL(SP) :: tm
    tm = time
    days = floor(tm/86400.0_sp)
    tm = modulo(tm,86400.0_sp)
    hours = floor(tm/3600.0_sp)
    tm = modulo(tm,3600.0_sp)
    minutes = floor(tm/60.0_sp)
    seconds = modulo(tm,60.0_sp)
  END SUBROUTINE timeconvert

END MODULE amPhysConst
