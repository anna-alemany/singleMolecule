! Module que conté els paràmetres usats per analitzar hairpins.

MODULE HairpinParameters
  USE nrtype
  IMPLICIT NONE

! Temperatura
  REAL(SP), PARAMETER :: T=24._sp  !ºC
  REAL(SP), PARAMETER :: kb = 0.013806503_sp, kT = kb*(T+273.15)

! Informació general del hairpin
  REAL(SP), PARAMETER :: d0 = 2. !nm
  REAL(SP), PARAMETER :: a = 0.59 ! bases/nm
  INTEGER(SP), PARAMETER :: nloop = 4 ! bases loop
  INTEGER(SP), PARAMETER :: bp = 30 ! nombre base-pairs in stem
  REAL(SP), PARAMETER :: lcont = 2*a*real(bp,sp)+a*real(nloop,sp)

! Model de WLC
  REAL(SP), PARAMETER :: P = 1. !nm, Persistence length
! Model de FJC
  REAL(SP), PARAMETER :: Y = 1000. ! pN; Modul Young
  REAL(SP), PARAMETER :: b = 1.105 ! nm, Kuhn length

END MODULE HairpinParameters
