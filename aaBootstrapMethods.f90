! Utilities for evaluating errors with the Bootstrap method
! Written by Anna Alemany

! creat: 11/10/12

MODULE aaBootstrapMethods
  USE nrtype; USE nrutil, ONLY : nrerror, assert_eq
  USE ran_state!, ONLY: K4B,amm,lenran,ran_init, &
!   iran0,jran0,kran0,nran0,mran0,rans
  IMPLICIT NONE

  INTERFACE BSran1
    MODULE PROCEDURE BSran1_s, BSran1_v
  END INTERFACE

  CONTAINS

! Average and variances

FUNCTION av(f)
  REAL(SP), DIMENSION(:), INTENT(IN) :: f
  REAL(SP) :: av
  INTEGER(I4B) :: sz
  sz = size(f)
  av = sum(f)/real(sz,sp)
END FUNCTION av

FUNCTION var(f)
  REAL(SP), DIMENSION(:), INTENT(IN) :: f
  REAL(SP) :: var
  REAL(SP) :: p
  INTEGER(I4B) :: sz
  sz = size(f)
  p = av(f)
  var = sum((f-p)**2)/real(sz-1,sp) ! bootstrap
END FUNCTION var

! Bootstrap Sampling 

SUBROUTINE BSSampling(f,fsort)
  REAL(SP), DIMENSION(:), INTENT(IN) :: f
  REAL(SP), DIMENSION(:), INTENT(OUT) :: fsort
! Reordena el vector f en el fsort
  INTEGER(I4B) :: sz, i
  REAL(SP), DIMENSION(:), ALLOCATABLE :: rnd
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: rndint
  sz = assert_eq(size(f), size(fsort), 'Err BSSampling')
  allocate(rnd(sz),rndint(sz))
  call BSran1(rnd)
  rndint = sz*rnd+1
!   print *, rndint
  do i = 1, sz
    fsort(i) = f(rndint(i))
  end do
END SUBROUTINE BSSampling

! Filtrar dades: promitjar dades per fer vectors més petits
FUNCTION Filter(vector, filtre, sz1)
  REAL(SP), DIMENSION(:), INTENT(IN) :: vector
  INTEGER(I4B), INTENT(IN) :: filtre, sz1
  REAL(SP), DIMENSION(sz1) :: Filter
! ! !   
  INTEGER(I4B) :: sz0, i
  sz0 = size(vector)
  do i = 1, sz1
    if (filtre*i <= sz0) then
      Filter(i) = av(vector(filtre*(i-1)+1:filtre*i))
    else
      Filter(i) = av(vector(filtre*(i-1)+1:sz0))
    end if
  end do
  
END FUNCTION

! Numerical Recipes, adaptations

SUBROUTINE BSran1_s(harvest)
!   USE nrtype
!   USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
!   iran0,jran0,kran0,nran0,mran0,rans
!   IMPLICIT NONE
  REAL(SP), INTENT(OUT) :: harvest
  if (lenran < 1) call ran_init(1)
  rans=iran0-kran0
  if (rans < 0) rans=rans+2147483579_k4b
  iran0=jran0
  jran0=kran0
  kran0=rans
  nran0=ieor(nran0,ishft(nran0,13))
  nran0=ieor(nran0,ishft(nran0,-17))
  nran0=ieor(nran0,ishft(nran0,5))
  if (nran0 == 1) nran0=270369_k4b
  mran0=ieor(mran0,ishft(mran0,5))
  mran0=ieor(mran0,ishft(mran0,-13))
  mran0=ieor(mran0,ishft(mran0,6))
  rans=ieor(nran0,rans)+mran0
  harvest=amm*merge(rans,not(rans), rans<0 )
END SUBROUTINE BSran1_s

SUBROUTINE BSran1_v(harvest)
!   USE nrtype
!   USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
!   iran,jran,kran,nran,mran,ranv
!   IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
  INTEGER(K4B) :: n
  n=size(harvest)
  if (lenran < n+1) call ran_init(n+1)
  ranv(1:n)=iran(1:n)-kran(1:n)
  where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
  iran(1:n)=jran(1:n)
  jran(1:n)=kran(1:n)
  kran(1:n)=ranv(1:n)
  nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
  nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
  nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
  where (nran(1:n) == 1) nran(1:n)=270369_k4b
  mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),5))
  mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),-13))
  mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),6))
  ranv(1:n)=ieor(nran(1:n),ranv(1:n))+mran(1:n)
  harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
END SUBROUTINE BSran1_v

END MODULE aaBootstrapMethods