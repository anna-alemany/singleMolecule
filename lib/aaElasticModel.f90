! Utilities for evaluating elastic models
! Written by Anna Alemany

! creat: 08/07/10
! modificat: 16/12/10
! modificat: 14/01/11 (caràcter vectorial a algunes funcions)
! modificat: 26/04/11 (caràcter vectorial de FJC)

! WLC:
! xWLC -> Donada f, P i l, troba la x pel model de WLC
! dfdxWLC

! TC:
! fTC
! xTC

MODULE aaElasticModel
  USE nrtype; USE nrutil, ONLY : nrerror, arth
  IMPLICIT NONE

  INTERFACE xFJC
     MODULE PROCEDURE xFJC_s, xFJC_v
  END INTERFACE

  INTERFACE int_fdxFJC
    MODULE PROCEDURE int_fdxFJC_s, int_fdxFJC_v
  END INTERFACE

  INTERFACE int_xdfFJC
    MODULE PROCEDURE int_xdfFJC_s, int_xdfFJC_v
  END INTERFACE

  INTERFACE xWLC
    MODULE PROCEDURE xWLC_s, xWLC_v
  END INTERFACE 

  INTERFACE int_fdxWLC
    MODULE PROCEDURE int_fdxWLC_s, int_fdxWLC_v
  END INTERFACE

  INTERFACE int_xdfWLC
    MODULE PROCEDURE int_xdfWLC_s, int_xdfWLC_v
  END INTERFACE

  CONTAINS

! Worm-Like-Chain (WLC); Bustamante: f=kT/P*(0.25/(1-x/l)**2-0.25+x/l)

FUNCTION xWLC_s(f,P,l,kT) ! nm
  REAL(SP), INTENT(IN) :: f, P, l, kT
  REAL(SP) :: xWLC_s
! funció inversa de f=kT/P*(0.25/(1-x/l)**2-0.25+x/l)
  xWLC_s = rtsafeWLC(WLC,P,kT,f,0.,1.,1E-5)
  xWLC_s = l*(1.-xWLC_s)
END FUNCTION xWLC_s

FUNCTION xWLC_v(f,P,l,kT) ! nm
  REAL(SP), DIMENSION(:), INTENT(IN) :: f
  REAL(SP), INTENT(IN) :: P, l, kT
  REAL(SP), DIMENSION(size(f)) :: xWLC_v
  INTEGER(I4B) :: i, sz
! funció inversa de f=kT/P*(0.25/(1-x/l)**2-0.25+x/l)
  sz = size(f)
  do i = 1, sz
    xWLC_v(i) = rtsafeWLC(WLC,P,kT,f(i),0.,1.,1E-5)
    xWLC_v(i) = l*(1.-xWLC_v(i))
  end do
END FUNCTION xWLC_v

FUNCTION xWLC_Bouchiat(f,P,l,kT) ! nm
  REAL(SP), INTENT(IN) :: f, P, l, kT
  REAL(SP) :: xWLC_Bouchiat
! funció inversa de f=kT/P*(0.25/(1-x/l)**2-0.25+x/l)+an(x/l)^n
  xWLC_Bouchiat = rtsafeWLC_Bouchiat(f,P,l,kT,funcd_WLC_Bouchiat,0.,1.,1.e-5)
  xWLC_Bouchiat = l*xWLC_Bouchiat
END FUNCTION xWLC_Bouchiat

FUNCTION dfdxWLC(x,P,l,kT)
  REAL(SP), INTENT(IN) :: x, P, l, kT
  REAL(SP) :: dfdxWLC
  dfdxWLC = kT*(0.5/(1.-x/l)**3+1.)/(P*l)
END FUNCTION dfdxWLC

FUNCTION dfdxWLC_Bouchiat(x,P,l,kT)
  REAL(SP), INTENT(IN) :: x, P, l, kT
  REAL(SP) :: y
  REAL(SP) :: a2 = -0.5164228, a3 = -2.737418, a4 = 16.07497
  REAL(SP) :: a5 = -38.87607, a6 = 39.49944, a7 = -14.17718
  REAL(SP) :: dfdxWLC_Bouchiat
  y = x/l
  dfdxWLC_Bouchiat = 0.5/(1.-y)**3+1.+2*a2*y + 3*a3*y**2 + 4*a4*y**3 + 5*a5*y**4 + 6*a6*y**5 + 7*a7*y**6
  dfdxWLC_Bouchiat = kT*dfdxWLC_Bouchiat/(P*l)
END FUNCTION dfdxWLC_Bouchiat

FUNCTION int_xdfWLC_s(f,P,l,kT) ! pN*nm
  REAL(SP), INTENT(IN) :: f, P, l, kT
  REAL(SP) :: int_xdfWLC_s
  int_xdfWLC_s = f*xWLC(f,P,l,kT) - int_fdxWLC(f,P,l,kT)
END FUNCTION int_xdfWLC_s

FUNCTION int_xdfWLC_v(f,P,l,kT) ! pN*nm
  REAL(SP), DIMENSION(:), INTENT(IN) :: f
  REAL(SP), INTENT(IN) :: P, l, kT
  REAL(SP), DIMENSION(size(f)) :: int_xdfWLC_v
  int_xdfWLC_v = f*xWLC(f,P,l,kT) - int_fdxWLC(f,P,l,kT)
END FUNCTION int_xdfWLC_v

FUNCTION int_xdfWLC_Bouchiat(f,P,l,kT) ! pN*nm
  REAL(SP), INTENT(IN) :: f, P, l, kT
  REAL(SP) :: int_xdfWLC_Bouchiat
  int_xdfWLC_Bouchiat = f*xWLC_Bouchiat(f,P,l,kT) - int_fdxWLC_Bouchiat(f,P,l,kT)
END FUNCTION int_xdfWLC_Bouchiat

FUNCTION int_fdxWLC_s(f,P,l,kT) ! pN*nm
  REAL(SP), INTENT(IN) :: f, p, l, kT
  REAL(SP) :: x
  REAL(SP) :: int_fdxWLC_s
  x = xWLC(f,P,l,kT)
  if (l == 0.) then
    int_fdxWLC_s = 0._sp
  else 
    int_fdxWLC_s = 0.25*l/(1-x/l)-0.25*l-0.25*x+0.5*x*x/l
    int_fdxWLC_s = kT*int_fdxWLC_s/P
  end if
END FUNCTION int_fdxWLC_s

FUNCTION int_fdxWLC_v(f,P,l,kT) ! pN*nm
  REAL(SP), DIMENSION(:), INTENT(IN) :: f
  REAL(SP), INTENT(IN) :: P, l, kT
  REAL(SP), DIMENSION(size(f)) :: x
  REAL(SP), DIMENSION(size(f)) :: int_fdxWLC_v
  x = xWLC(f,P,l,kT)
  if (l == 0.) then
    int_fdxWLC_v = 0._sp
  else 
    int_fdxWLC_v = 0.25*l/(1-x/l)-0.25*l-0.25*x+0.5*x*x/l
    int_fdxWLC_v = kT*int_fdxWLC_v/P
  end if
END FUNCTION int_fdxWLC_v

FUNCTION int_fdxWLC_Bouchiat(f,P,l,kT)
  REAL(SP), INTENT(IN) :: f
  REAL(SP), INTENT(IN) :: P, l, kT
  REAL(SP) :: x
  REAL(SP) :: int_fdxWLC_Bouchiat
  REAL(SP) :: a2 = -0.5164228, a3 = -2.737418, a4 = 16.07497
  REAL(SP) :: a5 = -38.87607, a6 = 39.49944, a7 = -14.17718
  x = xWLC_Bouchiat(f,P,l,kT)
  if (l == 0.) then
    int_fdxWLC_Bouchiat = 0._sp
  else 
    int_fdxWLC_Bouchiat = 0.25*l/(1-x/l)-0.25*l-0.25*x+0.5*x*x/l
    int_fdxWLC_Bouchiat = int_fdxWLC_Bouchiat + a2*x**3/(3.*l**2) + a3*x**4/(4.*l**3) + a4*x**5/(5.*l**4) + &
		    a5*x**6/(6.*l**5) + a6*x**7/(7.*l**6) + a7*x**8/(8.*l**7)
    int_fdxWLC_Bouchiat = kT*int_fdxWLC_Bouchiat/P
  end if
END FUNCTION int_fdxWLC_Bouchiat

FUNCTION int_fdxWLCext_Bouchiat(f,P,l,S,kT)
  REAL(SP), INTENT(IN) :: f
  REAL(SP), INTENT(IN) :: P, l, S, kT
  REAL(SP) :: int_fdxWLCext_Bouchiat
  REAL(SP) :: a2 = -0.5164228, a3 = -2.737418, a4 = 16.07497
  REAL(SP) :: a5 = -38.87607, a6 = 39.49944, a7 = -14.17718
  REAL(SP) :: xf, dx = 1e-3, df = 1e-3
  REAL(SP), DIMENSION(:), ALLOCATABLE :: x, fint
  INTEGER(I4B) :: i,n 
  xf = xWLCext_Bouchiat(f,P,l,S,kT)

!   x = dx
!   int_fdxWLCext_Bouchiat = 0._sp
!   do while (x<=xf)
!     fint = fWLCext_Bouchiat(x,f,P,l,S,kT)
!     write(222,*), x, fint ! AIXO ESTA BE
!     int_fdxWLCext_Bouchiat = int_fdxWLCext_Bouchiat + f*dx
!     x = x+dx
!   end do
!   int_fdxWLCext_Bouchiat = 0._sp
!   do while (fint<=f)
!     fint = fint+df
!     x = xWLCext_Bouchiat(fint,P,l,S,kT)
!     int_fdxWLCext_Bouchiat = int_fdxWLCext_Bouchiat  + x*df
!     write(222,*) x, fint, int_fdxWLCext_Bouchiat  fint AIXO ESTA BE
!   end do
!   int_fdxWLCext_Bouchiat = -int_fdxWLCext_Bouchiat + f*xf
! print *, int_fdxWLCext_Bouchiat,f, 'AQUI'

  n = ceiling(xf/dx)
  allocate(fint(n), x(n))
  int_fdxWLCext_Bouchiat = 0._sp
  do i = 1, n
    x(i) = i*dx
    fint(i) = fWLCext_Bouchiat(x(i),f+df,P,l,S,kT)
!     write(222,*), x(i), fint(i)
  end do

  int_fdxWLCext_Bouchiat = (sum(fint(5:n-4)) + 55.*fint(2)/24.-fint(3)/6.+11.*fint(4)/8.+&
	  11.*fint(n-3)/8.-fint(n-2)/6.+55.*fint(n-1)/24.)*dx
!   print *, fint(n), f
  deallocate(fint, x)
END FUNCTION int_fdxWLCext_Bouchiat

FUNCTION xWLCext_Bouchiat(f,P,l,S,kT)
  REAL(SP), INTENT(IN) :: f
  REAL(SP), INTENT(IN) :: P, l, S, kT
  REAL(SP) :: xWLCext_Bouchiat
  xWLCext_Bouchiat = rtsafe_xWLCext_Bouchiat(f,P,l,S,kT,funcd_xWLCext_Bouchiat,0.+f/S,1.+f/S,1e-5)
  xWLCext_Bouchiat = xWLCext_Bouchiat*l
END FUNCTION xWLCext_Bouchiat

FUNCTION fWLCext_Bouchiat(x,fmax,P,l,S,kT)
  REAL(SP), INTENT(IN) :: x,fmax
  REAL(SP), INTENT(IN) :: P, l, S, kT
  REAL(SP) :: fWLCext_Bouchiat
  REAL(SP) :: a2 = -0.5164228, a3 = -2.737418, a4 = 16.07497
  REAL(SP) :: a5 = -38.87607, a6 = 39.49944, a7 = -14.17718
  fWLCext_Bouchiat = rtsafe_fWLCext_Bouchiat(x,P,l,S,kT,funcd_fWLCext_Bouchiat,0.,fmax,1.e-5)
END FUNCTION fWLCext_Bouchiat

SUBROUTINE WLC(z,f,P,kT,fval,fderiv)
! Proporciona el valor de la funcio 4*x**3+(4*P*f/kT)*x**2-1 i la seva derivada, avaluats a x
! Seveix per la subrutina xWLC
  REAL(SP), INTENT(IN) :: P, kT
  REAL(SP), INTENT(IN) :: z,f
  REAL(SP), INTENT(OUT) :: fval, fderiv
  fval = 4.*z**3+(4.*P*f/kT-3.)*z**2-1.
  fderiv = 12*z**2+2*(4*P*f/kT-3)*z;
END SUBROUTINE WLC

! Freely Jointed Chain, x=l*(1+f/Y)*(1/tanh(b*f/kT)-kT/(b*f))

FUNCTION xFJC_s(f,Y,b,l,kT) ! nm
  REAL(SP), INTENT(IN) :: f
  REAL(SP), INTENT(IN) :: Y, b, l, kT
  REAL(SP) :: xFJC_s
! equació x=l*(1+f/Y)*(1/tanh(b*f/kT)-kT/(b*f))
  if (l == 0) then
    xFJC_s = 0.
  else if ( Y >= 812000 ) then
    xFJC_s = l*(1./tanh(b*f/kT)-kT/(b*f))
  else
    xFJC_s = l*(1.+f/Y)*(1./tanh(b*f/kT)-kT/(b*f))
  end if
END FUNCTION xFJC_s

FUNCTION xFJC_v(f,Y,b,l,kT) ! nm
  REAL(SP), DIMENSION(:), INTENT(IN) :: f
  REAL(SP), INTENT(IN) :: Y, b, l, kT
  REAL(SP), DIMENSION(size(f)) :: xFJC_v
! equació x=l*(1+f/Y)*(1/tanh(b*f/kT)-kT/(b*f))
  if (l == 0.) then
    xFJC_v = 0.
  else if ( Y >= 812000 ) then
    xFJC_v = l*(1./tanh(b*f/kT)-kT/(b*f))
  else
    xFJC_v = l*(1.+f/Y)*(1./tanh(b*f/kT)-kT/(b*f))
  end if
END FUNCTION xFJC_v

FUNCTION int_xdfFJC_s(f,Y,b,l,kT)
  REAL(SP), INTENT(IN) :: f
  REAL(SP), INTENT(IN) :: Y, b, l, kT
  REAL(SP) :: int_xdfFJC_s
  INTEGER(I4B) :: j
  if ( Y < 812000.) then
    integral: do j = 1,8
      call midsquFJC(xFJC_v,0.,f,int_xdfFJC_s,j,Y,b,l,kT)
    end do integral
  else if (Y >= 812000.) then
    int_xdfFJC_s = l*kT*log(kT*sinh(f*b/kT)/(f*b))/b ! = \int xdf (sense mòdul de Young)
  end if
END FUNCTION int_xdfFJC_s

FUNCTION int_xdfFJC_v(f,Y,b,l,kT)
  REAL(SP), DIMENSION(:), INTENT(IN) :: f
  REAL(SP), INTENT(IN) :: Y, b, l, kT
  REAL(SP), DIMENSION(size(f)) :: int_xdfFJC_v
  INTEGER(I4B) :: j, i
  do i = 1, size(f)
    int_xdfFJC_v(i) = int_xdfFJC_s(f(i),Y,b,l,kT)
  end do
END FUNCTION int_xdfFJC_v

FUNCTION int_fdxFJC_s(f,Y,b,l,kT)
  REAL(SP), INTENT(IN) :: f, Y, b, l, kT
  REAL(SP) :: int_fdxFJC_s
  if ( (b == 0.).or.(l == 0.) ) then
    int_fdxFJC_s = 0.
  else
    int_fdxFJC_s = f*xFJC(f,Y,b,l,kT)-int_xdfFJC(f,Y,b,l,kT)
  end if
END FUNCTION int_fdxFJC_s

FUNCTION int_fdxFJC_v(f,Y,b,l,kT)
  REAL(SP), DIMENSION(:), INTENT(IN) :: f
  REAL(SP), INTENT(IN) :: Y, b, l, kT
  REAL(SP), DIMENSION(size(f)) :: int_fdxFJC_v
  int_fdxFJC_v = f*xFJC(f,Y,b,l,kT)-int_xdfFJC(f,Y,b,l,kT)
END FUNCTION int_fdxFJC_v

FUNCTION dfdxFJC(f,Y,b,l,kT)
  REAL(SP), INTENT(IN) :: f
  REAL(SP), INTENT(IN) :: Y, b, l, kT
  REAL(SP) :: dfdxFJC
!   dfdxFJC = l*( kT/(b*f**2) - b/(kT*sinh(b*f/kT)*sinh(b*f/kT)))
  dfdxFJC = (l*b/kT)*( (kT/(f*b))**2 - 1/(sinh(b*f/kT)*sinh(b*f/kT)) )
  dfdxFJC = 1._sp/dfdxFJC
END FUNCTION dfdxFJC

! Thick Chain (TC); A. Rosa, Micheletti

FUNCTION fTC(x,K,D,a,kT)
  REAL(SP), INTENT(IN) :: x, K, D, a, kT
  REAL(SP) :: part1, part2
  REAL(SP) :: fTC
! Dona la força donada una x/l0 per model de TC. Necessita la funció yTC(K,D/a)
  if ((x > 1.).or.(x < 0.)) call nrerror('fTC: x ha d''estar normalitzada per long. contorn')
  part1 = 2.*K*(sqrt(1.+1./(2.*K*(1.-x))*2)-sqrt(1.+1./(2.*K)**2))
  part2 = (3.*(1.-yTC(K,D/a))/(1.+yTC(K,D/a))-1./(2.*K*sqrt(1.+(1./2.*K)**2)))*x
  fTC = (kT/a)*(part1+part2)
END FUNCTION fTC

FUNCTION xTC(f,K,D,a,kT)
  REAL(SP), INTENT(IN) :: f, K, D, a, kT
  REAL(SP) :: xTC
! funció inversa de fTC
  xTC = rtsafeTC(TC,K,D,a,kT,f,0.,0.99,1E-5)
!   xTC = l*xTC
END FUNCTION xTC

FUNCTION yTC(x1,x2)
  REAL(SP), INTENT(IN) :: x1, x2
  REAL(SP) :: yTC, z
! Necessària per avaluar f(x) per a TC model. Conté els efectes de volum exclós
  if (x2 > 0.5) then
    z = 0.5*x1/(x2*x2)
    yTC = 0.5/(x2*x2)
    yTC = yTC*(1./(1-exp(z))+1/z)
    yTC = 1. - yTC
  else if (x2 <= 0.5) then
    yTC = 1./tanh(x1)-1./x1
  end if
END FUNCTION yTC

SUBROUTINE TC(z,f,K,D,a,kT,fval,fderiv)
  REAL(SP), INTENT(IN) :: z,f,K,D,a,kT
  REAL(SP), INTENT(OUT) :: fval, fderiv
  REAL(SP) :: c1, c2, c3
! Ajuda per trobar x donada f
  fval = fTC(z,K,D,a,kT)-f
  c1 = (3.*(1.-yTC(K,D/a))/(1.+yTC(K,D/a))-1./(2.*K*sqrt(1.+(1./2.*K)**2)))*kT/a
  c2 = 2*K*kT/a
  c3 = 0.5/K
  fderiv = c1+c3*c3*c2/((1-z)**3*sqrt(1+(c3/(1-z))**2))
END SUBROUTINE TC

! Adaptacions NR

SUBROUTINE midsquFJC(funk,aa,bb,s,n,Y,kuhn,l,kT)
!   USE nrtype; USE nrutil, ONLY : arth
!   IMPLICIT NONE
  REAL(SP), INTENT(IN) :: aa,bb
  REAL(SP), INTENT(IN) :: Y,kuhn,l,kT
  REAL(SP), INTENT(INOUT) :: s
  INTEGER(I4B), INTENT(IN) :: n
  INTERFACE
    FUNCTION funk(x,Y,kuhn,l,kT)
    USE nrtype
    REAL(SP), DIMENSION(:), INTENT(IN) :: x
    REAL(SP), INTENT(IN) :: Y, kuhn, l, kT
    REAL(SP), DIMENSION(size(x)) :: funk
    END FUNCTION funk
  END INTERFACE
  REAL(SP) :: a,b,del
  INTEGER(I4B) :: it
  REAL(SP), DIMENSION(2*3**(n-2)) :: x
  b=sqrt(bb-aa)
  a=0.0
  if (n == 1) then
    s=(b-a)*sum(func( (/0.5_sp*(a+b)/),Y,b,l,kT ))
  else
    it=3**(n-2)
    del=(b-a)/(3.0_sp*it)
    x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
    x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
    s=s/3.0_sp+del*sum(func(x,Y,b,l,kT))
  end if
CONTAINS
!BL
FUNCTION func(x,Y,kuhn,l,kT)
  REAL(SP), DIMENSION(:), INTENT(IN) :: x
  REAL(SP), INTENT(IN) :: Y,kuhn,l,kT
  REAL(SP), DIMENSION(size(x)) :: func
  func=2.0_sp*x*funk(bb-x**2,Y,kuhn,l,kT)
END FUNCTION func
END SUBROUTINE midsquFJC

FUNCTION rtsafeWLC(funcd,P,kT,force,x1,x2,xacc)
  REAL(SP), INTENT(IN) :: P, kT
  REAL(SP), INTENT(IN) :: force,x1,x2,xacc
  REAL(SP) :: rtsafeWLC
  INTERFACE
    SUBROUTINE funcd(x,force,P,kT,fval,fderiv)
      USE nrtype
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: P, kT
      REAL(SP), INTENT(IN) :: x,force
      REAL(SP), INTENT(OUT) :: fval,fderiv
    END SUBROUTINE funcd
  END INTERFACE
! Using the combination fo Newton-Raphson and bisection, find the root of a function bracketed between x1 and x2.
! The root, returned as the function value rtsafe, will be refined unit its accuracy is known within +/-xacc. f
! funcd is a user-supplied subroutien which returns both the function value and the first derivative of the function.
  INTEGER(I4B), PARAMETER :: MAXIT=100
  INTEGER(I4B) :: j
  REAL(SP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
  call funcd(x1,force,P,kT,fl,df);
  call funcd(x2,force,P,kT,fh,df);
  if ((fl > 0.0 .and. fh > 0.0) .or. &
    (fl < 0.0 .and. fh < 0.0)) &
    call nrerror('root must be bracketed in rtsafeWLC')
  if (fl == 0.0) then
    rtsafeWLC=x1
    RETURN
  else if (fh == 0.0) then
    rtsafeWLC=x2
    RETURN
  else if (fl < 0.0) then
    xl=x1
    xh=x2
  else
    xh=x1
    xl=x2
  end if
  rtsafeWLC=0.5_sp*(x1+x2)
  dxold=abs(x2-x1)
  dx=dxold
  call funcd(rtsafeWLC,force,P,kT,f,df)
  do j=1,MAXIT
    if (((rtsafeWLC-xh)*df-f)*((rtsafeWLC-xl)*df-f) > 0.0 .or. &
      abs(2.0_sp*f) > abs(dxold*df) ) then
      dxold=dx
      dx=0.5_sp*(xh-xl)
      rtsafeWLC=xl+dx
      if (xl == rtsafeWLC) RETURN
    else
      dxold=dx
      dx=f/df
      temp=rtsafeWLC
      rtsafeWLC=rtsafeWLC-dx
      if (temp == rtsafeWLC) RETURN
    end if
    if (abs(dx) < xacc) RETURN
    call funcd(rtsafeWLC,force,P,kT,f,df)
    if (f < 0.0) then
	    xl=rtsafeWLC
    else
	    xh=rtsafeWLC
	  end if
  end do
  call nrerror('rtsafeWLC: exceeded maximum iterations')
END FUNCTION rtsafeWLC

FUNCTION rtsafeTC(funcd,K,D,a,kT,force,x1,x2,xacc)
  REAL(SP), INTENT(IN) :: K, D, a, kT
  REAL(SP), INTENT(IN) :: force,x1,x2,xacc
  REAL(SP) :: rtsafeTC
  INTERFACE
    SUBROUTINE funcd(x,force,K,D,a,kT,fval,fderiv)
      USE nrtype
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: K, D, a, kT
      REAL(SP), INTENT(IN) :: x,force
      REAL(SP), INTENT(OUT) :: fval,fderiv
    END SUBROUTINE funcd
  END INTERFACE
! Using the combination fo Newton-Raphson and bisection, find the root of a function bracketed between x1 and x2.
! The root, returned as the function value rtsafeTC, will be refined unit its accuracy is known within +/-xacc. f
! funcd is a user-supplied subroutien which returns both the function value and the first derivative of the function.
  INTEGER(I4B), PARAMETER :: MAXIT=100
  INTEGER(I4B) :: j
  REAL(SP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
  call funcd(x1,force,K,D,a,kT,fl,df);
  call funcd(x2,force,K,D,a,kT,fh,df);
  if ((fl > 0.0 .and. fh > 0.0) .or. &
    (fl < 0.0 .and. fh < 0.0)) &
    call nrerror('root must be bracketed in rtsafeTC')
  if (fl == 0.0) then
    rtsafeTC=x1
    RETURN
  else if (fh == 0.0) then
    rtsafeTC=x2
    RETURN
  else if (fl < 0.0) then
    xl=x1
    xh=x2
  else
    xh=x1
    xl=x2
  end if
  rtsafeTC=0.5_sp*(x1+x2)
  dxold=abs(x2-x1)
  dx=dxold
  call funcd(rtsafeTC,force,K,D,a,kT,f,df)
  do j=1,MAXIT
    if (((rtsafeTC-xh)*df-f)*((rtsafeTC-xl)*df-f) > 0.0 .or. &
      abs(2.0_sp*f) > abs(dxold*df) ) then
      dxold=dx
      dx=0.5_sp*(xh-xl)
      rtsafeTC=xl+dx
      if (xl == rtsafeTC) RETURN
    else
      dxold=dx
      dx=f/df
      temp=rtsafeTC
      rtsafeTC=rtsafeTC-dx
      if (temp == rtsafeTC) RETURN
    end if
    if (abs(dx) < xacc) RETURN
    call funcd(rtsafeTC,force,K,D,a,kT,f,df)
    if (f < 0.0) then
	    xl=rtsafeTC
    else
	    xh=rtsafeTC
	  end if
  end do
  call nrerror('rtsafeTC: exceeded maximum iterations')
END FUNCTION rtsafeTC


FUNCTION rtsafeWLC_Bouchiat(force,P,l,kT,funcd,x1,x2,xacc)
  REAL(SP), INTENT(IN) :: P, kT,l
  REAL(SP), INTENT(IN) :: force,x1,x2,xacc
  REAL(SP) :: rtsafeWLC_Bouchiat
  INTERFACE
    SUBROUTINE funcd(x,force,P,l,kT,fval,fderiv)
      USE nrtype
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: P, kT, l
      REAL(SP), INTENT(IN) :: x,force
      REAL(SP), INTENT(OUT) :: fval,fderiv
    END SUBROUTINE funcd
  END INTERFACE
  INTEGER(I4B), PARAMETER :: MAXIT=100
  INTEGER(I4B) :: j
  REAL(SP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
  call funcd(x1,force,P,l,kT,fl,df);
  call funcd(x2,force,P,l,kT,fh,df);
  if ((fl > 0.0 .and. fh > 0.0) .or. (fl < 0.0 .and. fh < 0.0)) &
    call nrerror('root must be bracketed in rtsafeWLC_Bouchiat')
  if (fl == 0.0) then
    rtsafeWLC_Bouchiat=x1
    RETURN
  else if (fh == 0.0) then
    rtsafeWLC_Bouchiat=x2
    RETURN
  else if (fl < 0.0) then
    xl=x1
    xh=x2
  else
    xh=x1
    xl=x2
  end if
  rtsafeWLC_Bouchiat=0.5_sp*(x1+x2)
  dxold=abs(x2-x1)
  dx=dxold
  call funcd(rtsafeWLC_Bouchiat,force,P,l,kT,f,df)
  do j=1,MAXIT
    if (((rtsafeWLC_Bouchiat-xh)*df-f)*((rtsafeWLC_Bouchiat-xl)*df-f) > 0.0 .or. &
      abs(2.0_sp*f) > abs(dxold*df) ) then
      dxold=dx
      dx=0.5_sp*(xh-xl)
      rtsafeWLC_Bouchiat=xl+dx
      if (xl == rtsafeWLC_Bouchiat) RETURN
    else
      dxold=dx
      dx=f/df
      temp=rtsafeWLC_Bouchiat
      rtsafeWLC_Bouchiat=rtsafeWLC_Bouchiat-dx
      if (temp == rtsafeWLC_Bouchiat) RETURN
    end if
    if (abs(dx) < xacc) RETURN
    call funcd(rtsafeWLC_Bouchiat,force,P,l,kT,f,df)
    if (f < 0.0) then
	    xl=rtsafeWLC_Bouchiat
    else
	    xh=rtsafeWLC_Bouchiat
	  end if
  end do
  call nrerror('rtsafeWLC_Bouchiat: exceeded maximum iterations')
END FUNCTION rtsafeWLC_Bouchiat

SUBROUTINE funcd_WLC_Bouchiat(x,f,P,l,kT,fval,fderiv)
  REAL(SP), INTENT(IN) :: f, P, l, kT, x
  REAL(SP), INTENT(OUT) :: fval, fderiv
  REAL(SP) :: y
  REAL(SP) :: a2 = -0.5164228, a3 = -2.737418, a4 = 16.07497
  REAL(SP) :: a5 = -38.87607, a6 = 39.49944, a7 = -14.17718
  y = x
  fval = 1+((1-y)**2)*(-1+4*y+4*(a2*y**2+a3*y**3+a4*y**4+a5*y**5+a6*y**6+a7*y**7)-4*P*f/kT)
  fderiv = -2*(1-y)*(-1+4*y+4*(a2*y**2+a3*y**3+a4*y**4+a5*y**5+a6*y**6+a7*y**7)-4*P*f/kT) &
+ 4*((1-y)**2)*(1+2*a2*y+3*a3*y**2+4*a4*y**3+5*a5*y**4+6*a6*y**5+7*a7*y**6)
END SUBROUTINE funcd_WLC_Bouchiat

FUNCTION rtsafe_xWLCext_Bouchiat(force,P,l,S,kT,funcd,x1,x2,xacc)
  REAL(SP), INTENT(IN) :: P, kT,l,S
  REAL(SP), INTENT(IN) :: force,x1,x2,xacc
  REAL(SP) :: rtsafe_xWLCext_Bouchiat
  INTERFACE
    SUBROUTINE funcd(x,force,P,l,S,kT,fval,fderiv)
      USE nrtype
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: P, kT, l, S
      REAL(SP), INTENT(IN) :: x,force
      REAL(SP), INTENT(OUT) :: fval,fderiv
    END SUBROUTINE funcd
  END INTERFACE
  INTEGER(I4B), PARAMETER :: MAXIT=100
  INTEGER(I4B) :: j
  REAL(SP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
  call funcd(x1,force,P,l,S,kT,fl,df);
  call funcd(x2,force,P,l,S,kT,fh,df);
  if ((fl > 0.0 .and. fh > 0.0) .or. (fl < 0.0 .and. fh < 0.0)) &
    call nrerror('root must be bracketed in rtsafeWLC_Bouchiat')
  if (fl == 0.0) then
    rtsafe_xWLCext_Bouchiat=x1
    RETURN
  else if (fh == 0.0) then
    rtsafe_xWLCext_Bouchiat=x2
    RETURN
  else if (fl < 0.0) then
    xl=x1
    xh=x2
  else
    xh=x1
    xl=x2
  end if
  rtsafe_xWLCext_Bouchiat=0.5_sp*(x1+x2)
  dxold=abs(x2-x1)
  dx=dxold
  call funcd(rtsafe_xWLCext_Bouchiat,force,P,l,S,kT,f,df)
  do j=1,MAXIT
    if (((rtsafe_xWLCext_Bouchiat-xh)*df-f)*((rtsafe_xWLCext_Bouchiat-xl)*df-f) > 0.0 .or. &
      abs(2.0_sp*f) > abs(dxold*df) ) then
      dxold=dx
      dx=0.5_sp*(xh-xl)
      rtsafe_xWLCext_Bouchiat=xl+dx
      if (xl == rtsafe_xWLCext_Bouchiat) RETURN
    else
      dxold=dx
      dx=f/df
      temp=rtsafe_xWLCext_Bouchiat
      rtsafe_xWLCext_Bouchiat=rtsafe_xWLCext_Bouchiat-dx
      if (temp == rtsafe_xWLCext_Bouchiat) RETURN
    end if
    if (abs(dx) < xacc) RETURN
    call funcd(rtsafe_xWLCext_Bouchiat,force,P,l,S,kT,f,df)
    if (f < 0.0) then
	    xl=rtsafe_xWLCext_Bouchiat
    else
	    xh=rtsafe_xWLCext_Bouchiat
	  end if
  end do
  call nrerror('rtsafe_xWLCext_Bouchiat: exceeded maximum iterations')
END FUNCTION rtsafe_xWLCext_Bouchiat

SUBROUTINE funcd_xWLCext_Bouchiat(y,f,P,l,S,kT,fval,fderiv)
  REAL(SP), INTENT(IN) :: y, f
  REAL(SP), INTENT(IN) :: P, l, kT, S
  REAL(SP), INTENT(OUT) :: fval, fderiv
  REAL(SP) :: z
  REAL(SP) :: a2 = -0.5164228, a3 = -2.737418, a4 = 16.07497
  REAL(SP) :: a5 = -38.87607, a6 = 39.49944, a7 = -14.17718
  z = y-f/S
  fval = 1-4*((1-z)**2)*(f*P/kT+0.25-z-(a2*z**2+a3*z**3+a4*z**4+a5*z**5+a6*z**6+a7*z**7))
  fderiv = 8*(1-z)*(f*P/kT+0.25-z-(a2*z**2+a3*z**3+a4*z**4+a5*z**5+a6*z**6+a7*z**7)) &
      - 4*((1-z)**2)*(-1-(2*a2*z+3*a3*z**2+4*a4*z**3+5*a5*z**4+6*a6*z**5+7*a7*z**6))
END SUBROUTINE funcd_xWLCext_Bouchiat

FUNCTION rtsafe_fWLCext_Bouchiat(x,P,l,S,kT,funcd,force1,force2,facc)
  REAL(SP), INTENT(IN) :: P, kT,l,S
  REAL(SP), INTENT(IN) :: x,force1,force2,facc
  REAL(SP) :: rtsafe_fWLCext_Bouchiat
  INTERFACE
    SUBROUTINE funcd(force,x,P,l,S,kT,fval,fderiv)
      USE nrtype
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: P, kT, l, S
      REAL(SP), INTENT(IN) :: x,force
      REAL(SP), INTENT(OUT) :: fval,fderiv
    END SUBROUTINE funcd
  END INTERFACE
  INTEGER(I4B), PARAMETER :: MAXIT=100
  INTEGER(I4B) :: j
  REAL(SP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
  call funcd(force1,x,P,l,S,kT,fl,df);
  call funcd(force2,x,P,l,S,kT,fh,df);
  if ((fl > 0.0 .and. fh > 0.0) .or. (fl < 0.0 .and. fh < 0.0)) &
    call nrerror('root must be bracketed in rtsafeWLC_Bouchiat')
  if (fl == 0.0) then
    rtsafe_fWLCext_Bouchiat=force1
    RETURN
  else if (fh == 0.0) then
    rtsafe_fWLCext_Bouchiat=force2
    RETURN
  else if (fl < 0.0) then
    xl=force1
    xh=force2
  else
    xh=force1
    xl=force2
  end if
  rtsafe_fWLCext_Bouchiat=0.5_sp*(force1+force2)
  dxold=abs(force2-force1)
  dx=dxold
  call funcd(rtsafe_fWLCext_Bouchiat,x,P,l,S,kT,f,df)
  do j=1,MAXIT
    if (((rtsafe_fWLCext_Bouchiat-xh)*df-f)*((rtsafe_fWLCext_Bouchiat-xl)*df-f) > 0.0 .or. &
      abs(2.0_sp*f) > abs(dxold*df) ) then
      dxold=dx
      dx=0.5_sp*(xh-xl)
      rtsafe_fWLCext_Bouchiat=xl+dx
      if (xl == rtsafe_fWLCext_Bouchiat) RETURN
    else
      dxold=dx
      dx=f/df
      temp=rtsafe_fWLCext_Bouchiat
      rtsafe_fWLCext_Bouchiat=rtsafe_fWLCext_Bouchiat-dx
      if (temp == rtsafe_fWLCext_Bouchiat) RETURN
    end if
    if (abs(dx) < facc) RETURN
    call funcd(rtsafe_fWLCext_Bouchiat,x,P,l,S,kT,f,df)
    if (f < 0.0) then
	    xl=rtsafe_fWLCext_Bouchiat
    else
	    xh=rtsafe_fWLCext_Bouchiat
	  end if
  end do
  call nrerror('rtsafe_fWLCext_Bouchiat: exceeded maximum iterations')
END FUNCTION rtsafe_fWLCext_Bouchiat

SUBROUTINE funcd_fWLCext_Bouchiat(f,x,P,l,S,kT,fval,fderiv)
  REAL(SP), INTENT(IN) :: f, x
  REAL(SP), INTENT(IN) :: P, l, S, kT
  REAL(SP), INTENT(OUT) :: fval, fderiv
  REAL(SP) :: z, y
  REAL(SP) :: a2 = -0.5164228, a3 = -2.737418, a4 = 16.07497
  REAL(SP) :: a5 = -38.87607, a6 = 39.49944, a7 = -14.17718
  y = x/l
  z = y-f/S
  fval = 1-4*((1-z)**2)*(f*P/kT+0.25-z-(a2*z**2+a3*z**3+a4*z**4+a5*z**5+a6*z**6+a7*z**7))  
  fderiv = -8*(1-z)*(f*P/kT+0.25-z-(a2*z**2+a3*z**3+a4*z**4+a5*z**5+a6*z**6+a7*z**7))/S & 
      - 4*(1-z)*(P/kT + 1./S + (2*a2*z+3*a3*z**2+4*a4*z**3+5*a5*z**4+6*a6*z**5+7*a7*z**6)/S)
END SUBROUTINE funcd_fWLCext_Bouchiat


END MODULE aaElasticModel
