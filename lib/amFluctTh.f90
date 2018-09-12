! Routines for performing tests of Fluctuation Theorems
! Written by Alessandro Mossa.

MODULE amFluctTh
  USE nrtype; USE nrutil, ONLY : nrerror
  IMPLICIT NONE

  INTERFACE zdiff
     MODULE PROCEDURE zdiff_s, zdiff_v
  END INTERFACE

  INTERFACE average
     MODULE PROCEDURE average_s, average_d
  END INTERFACE

CONTAINS

  FUNCTION Ben(WU,WR,logP)
    REAL(DP), DIMENSION(:), INTENT(IN) :: WU, WR
    REAL(DP), INTENT(IN), OPTIONAL :: logP
    REAL(DP), DIMENSION(2) :: Ben
    ! Given the Unfolding work distribution WU and the 
    ! Refolding work distribution WR, outputs the estimate
    ! of the free energy according to Bennett's method.
    ! The second component of the answer is the associated
    ! error, estimated according to Shirts, Bair, Hooker and
    ! Pande, Phys. Rev. Lett. 91 (2003) 140601.
! L'error només té sentit quan les distribucions de treball U i R estan superposades.
    REAL(SP) ::xR, xU, x1, x2
    REAL(SP), PARAMETER :: tol = 1.0E-4
    xR = -sum(WR)/real(size(WR),sp)
    xU =  sum(WU)/real(size(WU),sp)
    x1 = min(xU,xR)-10.0_sp
    x2 = max(xU,xR)+10.0_sp
    if ( present(logP) ) then
       Ben(1) = zbrent(WU,WR,x1,x2,tol,logP)
    else
       Ben(1) = zbrent(WU,WR,x1,x2,tol)
    end if
    Ben(2) = BenError(WU,WR,Ben(1))
  END FUNCTION Ben

  FUNCTION BenError(WU,WR,DF)
    REAL(DP), DIMENSION(:), INTENT(IN) :: WU, WR
    REAL(DP), INTENT(IN) :: DF
    REAL(DP) :: BenError
    ! Compute the statistical error of the Bennett estimator
    ! according to Shirts, Bair, Hooker and Pande, 
    ! Phys. Rev. Lett. 91 (2003) 140601.
    INTEGER(I4B) :: nU, nR, nT
    REAL(DP) :: M
    REAL(DP), DIMENSION(size(WU)+size(WR)) :: aux
    nU = size(WU); nR = size(WR); nT = nU + nR
    M = log(real(nU,dp)/real(nR,dp))
    aux(1:nU)    = M-DF+WU(1:nU)
    aux(nU+1:nT) = M-DF+WR(1:nR)
    aux = 2.0_dp+2.0_dp*cosh(aux)
    aux = 1.0_dp/aux
    BenError = sum(aux)/real(nT,dp)
    BenError = 1.0_dp/BenError-(real(nT,dp)/real(nR,dp)+real(nT,dp)/real(nU,dp))
    BenError = sqrt(BenError/real(nT,dp))
  END FUNCTION BenError

  FUNCTION func(x,W,nU,nR)
    REAL(DP), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: W
    INTEGER(I4B), INTENT(IN) :: nU, NR
    REAL(DP), DIMENSION(size(W)) :: func
    ! func = [1+r*exp(W-x)]^-1 where r = nU/nR
    REAL(DP) :: r
    r = real(nU,dp)/real(nR,dp)
    func = 1.0_dp+r*exp(W-x)
    func = 1.0_dp/func
  END FUNCTION func

  FUNCTION Jarz(W)
    REAL(DP), DIMENSION(:), INTENT(IN) :: W
    REAL(DP) :: Jarz
    ! Jarzynski estimator. It is understood that the values of W
    ! are already expressed in units of kBT.
    REAL(DP) :: W0
    REAL(DP), DIMENSION(size(W)) :: aux
    INTEGER(I4B) :: sz
    sz = size(W)
    W0 = sum(W)/real(sz,dp)
    aux = W-W0
    aux = dexp(-aux)
    Jarz = W0-dlog(sum(aux)/real(sz,dp))
  END FUNCTION Jarz

FUNCTION JarzLargeNumbers(W)
  REAL(DP), DIMENSION(:), INTENT(IN) :: W
  REAL(DP) :: JarzLargeNumbers
! Jarzynski estimator. 
! We need to use it in case work values are too large (over 1000kT)
! The exponential is evaluated in a tricky way
  REAL(DP) :: W0, logW
  REAL(DP), DIMENSION(size(W)) :: Waux, aux
  REAL(DP), DIMENSION(size(W)) :: A
  INTEGER(I4B), DIMENSION(size(W)) :: B, Baux
  INTEGER(I4B) :: B0, Bmax
  INTEGER(I4B) :: sz, i
  INTEGER(I4B), DIMENSION(1) :: imax

  sz = size(W)
  W0 = sum(W)/real(sz,dp)
  Waux = W!-W0
  B = -int(Waux*dlog10(dexp(1._dp)),sp)
  A = 10._dp**(-Waux*dlog10(dexp(1._dp))-real(B,dp))
  B0 = int(sum(real(B,sp))/sz, sp)
  Baux = B-B0;
  if ( maxval(Baux)>308) then
    i = maxval(Baux)-308
    B0 = B0+i
    Baux = B-B0
  end if
  imax = maxloc(Baux)
  Bmax = maxval(B)+100
!   do i = 1, sz
!     print *, Waux(i), B(i), B(i)-Bmax, B(i)-B0, A(i)*10._dp**B(i), A(i)*10._dp**(B(i)-Bmax), &
! 	    A(i)*10._dp**(B(i)-B0)
!   end do
!   print *, B0, Bmax, W0, maxval(W), minval(W)

! Mètode 1  
!   if ( Baux(imax(1))<=308) then
!     logW = B0*dlog(10._dp)+dlog(sum(A*10._dp**Baux)/real(sz,dp))
!     JarzLargeNumbers = -logW !+ W0
!   else 
! !     logW = dlog(A(imax(1))/real(sz,dp)) + B(imax(1))*dlog(10._dp)
!     print *, 'Aproximació a treball mínim'
!     JarzLargeNumbers = minval(W)
!   end if

! Mètode 2
!   aux = dexp(-Waux)
!   JarzLargeNumbers = -dlog(sum(aux)/real(sz,dp))

! Mètode 3
  aux = A*10._dp**(Baux)
  JarzLargeNumbers = -B0*dlog(10._dp)-dlog(sum(aux)/real(sz,dp))

END FUNCTION JarzLargeNumbers

  SUBROUTINE logplot(Wfor,Wrev,lp,unitF,unitR)
    USE nrutil, ONLY : arth, nrerror
    USE amProbDistrib, ONLY : histogram
    REAL(SP), DIMENSION(:), INTENT(IN) :: Wfor, Wrev
    REAL(SP), DIMENSION(:,:), INTENT(OUT) :: lp
    INTEGER(I4B), INTENT(IN), OPTIONAL :: unitF, unitR
    ! Given a forward work distribution Wfor and a reverse work distribution Wrev,
    ! and a fixed Nbins, build a histogram and find the linear fit of the logplot
    ! used for the verification of the fluctuation theorem. NOTE: it is understood 
    ! that Wrev is in fact the distribution of -W along the reverse trajectory.
    REAL(SP) :: Wmin, Wmax, maxRev, minRev, maxFor, minFor, binwidth
    REAL(SP), DIMENSION(:,:), ALLOCATABLE :: histoFor, histoRev
    INTEGER(I4B) :: NbinsFor, NbinsRev, iFor, iRev, Nbins, Ncol, j
    Nbins = size(lp,1); Ncol = size(lp,2)
    if ( Nbins < 2 ) call nrerror("logplot: lp must have at least 2 rows") 
    if ( Ncol < 3 ) call nrerror("logplot: lp must have at least 3 columns")
    maxRev = maxval(Wrev); minRev = minval(Wrev)
    maxFor = maxval(Wfor); minFor = minval(Wfor)
    Wmin = max(minFor,minRev)
    Wmax = min(maxFor,maxRev)
    binwidth = (Wmax-Wmin)/real(Nbins-1,sp) 
    iFor = nint((minFor-Wmin)/binwidth)
    iRev = nint((minRev-Wmin)/binwidth)
    NbinsFor = Nbins-iFor+nint((maxFor-Wmax)/binwidth)
    NbinsRev = Nbins-iRev+nint((maxRev-Wmax)/binwidth)
    allocate(histoFor(NbinsFor,3),histoRev(NbinsRev,3))
    histoFor(:,1) = arth(Wmin+iFor*binwidth,binwidth,NbinsFor)
    histoRev(:,1) = arth(Wmin+iRev*binwidth,binwidth,NbinsRev)
    histoFor(:,2:3) = 0.0
    histoRev(:,2:3) = 0.0
    call histogram(Wfor,histoFor)
    call histogram(Wrev,histoRev)
    lp(:,1) = arth(Wmin,binwidth,Nbins)
    do j = 1, Nbins
       lp(j,2) = log(histoFor(j-iFor,2)/histoRev(j-iRev,2))
       lp(j,3) = sqrt((histoFor(j-iFor,3)/histoFor(j-iFor,2))**2+(histoRev(j-iRev,3)/histoRev(j-iRev,2))**2)
    end do
    if ( present(unitF) ) then
       do j = 1, NbinsFor
          write(unitF,*) histoFor(j,1:3)
       end do
    end if
    if ( present(unitR) ) then
       do j = 1, NbinsRev
          write(unitR,*) histoRev(j,1:3)
       end do
    end if
    deallocate(histoFor,histoRev)
  END SUBROUTINE logplot

  FUNCTION zbrent(WU,WR,x1,x2,tol,logP)
    REAL(DP), DIMENSION(:), INTENT(IN) :: WU, WR
    REAL(SP), INTENT(IN) :: x1,x2,tol
    REAL(DP), INTENT(IN), OPTIONAL :: logP
    REAL(SP) :: zbrent
    INTEGER(I4B), PARAMETER :: ITMAX=100
    REAL(SP), PARAMETER :: EPS=epsilon(x1)
    ! Using Brent's method, find the root of the function zdiff 
    ! known to lie between x1 and x2. The root, returned as 
    ! zbrent, will be refined until its accuracy is tol.
    ! Parameters: Maximum allowed number of iterations, and 
    ! machine floating-point precision.
    ! Adapted from Numerical Recipes.
    INTEGER(I4B) :: iter
    REAL(SP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    a=x1
    b=x2
    if ( present(logP) ) then
       fa=zdiff(real(a,dp),WU,WR,logP)-a
       fb=zdiff(real(b,dp),WU,WR,logP)-b
    else
       fa=zdiff(real(a,dp),WU,WR)-a
       fb=zdiff(real(b,dp),WU,WR)-b
    end if
! print *,fa,fb,'fa,fb del zbrent per trobar l''arrel...'
    if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
	 print *, fa, fb, '-> error a zbrent'
!          call nrerror('root must be bracketed for zbrent')
	fb = -fa
    end if
    c=b
    fc=fb
    do iter=1,ITMAX
       if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
          c=a
          fc=fa
          d=b-a
          e=d
       end if
       if (abs(fc) < abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       end if
       tol1=2.0_sp*EPS*abs(b)+0.5_sp*tol
       xm=0.5_sp*(c-b)
       if (abs(xm) <= tol1 .or. fb == 0.0) then
          zbrent=b
          RETURN
       end if
       if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s=fb/fa
          if (a == c) then
             p=2.0_sp*xm*s
             q=1.0_sp-s
          else
             q=fa/fc
             r=fb/fc
             p=s*(2.0_sp*xm*q*(q-r)-(b-a)*(r-1.0_sp))
             q=(q-1.0_sp)*(r-1.0_sp)*(s-1.0_sp)
          end if
          if (p > 0.0) q=-q
          p=abs(p)
          if (2.0_sp*p  <  min(3.0_sp*xm*q-abs(tol1*q),abs(e*q))) then
             e=d
             d=p/q
          else
             d=xm
             e=d
          end if
       else
          d=xm
          e=d
       end if
       a=b
       fa=fb
       b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
       if ( present(logP) ) then
          fb=zdiff(real(b,dp),WU,WR,logP)-b
       else
          fb=zdiff(real(b,dp),WU,WR)-b
       end if
    end do
    call nrerror('zbrent: exceeded maximum iterations')
    zbrent=b
  END FUNCTION zbrent

  FUNCTION zdiff_s(x,WU,WR,logP)
    REAL(DP), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: WU, WR
    REAL(DP), INTENT(IN), OPTIONAL :: logP 
    REAL(SP) :: zdiff_s
    ! Function zR-zU, scalar version.
    INTEGER(I4B) :: nU, nR
    nU = size(WU)
    nR = size(WR)
! print *, WR
    if ( present(logP) ) then
       zdiff_s = real(zR(x,WR,nU,nR)-zU(x,WU,nU,nR)-logP,sp)
    else
       zdiff_s = real(zR(x,WR,nU,nR)-zU(x,WU,nU,nR),sp)
    end if
!     print *,zdiff_s,'zR(x)-zF(x) escalar (zdiff)'
  END FUNCTION zdiff_s

  FUNCTION zdiff_v(x,WU,WR,logP)
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: WU, WR
    REAL(DP), INTENT(IN), OPTIONAL :: logP
    REAL(SP), DIMENSION(size(x)) :: zdiff_v
    ! difference zR-zU. The optional argument is the logarithm of 
    ! the probability P(A->B)/P(A)
    INTEGER(I4B) :: sz, i, nU, nR
    nU = size(WU)
    nR = size(WR)
! print *, WR
    sz = size(x)
    if ( present(logP) ) then
       do i = 1, sz
          zdiff_v(i) = real(zR(x(i),WR,nU,nR)-zU(x(i),WU,nU,nR)-logP,sp)
       end do
    else
       do i = 1, sz
          zdiff_v(i) = real(zR(x(i),WR,nU,nR)-zU(x(i),WU,nU,nR),sp)
       end do
    end if
  END FUNCTION zdiff_v

 FUNCTION zR(x,W,nU,nR)
   REAL(DP), INTENT(IN) :: x
   REAL(DP), DIMENSION(:), INTENT(IN) :: W
   INTEGER(I4B), INTENT(IN) :: nU, nR
   REAL(DP) :: zR
   ! Function zR, with tricks to prevent numerical overflowing
   REAL(DP) :: W0
   REAL(DP), DIMENSION(size(W)) :: aux
   INTEGER(I4B) :: sz
   sz = size(W)
   W0 = sum(W)/real(sz,dp)
   aux = W-W0
   zR = sum(func(x,-W,nU,nR)*exp(-aux))/real(sz,dp)
   zR = log(zR)-W0
 END FUNCTION zR

  FUNCTION zU(x,W,nU,nR)
    REAL(DP), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: W
    INTEGER(I4B), INTENT(IN) :: nU, nR
    REAL(DP) :: zU
    ! Function zU 
    REAL(DP), DIMENSION(size(W)) :: aux
    aux = func(x,W,nU,nR)
    zU = average(aux)
!     print *, zU, 'zU exponencial'
    zU = log(zU)
!     print *, zU, 'zU logaritme'
  END FUNCTION zU

  FUNCTION average_s(x)
    REAL(SP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP) :: average_s
    ! Smart way of averaging the vector of data x
    INTEGER(I4B) :: sz, i
    sz = size(x)
    if (sz < 1 ) call nrerror("average_s: wrong dimension of x")
    average_s = x(1)
    do i = 2, sz
       average_s = average_s*real(i-1,dp)/real(i,dp)+x(i)/real(i,dp)
    end do
  END FUNCTION average_s

  FUNCTION average_d(x)
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP) :: average_d
    ! Smart way of averaging the vector of data x
    INTEGER(I4B) :: sz, i  
    sz = size(x)
    if (sz < 1) call nrerror("average_d: wrong dimension of x")
    average_d = x(1)
    do i = 2, sz
       average_d = average_d*real(i-1,dp)/real(i,dp)+x(i)/real(i,dp)
    end do
  END FUNCTION average_d
	
END MODULE amFluctTh
