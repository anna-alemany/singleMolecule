! Utilities for analysing FDC
! Written by Anna Alemany

! creat: 17/07/12
! modificat: 28/06/13

MODULE aaFDCanalysis

  USE nrtype; USE nrutil, ONLY : nrerror, assert_eq, imaxloc,outerprod,swap, arth, assert,poly

  IMPLICIT NONE

!   INTERFACE int_fdxFJC
!     MODULE PROCEDURE int_fdxFJC_s, int_fdxFJC_v
!   END INTERFACE

  CONTAINS

! Find Rate

FUNCTION Rate(y,t)
! Fa un fit lineal de les dades i agafa el pendent com a velocitat
  REAL(SP), DIMENSION(:), INTENT(IN) :: y, t
  REAl(SP), DIMENSION(size(t)) :: t1
  REAL(SP) :: Rate
  REAL(SP) :: v_naif, t0
  REAL(SP) :: syx, sx, sy, sxx
  INTEGER(I4B) :: sz, i
  
  sz = assert_eq(size(y), size(t), 'aaFDCanalysis: Rate')
  if (sz < 3) call nrerror('aaFDCanalysis: too few data')
  t0 = t(1); t1 = t-t0

  syx = sum(y(1:sz)*t1(1:sz)); sxx = sum(t1(1:sz)*t1(1:sz))
  sx = sum(t1(1:sz)); sy = sum(y(1:sz))

  v_naif = (y(sz)-y(1))/(t(sz)-t(1))
  Rate = (sz*syx-sx*sy)/(sz*sxx-sx*sx) 
!   if (abs(v_naif - Rate)>3._sp) call nrerror('aaFDCanalysis: not linear')

END FUNCTION Rate

! find 1st rupture force, lambda and time; hopping number

SUBROUTINE RutpureForce1(tipus,s0,fv,yv,tv,a,b,c, fr,yr,tr,dfr, hop,di)
! Adaptació de ForcesRutpure-v2
  CHARACTER(1), INTENT(IN) :: tipus
  REAL(SP), DIMENSION(:), INTENT(IN) :: fv, yv, tv
  REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c
  INTEGER(I4B), INTENT(IN) :: s0
  REAL(SP), INTENT(OUT) :: fr, yr, tr,dfr
  INTEGER(I4B), INTENT(OUT) :: hop, di

  INTEGER(I4B), PARAMETER :: av = 5, dime = 10

  REAL(SP) :: f1, y1, t1, f2, y2, t2
  REAL(SP), DIMENSION(size(a)) :: f, dist
  REAL(SP), DIMENSION(av) :: fav, yav
  INTEGER(I4B) :: i, j, k
  INTEGER(I4B) :: sz, X, n
  INTEGER(I4B) :: fila, s1, i1
  INTEGER(I4B), DIMENSION(1) :: minDistInd

  X = assert_eq(size(a), size(b), size(c), 'aaFDCanalysis: RuptureForce1, dimension aU')
  sz = assert_eq(size(fv), size(yv), size(tv), 'aaFDCanalysis: RuptureForce1, dimension f,y,t')

!   select case(tipus)
!   case ('u')
!     s0 = 1
!   case ('f')
!     s0 = X
!   end select
  hop = 0

  n = 0; i = 0
  do while ( i <= sz-1 )
    i = i+1
    f = a+b*yv(i)+c*yv(i)*yv(i)
    dist = abs(f-fv(i)); minDistInd = minloc(dist); s1 = minDistInd(1)
    if ( s1 == s0 ) n = 1
    if ( s1 /= s0 ) n = n+1
    if ( n == 2) then
      if (i >= av) then
	f1 = sum(fav)/real(av,sp)
	y1 = sum(yav)/real(av,sp)
      else 
	f1 = sum(fav)/real(i-1,sp)
	y1 = sum(yav)/real(i-1,sp)
      end if
      t1 = tv(i)
      i1 = i        
    end if   
    if ( n == dime ) then
      fav(mod(i-1,av)+1) = fv(i)
      do j = dime,av
	i = i+1
	fav(mod(i-1,av)+1) = fv(i)
	yav(mod(i-1,av)+1) = yv(i)
      end do
      f2 = sum(fav)/real(av,sp)
      y2 = sum(yav)/real(av,sp)
      t2 = tv(i)
      di = i
      i = i1

      hop = hop+1; 
      if (hop == 1) then
	yr = 0.5*(y1+y2)
	fr = f1
	tr = 0.5*(t1+t2)
	dfr = f2-f1
	exit
      end if
!       s0 = s1
    end if

    fav(mod(i-1,av)+1) = fv(i) 
    yav(mod(i-1,av)+1) = yv(i)

  end do

!   if (hop == 0) call nrerror('aaFDCanalysis: RuptureForce1, hop=0')

END SUBROUTINE RutpureForce1

FUNCTION DetectState(fav, yav, a,b,c)
  REAL(SP), INTENT(IN) :: fav, yav
  REAL(SP), DIMENSION(:) :: a,b,c
! ! ! 
  REAL(SP), DIMENSION(1) :: minDistInd
  REAL(SP), DIMENSION(:), ALLOCATABLE :: fU, dist
  INTEGER(I4B) :: X
  INTEGER(I4B) :: DetectState
  X = assert_eq(size(a),size(b),size(c),'DetectState')
  allocate(fU(X), dist(X))
  fU = a+b*yav+c*yav*yav
  dist = abs(fav-fU); minDistInd = minloc(dist); DetectState = minDistInd(1)
  deallocate(fU, dist)
END FUNCTION DetectState

! find last rupture forces and hopping number

SUBROUTINE lastRutpureForce(tipus,fv,yv,a,b,c,av,dime,fmin,fmax,fr,yr, hop)
! Adaptació de ForcesRutpure-v2
  CHARACTER(1), INTENT(IN) :: tipus
  REAL(SP), INTENT(IN) :: fmin, fmax
  REAL(SP), DIMENSION(:), INTENT(IN) :: fv, yv
  REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c
  REAL(SP), INTENT(OUT) :: fr, yr
  INTEGER(I4B), INTENT(OUT) :: hop

  INTEGER(I4B), INTENT(IN) :: av, dime

  REAL(SP) :: f1, y1, t1, f2, y2, t2
  REAL(SP), DIMENSION(size(a)) :: f, dist
  REAL(SP), DIMENSION(av) :: fav, yav
  INTEGER(I4B) :: i, j, k
  INTEGER(I4B) :: sz, X, n
  INTEGER(I4B) :: fila, s0, s1, i1
  INTEGER(I4B), DIMENSION(1) :: minDistInd

  X = assert_eq(size(a), size(b), size(c), 'aaFDCanalysis: lastRuptureForce, dimension fit')
  sz = assert_eq(size(fv), size(yv), 'aaFDCanalysis: lastRuptureForce, dimension raw data')

  select case(tipus)
  case ('u')
    s0 = 1
  case ('f')
    s0 = X
  end select
  hop = 0

  n = 0; i = 0
  do while ( i <= sz-1 )
    i = i+1
    if ((fv(i)>fmax).or.(fv(i)<fmin)) cycle
    f = a+b*yv(i)+c*yv(i)*yv(i)
    dist = abs(f-fv(i)); minDistInd = minloc(dist); s1 = minDistInd(1)
    if ( s1 == s0 ) n = 1
    if ( s1 /= s0 ) n = n+1
    if ( n == 2) then
      f1 = sum(fav)/real(av,sp)
      y1 = sum(yav)/real(av,sp)
      i1 = i        
    end if   
    if ( n == dime ) then
      fav(mod(i-1,av)+1) = fv(i)
      do j = dime,av
	i = i+1
	fav(mod(i-1,av)+1) = fv(i)
	yav(mod(i-1,av)+1) = yv(i)
      end do
      f2 = sum(fav)/real(av,sp)
      y2 = sum(yav)/real(av,sp)
      i = i1

      hop = hop+1; 
!       if (hop == 1) then
      yr = 0.5*(y1+y2)
      fr = f1
!       end if
      s0 = s1
    end if

    fav(mod(i-1,av)+1) = fv(i) 
    yav(mod(i-1,av)+1) = yv(i)

  end do

  if (hop == 0) call nrerror('aaFDCanalysis: RuptureForce1, hop=0')

END SUBROUTINE lastRutpureForce

! evaluate work

FUNCTION Work(traj,fv,yv,y0,y1)
  CHARACTER(5), INTENT(IN) :: traj
  REAL(SP), DIMENSION(:), INTENT(IN) :: fv, yv
  REAL(SP), INTENT(IN) :: y0, y1
  REAL(SP) :: Work
! 
  INTEGER(I4B) :: sz, cnt, c0, c1, j
  REAL(SP) :: yt, ym, fy
  REAL(SP) :: f0, f1
  REAL(SP), DIMENSION(:), ALLOCATABLE :: vy, vf

  call cutData(y0,y1,yv,traj(1:1),c0,c1)
  sz = c1-c0+1
  allocate(vy(sz), vf(sz))
  vy=yv(c0:c1); vf=fv(c0:c1)

  select case (traj(1:1))
  case('u')
    call normalize(traj,y0,vy(1:2),vf(1:2),f0)
    vy(1) = y0; vf(1) = f0
    call normalize(traj,y1,vy(sz-1:sz),vf(sz-1:sz),f1)
    vy(sz) = y1; vf(sz) = f1
  case('f')
    call normalize(traj,y1,vy(1:2),vf(1:2),f1)
    vy(1) = y1; vf(1) = f1
    call normalize(traj,y0,vy(sz-1:sz),vf(sz-1:sz),f0)
    vy(sz) = y0; vf(sz) = f0
  case default
    call nrerror("aaFDC, work: unrecognized traj.")
  end select
  ! Computation of the integrals and output of the results
  Work = integral(vy,vf)

END FUNCTION Work

SUBROUTINE cutData(x0,x1,data,flag,c0,c1)
  REAL(SP), INTENT(IN) :: x0, x1
  REAL(SP), DIMENSION(:), INTENT(IN) :: data
  CHARACTER(1), INTENT(IN) :: flag
  INTEGER(I4B), INTENT(OUT) :: c0, c1
  ! given a vector of data, cut outside the values x0 and x1.
  ! the output are the indices where the cut is to be applied
  INTEGER(I4B) :: sz
  sz = size(data)
  select case(flag)
  case('u')
    do c0 = 1, sz
      if ( data(c0+1) > x0 ) exit
    end do
    do c1 = sz, 1, -1
      if ( data(c1-1) < x1 ) exit
    end do
  case('f')
    do c0 = 1, sz
      if ( data(c0+1) < x1 ) exit
    end do
    do c1 = sz, 1, -1
      if ( data(c1-1) > x0 ) exit
    end do
  case default
    call nrerror("aaFDC, cutData: flag must be either 'u' or 'f'")
  end select
END SUBROUTINE cutData

SUBROUTINE normalize(traj,y,yv,fv,f)
  CHARACTER(5), INTENT(IN) :: traj
  REAL(SP), INTENT(IN) :: y
  REAL(SP), DIMENSION(2), INTENT(IN) :: yv, fv!, lv
  REAL(SP), INTENT(OUT) :: f!, l
  ! Normalize the integration extrema
  if ( ( y < yv(1) .and. y < yv(2) ) .or. ( y > yv(1) .and. y > yv(2) ) ) then
    call nrerror("aaFDC, normalize: bad range in normalize")
  end if
  if ( yv(1) /= yv(2) ) then 
    f  = (fv(1)-fv(2))/(yv(1)-yv(2))*(y-yv(2))+fv(2)
  else
    f  = 0.5_sp*(fv(1)+fv(2))
  end if
END SUBROUTINE normalize

FUNCTION trap(x1,y1,x2,y2)
  REAL(SP), INTENT(IN) :: x1, y1, x2, y2
  REAL(DP) :: trap
  ! Computes the area of the trapezium defined by the x-axis and the points 
  ! (x1,y1) and (x2,y2).
  REAL(DP) :: xA, yA, xB, yB
  xA = real(x1,dp)
  xB = real(x2,dp)
  yA = real(y1,dp)
  yB = real(y2,dp)
  trap = 0.5_dp*(yA+yB)*(xB-xA)
END FUNCTION trap

FUNCTION integral(x,y)
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: x, y
  REAL(DP) :: integral
  ! Avaluació de l'area
  REAL(SP) :: m, q, y2, y1, x2, x1
  REAL(DP) :: int0, inttrap
  REAL(DP), DIMENSION(size(y)-1) :: z
  REAL(SP), DIMENSION(size(y)-1) :: w2
  REAL(SP), DIMENSION(size(y)) :: ajuda
  INTEGER(I4B) :: sz, i

  sz = assert_eq(size(x),size(y),'integral')
  x1 = x(1); x2 = x(sz)
  y1 = y(1); y2 = y(sz)
  m = (y2-y1)/(x2-x1)  ! pendent "global" de les dades
  q = y1-m*x1          ! terme independent "global"    
  int0 = real(q+m*(x(sz)+x(1))/real(2,sp),dp)
  int0 = int0*real(x(sz)-x(1),dp)  ! Integral de la recta, que és l'àrea del trapezi:
  do i = 1, sz-1
    z(i) = trap(x(i),y(i)-m*x(i)-q,x(i+1),y(i+1)-m*x(i+1)-q)
  end do
  integral = sum(z)
  integral = integral+int0
END FUNCTION integral

! SGsmooth

SUBROUTINE SGsmooth(nrr,nl,m,ld,x,y,yf)
  INTEGER(I4B), INTENT(IN) :: nrr, nl, m, ld	! Paràmetres de filtratge
  REAl(SP), DIMENSION(:), INTENT(IN) :: x, y	! dades a filtrar (suposo x equiespaiat)
  REAl(SP), DIMENSION(:), INTENT(OUT):: yf 	! dades filtrades

  REAL(SP), DIMENSION(nrr+nl+1) :: c ! coeficients de Savitzky-Golay

  REAL(SP) :: dx, xx, yy
  INTEGER(I4B) :: n, i, j
  INTEGER(I4B) :: nrrr, nll
  INTEGER(I4B) :: OpenStatus, ReadStatus
  CHARACTER(50) :: input, output 

  n = assert_eq(size(y), size(y), size(yf), 'SGsmooth: dimensions input/output vectors')
  c = savgol(nl,nrr,ld,m)

  dx = sum(x(2:n)-x(1:n-1))/real(n,sp)
  yf = y
  do i = 1+nl, n-nrr
    yf(i) = 0._sp
    do j = 1, nl+1
      yf(i) = yf(i)+c(j)*y(i-j+1)
    end do
    do j = 1, nrr
      yf(i) = yf(i)+c(nl+1+j)*y(i+j)
    end do
  end do

  do i = 1, nl
    nll = i-1
    nrrr = nrr+nl-nll
    c = savgol(nll,nrrr,ld,m)
    yf(i) = 0._sp
    do j = 1, nll+1
      yf(i) = yf(i)+c(j)*y(i-j+1)
    end do
    do j = 1, nrrr
      yf(i) = yf(i)+c(nll+1+j)*y(i+j)
    end do
  end do

  do i = n-nrr+1, n
    nrrr = n-i
    nll = nrr+nl-nrrr
    c = savgol(nll,nrrr,ld,m)
    yf(i) = 0._sp
    do j = 1, nll+1
      yf(i) = yf(i)+c(j)*y(i-j+1)
    end do
    do j = 1, nrrr
      yf(i) = yf(i)+c(nll+1+j)*y(i+j)
    end do
  end do
  yf(2) = 0.5*(yf(1)+yf(3)) 
! print *, 'A:',yf
! print *, 'B:',real(fac(ld),sp)*yf(i)/dx**ld
  yf = real(fac(ld),sp)*yf/dx**ld
  
END SUBROUTINE SGsmooth

FUNCTION fac(n)
  INTEGER(I4B), INTENT(IN) :: n
  INTEGER(I4B) :: fac
  INTEGER(I4B) :: i
  if (n < 0) then
    call nrerror('No es pot calcular el factorial')
  else if (n == 0) then
    fac = 1
  else if (n > 0) then
    fac = 1
    do i = 1, n
      fac = fac*i
    end do
  end if
END FUNCTION fac

! Numerical Recipes - adapatations

FUNCTION savgol(nl,nrr,ld,m)
  INTEGER(I4B), INTENT(IN) :: nl,nrr,ld,m
  REAL(SP), DIMENSION(nl+nrr+1) :: savgol
  INTEGER(I4B) :: imj,ipj,mm,np
  INTEGER(I4B), DIMENSION(m+1) :: indx
  REAL(SP) :: d,sm
  REAL(SP), DIMENSION(m+1) :: b
  REAL(SP), DIMENSION(m+1,m+1) :: a
  INTEGER(I4B) :: irng(nl+nrr+1)
  call assert(nl >= 0, nrr >= 0, ld <= m, nl+nrr >= m, 'savgol args')
  do ipj=0,2*m
    sm=sum(arth(1.0_sp,1.0_sp,nrr)**ipj)+&
	sum(arth(-1.0_sp,-1.0_sp,nl)**ipj)
    if (ipj == 0) sm=sm+1.0_sp
    mm=min(ipj,2*m-ipj)
    do imj=-mm,mm,2
      a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sm
    end do
  end do
  call ludcmp(a(:,:),indx(:),d)
  b(:)=0.0
  b(ld+1)=1.0
  call lubksb(a(:,:),indx(:),b(:))
  savgol(:)=0.0
  irng(:)=arth(-nl,1,nrr+nl+1)
  np=nl+nrr+1
  savgol(mod(np-irng(:),np)+1)=poly(real(irng(:),sp),b(:))
END FUNCTION savgol

SUBROUTINE lubksb(a,indx,b)
  REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
  INTEGER(I4B) :: i,n,ii,ll
  REAL(SP) :: summ
  n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
  ii=0
  do i=1,n
    ll=indx(i)
    summ=b(ll)
    b(ll)=b(i)
    if (ii /= 0) then
      summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
    else if (summ /= 0.0) then
      ii=i
    end if
    b(i)=summ
  end do
  do i=n,1,-1
    b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
  end do
END SUBROUTINE lubksb

SUBROUTINE ludcmp(a,indx,d)
  REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
  INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
  REAL(SP), INTENT(OUT) :: d
  REAL(SP), DIMENSION(size(a,1)) :: vv
  REAL(SP), PARAMETER :: TINY=1.0e-20_sp
  INTEGER(I4B) :: j,n,imax
  n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
  d=1.0
  vv=maxval(abs(a),dim=2)
  if (any(vv == 0.0)) call nrerror('singular matrix in ludcmp')
  vv=1.0_sp/vv
  do j=1,n
    imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
    if (j /= imax) then
      call swap(a(imax,:),a(j,:))
      d=-d
      vv(imax)=vv(j)
    end if
    indx(j)=imax
    if (a(j,j) == 0.0) a(j,j)=TINY
    a(j+1:n,j)=a(j+1:n,j)/a(j,j)
    a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
  end do
END SUBROUTINE ludcmp


END MODULE aaFDCanalysis
