! Align the unfolding trajectories between themselves and the folding 
! ones likewise. User selects two regions F0a-F0b and F1a-F1b used for 
! the alignment. All the points in each region are selected, and the
! linear fits are compared to the first. Output a list of shift factors.
! Multiple datasets can be aligned together by setting the index integer 
! jj in the inputfile. 

! INPUT FILE = arxiu amb: adreça, f0a, f0b, f1a, f1b, index+
! OUTPUT FILE = arxiu anomenat SHIFT.txt amb: trajectoria, N0, sh0, N1, shi1, 2

! última modificació: 18/08/10 (defal i modul, per posar align default)
! última modificació: 24/08/10 (proba per alinear millor quan no és 2)
! última modificació: 25/09/12 (identificació adequada dels intervals on fer el fit lineal)

MODULE Parametres
  USE nrtype
  IMPLICIT NONE
  INTEGER(I4B), PARAMETER :: columns = 4
  INTEGER(I4B), PARAMETER :: Nmin = 6 ! minimum number of points to draw a fit 
  INTEGER(I4B), PARAMETER :: defal = 0 ! 0,1,2, valors possibles
  REAL(SP), PARAMETER :: dshMAX = 10.!7.
END MODULE Parametres

PROGRAM dumbAlign
  USE nrtype; USE nrutil, ONLY : nrerror, assert_eq
  USE nr, ONLY : fit
  USE amFileUtil, ONLY : CountRows
  USE Parametres
  IMPLICIT NONE
  REAL(SP) :: F0a, F0b, F1a, F1b, F0m, F1m, y0m, y1m
  REAL(SP) :: sh0, sh1, dsh, sh0old, sh1old
  REAL(SP) :: a, b, siga, sigb, chi2, q
  REAL(SP), DIMENSION(:), ALLOCATABLE :: yv, Fv
  REAL(SP), DIMENSION(:), ALLOCATABLE :: y0Uref, y1Uref, y0Fref, y1Fref
  INTEGER(I4B) :: N0, N1, sz, j, jj, i, control = 0
  INTEGER(I4B) :: counter0, counter1, dim
  INTEGER(I4B) :: OpenStatus, ReadStatus
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: idx
  CHARACTER(1) :: flag
  CHARACTER(256) :: inputfile, molecula, traj
  
  print *, "Please enter the path of the input file"
  print *, '(no spaces allowed, enclose between quotes "" if slashes / are required):'
  read *, inputfile
  open(12, file=trim(inputfile), status='old', iostat=OpenStatus)
  if (OpenStatus /= 0) call nrerror("ERR001: Inputfile not found.")
  sz = CountRows(12,'#')
  allocate(idx(sz),y0Uref(sz),y0Fref(sz),y1Uref(sz),y1Fref(sz))
  idx = 0; y0Uref = 0.0_sp; y0Fref = 0.0_sp; y1Uref = 0.0_sp; y1Fref = 0.0_sp
  j = 0
  outerloop: do 
    !!! Continugut inputfile !!!
    read(12,*,iostat=ReadStatus) molecula, F0a, F0b, F1a, F1b, jj
    !!! !!!!!!!!!!!!!!!!!!!! !!!
    if (ReadStatus < 0) exit outerloop
    if (ReadStatus > 0) call nrerror("ERR002: Error reading the Inputfile.")
    if (molecula(:1)=='#') cycle 
    j = j+1; print *, 'Alineant '//trim(molecula)
    F0m = 0.5_sp*(F0a+F0b); F1m = 0.5_sp*(F1a+F1b)
    idx(j) = jj
    open(17, file=trim(molecula)//"/SHIFT.txt", status='replace')
    open(16, file=trim(molecula)//"/LOG.txt", status='old', iostat=OpenStatus)
    if (OpenStatus /= 0) call nrerror("ERR003: LOGfile not found.")
    if ( j == 1 .or. all(idx(:j-1) /= jj) ) then
      print *, 'New set of alignments!'
      y0Uref(j) = 0.0_sp; y1Uref(j) = 0.0_sp
      y0Fref(j) = 0.0_sp; y1Fref(j) = 0.0_sp
    else 
      print *, 'Old alignment'
      do i = 1, j-1
	if ( idx(i) == jj ) then
	  y0Uref(j) = y0Uref(i); y0Fref(j) = y0Fref(i)
	  y1Uref(j) = y1Uref(i); y1Fref(j) = y1Fref(i)
	  exit
	end if
      end do
    end if
    innerloop: do
      read(16, *, iostat=ReadStatus) traj, counter0, counter1
      if ( ReadStatus < 0 ) exit innerloop
      if ( ReadStatus > 0 ) call nrerror("ERR004: Error reading the LOGfile.")
      if ( traj(1:1) == "#" ) cycle
      flag = traj(1:1)           ! u/f
      if ( flag /= 'f' .and. flag /= 'u' ) call nrerror("ERR009: Unrecognized flag.")
      open(14, file=trim(molecula)//"/events_raw/"//trim(traj)//".txt",status='old',iostat=OpenStatus)
      if ( OpenStatus /= 0 ) call nrerror("ERR005: Datafile not found.")
      dim = CountRows(14) ! counter1-counter0
      allocate ( yv(dim), Fv(dim) )
      call DataFill(flag,14,F0a,F0b,F1a,F1b,yv,Fv,N0,N1)
      close(14)
      ! fa els fits lineals
      if ( N0 < Nmin .or. N1 < Nmin ) then
	print *,trim(traj),N0,N1
	call nrerror("ERR012: N0 or N1 are too small.")
      end if
! Fit a la zona baixa
!     sh0 = 0; sh1 = 0
    if (N0 >= Nmin) then
      call fit(yv(1:N0),Fv(1:N0),a,b,siga,sigb,chi2,q)
      y0m = (F0m-a)/b ! y that corresponds to the middle value F0m
      select case (flag)
      case ('u')
	if ( y0Uref(j) == 0 ) y0Uref(j) = y0m
	sh0 = y0m-y0Uref(j);
      case ('f')
	if ( y0Fref(j) == 0 ) y0Fref(j) =y0m
	sh0 = y0m-y0Fref(j);
      end select
    end if
! Fit a al zona alta de forces
    if (N1 >= Nmin) then
      call fit(yv(dim-N1+1:dim),Fv(dim-N1+1:dim),a,b,siga,sigb,chi2,q)
      y1m = (F1m-a)/b ! y that corresponds to the middle value F0m
      select case (flag)
      case ('u')
	if ( y1Uref(j) == 0 ) y1Uref(j) = y1m
	sh1 = y1m-y1Uref(j)
      case ('f')
	if ( y1Fref(j) == 0 ) y1Fref(j) =y1m
	sh1 = y1m-y1Fref(j)
      end select
    end if
      dsh = abs(sh0-sh1); 
      if (defal == 2) then
	if ( (dsh < dshMAX) ) then
	  write(17,*) trim(traj), N0, sh0, N1, sh1, defal, '"'//trim(molecula)//'"'
	else if (dsh >= dshMAX) then
	  if (traj(1:1) == 'u') write(17,*) trim(traj), N0, sh0, N1, sh1, 0, '"'//trim(molecula)//'"'
	  if (traj(1:1) == 'f') write(17,*) trim(traj), N0, sh0, N1, sh1, 1, '"'//trim(molecula)//'"'
! 	else if (dsh == 0.) then ! vol dir que N1 = 0
! 	  print *, trim(traj), sh0, sh1, dsh, 'n1'
! 	  write(17,*) trim(traj), N0, sh0, N1, sh1, 0, '"'//trim(molecula)//'"'
! 	else if (dsh == sh1) then ! vol dir que N0 = 0
! 	  print *, trim(traj), sh0, sh1, dsh, 'no'
! 	  write(17,*) trim(traj), N0, sh0, N1, sh1, 1, '"'//trim(molecula)//'"'
	end if
      else
	write(17,*) trim(traj), N0, sh0, N1, sh1, defal, '"'//trim(molecula)//'"'
      end if
      deallocate(yv,Fv)
    end do innerloop
    close(17)
    close(16)
  end do outerloop
  close(12)

CONTAINS

SUBROUTINE DataFill(flag,unit,F0a,F0b,F1a,F1b,yv,Fv,N0,N1)
  CHARACTER(1), INTENT(IN) :: flag
  INTEGER(I4B), INTENT(IN) :: unit
  REAL(SP), INTENT(IN) :: F0a, F0b, F1a, F1b
  REAL(SP), DIMENSION(:), INTENT(OUT) :: yv, Fv
  INTEGER(I4B), INTENT(OUT) :: N0, N1
  ! Fill the vectors yv and Fv with data read from unit. The flag
  ! is 'u' or 'f' for unfolding or folding trajectory. F0a and F0b
  ! are the limits of the first region, F1a and F1b are the limits 
  ! of the second one. On output, N0 and N1 are the number of points 
  ! that have been selected.
  REAL(SP) :: yt, t, fy
  REAL(SP), DIMENSION(dim) :: f, y
  INTEGER(I4B) :: dim, ReadStatus, cnt, i
  CHARACTER(1) :: option
  CHARACTER(2) :: region
  !initialization phase
  N0 = 0; N1 = 0; yv = 0.0_sp; Fv = 0.0_sp
  dim = assert_eq(size(yv),size(Fv),'DataFill')
  readfile: do i = 1, dim
!     call ReadLine(unit,option,ReadStatus,cnt,yt,t,fy)
    read(unit,*,iostat=ReadStatus) cnt, yt, t, fy
    if (ReadStatus /= 0) call nrerror("ERR008: Error reading a datafile.")
    f(i) = fy; y(i) = yt
  end do readfile

  select case (flag)
  case ('u')
    i = 0
    readloopu00: do
      i = i+1; if (i > dim) call nrerror('Problem 00, '//trim(traj))
      if (f(i) >= F0a) then
	N0 = 1
	yv(N0) = y(i); Fv(N0) = f(i)
	exit readloopu00
      end if
    end do readloopu00
    readloopu01: do
      i = i+1; if (i > dim) call nrerror('Problem 01, '//trim(traj))
      N0 = N0+1
      yv(N0) = y(i); Fv(N0) = f(i)
      if (f(i) > F0b) exit readloopu01
    end do readloopu01
    i = dim+1
    readloopu10: do
      i = i-1; if (i < N0) call nrerror('Problem 10')
      if (f(i) < F1b) then
	N1 = 1
	yv(dim) = y(i); Fv(dim) = f(i)
	exit readloopu10
      end if
    end do readloopu10
    readloopu11: do
      i = i-1; if (i < N0) call nrerror('Problem 11')
      N1 = N1+1
      yv(dim-N1+1) = y(i); Fv(dim-N1+1) = f(i)
      if (f(i) < F1a) exit readloopu11
    end do readloopu11

  case ('f')
    i = dim+1
    readloopf00: do
      i = i-1; if (i <= 0) call nrerror('Problem 00f')
      if (f(i) >= F0a) then
	N0 = 1
	yv(N0) = y(i); Fv(N0) = f(i)
	exit readloopf00
      end if
    end do readloopf00
    readloopf01: do
      i = i-1; if (i <=0 ) call nrerror('Problem 01f, '//trim(traj))
      N0 = N0+1
      yv(N0) = y(i); Fv(N0) = f(i)
      if (f(i) > F0b) exit readloopf01
    end do readloopf01
    i = 0
    readloopf10: do
      i = i+1; if (i > dim-N0) call nrerror('Problem 10f')
      if (f(i) < F1b) then
	N1 = 1
	yv(dim) = y(i); Fv(dim) = f(i)
	exit readloopf10
      end if
    end do readloopf10
    readloopf11: do
      i = i+1; if (i > dim-N0) call nrerror('Problem 10f')
      N1 = N1+1
      yv(dim-N1+1) = y(i); Fv(dim-N1+1) = f(i)
      if (f(i) < F1a) exit readloopf11
    end do readloopf11

  case default
    call nrerror('ERR010: flag must be either "u" or "f"')
  end select

!   select case (flag)
!   case ('u')
!     region = 'n0'
!   case ( 'f' ) 
!     region = 'n2'
!   case default
!     call nrerror('ERR010: flag must be either "u" or "f"')
!   end select
!   readloop: do
! !        call ReadLine(unit,option,ReadStatus,cnt,yt,t,fy)
!     if (columns == 4 ) then
!       read(unit,*,iostat=ReadStatus) cnt, yt, t, fy
!     else if (columns == 3) then
!       read(unit,*,iostat=ReadStatus) cnt, yt, fy
!     else
!       call nrerror('Columnes a determinar')
!     end if
!     if ( ReadStatus < 0 ) exit readloop
!     if ( ReadStatus > 0 ) call nrerror("ERR008: Error reading a datafile.")
!     select case(region)
!     case ('n0') ! unfolding, inside low region
! !       if (( fy > F0a ).and.(yt > 180.)) then
!       if ( fy > F0a ) then
! 	N0 = 1
! 	yv(N0) = yt
! 	Fv(N0) = fy
! 	region = 'w0'
!       end if
!     case ('w0')
!       if ( flag == 'u' .and. fy > F0b ) region = 'n1'
!       if ( flag == 'f' .and. fy < F0a ) exit readloop
!       N0 = N0+1
!       yv(N0) = yt
!       Fv(N0) = fy
!     case ('n1')
!       if ( flag == 'u' .and. fy > F1a ) then
! 	N1 = 1
! 	yv(dim) = yt
! 	Fv(dim) = fy
! 	region = 'w1'
!       else if ( flag == 'f' .and. fy < F0b ) then
! 	N0 = 1
! 	yv(1) = yt
! 	Fv(1) = fy
! 	region = 'w0'
!       end if
!     case ('w1')
!       if ( flag == 'u' .and. fy > F1b ) exit readloop
!       if ( flag == 'f' .and. fy < F1a ) region = 'n1'
!       N1 = N1+1
!       yv(dim-N1+1) = yt
!       Fv(dim-N1+1) = fy
!     case ('n2')
!       if ( fy < F1b ) then
! 	N1 = 1
! 	yv(dim) = yt
! 	Fv(dim) = fy
! 	region = 'w1'
!       end if
!     case default
!       call nrerror("ERR011: unrecognized region code in DataFill.")
!     end select
!   end do readloop

END SUBROUTINE DataFill

!   SUBROUTINE ReadLine(unit,option,ReadStatus,cnt,yt,t,fy)
!     INTEGER(I4B), INTENT(IN) :: unit
!     CHARACTER(1), INTENT(IN) :: option
!     INTEGER(I4B), INTENT(OUT) :: ReadStatus, cnt
!     REAL(SP), INTENT(OUT) :: yt, t, fy
!     ! This routine is needed for the program to work equally well
!     ! with data in which the molecular extension is known or not.
!     select case(option)
!     case('3')
!        t = 0.0_sp
!        read(unit,*,iostat=ReadStatus) cnt, yt, fy
!     case('4')
!        read(unit,*,iostat=ReadStatus) cnt, yt, t, fy
!     end select
!   END SUBROUTINE ReadLine

END PROGRAM dumbAlign

!------------------------------------!
! From Numerical Recipes             !
! Routines for performing linear fit !
!------------------------------------!

  SUBROUTINE fit(x,y,a,b,siga,sigb,chi2,q,sig)
    USE nrtype; USE nrutil, ONLY : assert_eq
    USE nr, ONLY : gammq
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
    REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
    REAL(SP), DIMENSION(:), OPTIONAL, INTENT(IN) :: sig
    INTEGER(I4B) :: ndata
    REAL(SP) :: sigdat,ss,sx,sxoss,sy,st2
    REAL(SP), DIMENSION(size(x)), TARGET :: t
    REAL(SP), DIMENSION(:), POINTER :: wt
    if (present(sig)) then
       ndata=assert_eq(size(x),size(y),size(sig),'fit')
       wt=>t
       wt(:)=1.0_sp/(sig(:)**2)
       ss=sum(wt(:))
       sx=dot_product(wt,x)
       sy=dot_product(wt,y)
    else
       ndata=assert_eq(size(x),size(y),'fit')
       ss=real(size(x),sp)
       sx=sum(x)
       sy=sum(y)
    end if
    sxoss=sx/ss
    t(:)=x(:)-sxoss
    if (present(sig)) then
       t(:)=t(:)/sig(:)
       b=dot_product(t/sig,y)
    else
       b=dot_product(t,y)
    end if
    st2=dot_product(t,t)
    b=b/st2
    a=(sy-sx*b)/ss
    siga=sqrt((1.0_sp+sx*sx/(ss*st2))/ss)
    sigb=sqrt(1.0_sp/st2)
    t(:)=y(:)-a-b*x(:)
    q=1.0
    if (present(sig)) then
       t(:)=t(:)/sig(:)
       chi2=dot_product(t,t)
       if (ndata > 2) q=gammq(0.5_sp*(size(x)-2),0.5_sp*chi2)
    else
       chi2=dot_product(t,t)
       sigdat=sqrt(chi2/(size(x)-2))
       siga=siga*sigdat
       sigb=sigb*sigdat
    end if
  END SUBROUTINE fit

  FUNCTION gammq_s(a,x)
    USE nrtype; USE nrutil, ONLY : assert
    USE nr, ONLY : gcf,gser
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: a,x
    REAL(SP) :: gammq_s
    call assert( x >= 0.0,  a > 0.0, 'gammq_s args')
    if (x<a+1.0_sp) then
       gammq_s=1.0_sp-gser(a,x)
    else
       gammq_s=gcf(a,x)
    end if
  END FUNCTION gammq_s

  FUNCTION gammq_v(a,x)
    USE nrtype; USE nrutil, ONLY : assert,assert_eq
    USE nr, ONLY : gcf,gser
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
    REAL(SP), DIMENSION(size(a)) :: gammq_v
    LOGICAL(LGT), DIMENSION(size(x)) :: mask
    INTEGER(I4B) :: ndum
    ndum=assert_eq(size(a),size(x),'gammq_v')
    call assert( all(x >= 0.0),  all(a > 0.0), 'gammq_v args')
    mask = (x<a+1.0_sp)
    gammq_v=merge(1.0_sp-gser(a,merge(x,0.0_sp,mask)), &
         gcf(a,merge(x,0.0_sp,.not. mask)),mask)
  END FUNCTION gammq_v

  FUNCTION gcf_s(a,x,gln)
    USE nrtype; USE nrutil, ONLY : nrerror
    USE nr, ONLY : gammln
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: a,x
    REAL(SP), OPTIONAL, INTENT(OUT) :: gln
    REAL(SP) :: gcf_s
    INTEGER(I4B), PARAMETER :: ITMAX=100
    REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
    INTEGER(I4B) :: i
    REAL(SP) :: an,b,c,d,del,h
    if (x == 0.0) then
       gcf_s=1.0
       RETURN
    end if
    b=x+1.0_sp-a
    c=1.0_sp/FPMIN
    d=1.0_sp/b
    h=d
    do i=1,ITMAX
       an=-i*(i-a)
       b=b+2.0_sp
       d=an*d+b
       if (abs(d) < FPMIN) d=FPMIN
       c=b+an/c
       if (abs(c) < FPMIN) c=FPMIN
       d=1.0_sp/d
       del=d*c
       h=h*del
       if (abs(del-1.0_sp) <= EPS) exit
    end do
    if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_s')
    if (present(gln)) then
       gln=gammln(a)
       gcf_s=exp(-x+a*log(x)-gln)*h
    else
       gcf_s=exp(-x+a*log(x)-gammln(a))*h
    end if
  END FUNCTION gcf_s

  FUNCTION gcf_v(a,x,gln)
    USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
    USE nr, ONLY : gammln
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
    REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
    REAL(SP), DIMENSION(size(a)) :: gcf_v
    INTEGER(I4B), PARAMETER :: ITMAX=100
    REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
    INTEGER(I4B) :: i
    REAL(SP), DIMENSION(size(a)) :: an,b,c,d,del,h
    LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
    i=assert_eq(size(a),size(x),'gcf_v')
    zero=(x == 0.0)
    where (zero)
       gcf_v=1.0
    elsewhere
       b=x+1.0_sp-a
       c=1.0_sp/FPMIN
       d=1.0_sp/b
       h=d
    end where
    converged=zero
    do i=1,ITMAX
       where (.not. converged)
          an=-i*(i-a)
          b=b+2.0_sp
          d=an*d+b
          d=merge(FPMIN,d, abs(d)<FPMIN )
          c=b+an/c
          c=merge(FPMIN,c, abs(c)<FPMIN )
          d=1.0_sp/d
          del=d*c
          h=h*del
          converged = (abs(del-1.0_sp)<=EPS)
       end where
       if (all(converged)) exit
    end do
    if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_v')
    if (present(gln)) then
       if (size(gln) < size(a)) call &
            nrerror('gser: Not enough space for gln')
       gln=gammln(a)
       where (.not. zero) gcf_v=exp(-x+a*log(x)-gln)*h
    else
       where (.not. zero) gcf_v=exp(-x+a*log(x)-gammln(a))*h
    end if
  END FUNCTION gcf_v

  FUNCTION gser_s(a,x,gln)
    USE nrtype; USE nrutil, ONLY : nrerror
    USE nr, ONLY : gammln
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: a,x
    REAL(SP), OPTIONAL, INTENT(OUT) :: gln
    REAL(SP) :: gser_s
    INTEGER(I4B), PARAMETER :: ITMAX=100
    REAL(SP), PARAMETER :: EPS=epsilon(x)
    INTEGER(I4B) :: n
    REAL(SP) :: ap,del,summ
    if (x == 0.0) then
       gser_s=0.0
       RETURN
    end if
    ap=a
    summ=1.0_sp/a
    del=summ
    do n=1,ITMAX
       ap=ap+1.0_sp
       del=del*x/ap
       summ=summ+del
       if (abs(del) < abs(summ)*EPS) exit
    end do
    if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_s')
    if (present(gln)) then
       gln=gammln(a)
       gser_s=summ*exp(-x+a*log(x)-gln)
    else
       gser_s=summ*exp(-x+a*log(x)-gammln(a))
    end if
  END FUNCTION gser_s

  FUNCTION gser_v(a,x,gln)
    USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
    USE nr, ONLY : gammln
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
    REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
    REAL(SP), DIMENSION(size(a)) :: gser_v
    INTEGER(I4B), PARAMETER :: ITMAX=100
    REAL(SP), PARAMETER :: EPS=epsilon(x)
    INTEGER(I4B) :: n
    REAL(SP), DIMENSION(size(a)) :: ap,del,summ
    LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
    n=assert_eq(size(a),size(x),'gser_v')
    zero=(x == 0.0)
    where (zero) gser_v=0.0
    ap=a
    summ=1.0_sp/a
    del=summ
    converged=zero
    do n=1,ITMAX
       where (.not. converged)
          ap=ap+1.0_sp
          del=del*x/ap
          summ=summ+del
          converged = (abs(del) < abs(summ)*EPS)
       end where
       if (all(converged)) exit
    end do
    if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_v')
    if (present(gln)) then
       if (size(gln) < size(a)) call &
            nrerror('gser: Not enough space for gln')
       gln=gammln(a)
       where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gln)
    else
       where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gammln(a))
    end if
  END FUNCTION gser_v

  FUNCTION gammln_s(xx)
    USE nrtype; USE nrutil, ONLY : arth,assert
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: xx
    REAL(SP) :: gammln_s
    REAL(DP) :: tmp,x
    REAL(DP) :: stp = 2.5066282746310005_dp
    REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
         -86.50532032941677_dp,24.01409824083091_dp,&
         -1.231739572450155_dp,0.1208650973866179e-2_dp,&
         -0.5395239384953e-5_dp/)
    call assert(xx > 0.0, 'gammln_s arg')
    x=xx
    tmp=x+5.5_dp
    tmp=(x+0.5_dp)*log(tmp)-tmp
    gammln_s=tmp+log(stp*(1.000000000190015_dp+&
         sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
  END FUNCTION gammln_s

  FUNCTION gammln_v(xx)
    USE nrtype; USE nrutil, ONLY: assert
    IMPLICIT NONE
    INTEGER(I4B) :: i
    REAL(SP), DIMENSION(:), INTENT(IN) :: xx
    REAL(SP), DIMENSION(size(xx)) :: gammln_v
    REAL(DP), DIMENSION(size(xx)) :: ser,tmp,x,y
    REAL(DP) :: stp = 2.5066282746310005_dp
    REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
         -86.50532032941677_dp,24.01409824083091_dp,&
         -1.231739572450155_dp,0.1208650973866179e-2_dp,&
         -0.5395239384953e-5_dp/)
    if (size(xx) == 0) RETURN
    call assert(all(xx > 0.0), 'gammln_v arg')
    x=xx
    tmp=x+5.5_dp
    tmp=(x+0.5_dp)*log(tmp)-tmp
    ser=1.000000000190015_dp
    y=x
    do i=1,size(coef)
       y=y+1.0_dp
       ser=ser+coef(i)/y
    end do
    gammln_v=tmp+log(stp*ser/x)
  END FUNCTION gammln_v
