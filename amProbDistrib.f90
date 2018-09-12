! Routines for estimating the probability density function of some data.
! Written by Alessandro Mossa.

MODULE amProbDistrib
  USE nrtype; USE nrutil, ONLY : arth, assert_eq, nrerror
  IMPLICIT NONE
  
  INTERFACE F0
     MODULE PROCEDURE F0_s, F0_v
  END INTERFACE F0
  
  INTERFACE PDF
     MODULE PROCEDURE PDF_s, PDF_v
  END INTERFACE

  INTERFACE Percentile
     MODULE PROCEDURE Percentile_s, Percentile_v
  END INTERFACE
  
CONTAINS

  FUNCTION F0_s(a,b,x)
    REAL(SP), INTENT(IN) :: a, b, x
    REAL(SP) :: F0_s
    ! Initial approximation to the EPD
    if ( x <= a ) then
       F0_s = 0.0_sp
    else if ( x >= b ) then
       F0_s = 1.0_sp
    else
       F0_s = (x-a)/(b-a)
    end if
  END FUNCTION F0_s

  FUNCTION F0_v(a,b,x)
    REAL(SP), INTENT(IN) :: a, b
    REAL(SP), DIMENSION(:), INTENT(IN) :: x
    REAL(SP), DIMENSION(size(x)) :: F0_v
    ! Initial approximation to the EPD
    where ( x <= a ) 
       F0_v = 0.0_sp
    elsewhere ( x >= b )
       F0_v = 1.0_sp
    elsewhere
       F0_v = (x-a)/(b-a)
    end where
  END FUNCTION F0_v
    
  FUNCTION Fcoef(ECDF,xval,j)
    REAL(SP), DIMENSION(:), INTENT(IN) :: ECDF, xval
    INTEGER(I4B), INTENT(IN) :: j
    REAL(SP) :: Fcoef
    ! j-th coefficient of the Fourier series of the remainder
    REAL(SP) :: a, b
    REAL(SP), DIMENSION(size(ECDF)) :: aux
    INTEGER(I4B) :: N
    N = assert_eq(size(ECDF),size(xval),'Fcoef')
    a = xval(1); b = xval(N)
    if ( mod(j,2) == 1 ) then
       Fcoef = -1.0_sp
    else
       Fcoef = 1.0_sp
    end if
    aux = cos(PI*j*(xval-a)/(b-a))
    Fcoef = Fcoef-dot_product(ECDF(1:N-1),aux(2:N)-aux(1:N-1)) 
    Fcoef = Fcoef*2.0_sp/(PI*j)
  END FUNCTION Fcoef  

  SUBROUTINE histogram(data,histo)
    USE nrutil, ONLY : arth, nrerror
    REAL(SP), DIMENSION(:), INTENT(IN) :: data
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: histo
    ! Given a vector of data and a number of bins, compute the histogram.
    ! If histo(:,1) is empty, use the data to estimate the range. 
    ! Otherwise, use the given range.
    REAL(SP) :: datamin, datamax, binwidth, Ncols, area
    INTEGER(I4B) :: sz, Nbins, j, jj
    INTEGER(I4B), DIMENSION(size(histo,1)) :: tally
    sz = size(data) 
    Nbins = size(histo,1); Ncols = size(histo,2)
    if ( Nbins < 2 ) call nrerror("histogram: histo must have at least 2 rows")
    if ( Ncols < 2 ) call nrerror("histogram: histo must have at least 2 columns")
    if ( all(histo(:,1) == 0.0) ) then
       datamin = minval(data); datamax = maxval(data)
       binwidth = (datamax-datamin)/real(Nbins-1,sp)
       histo(:,1) = arth(datamin,binwidth,Nbins)
    else
       datamin = histo(1,1); datamax = histo(Nbins,1) 
       binwidth = (datamax-datamin)/real(Nbins-1,sp)
    end if
    tally = 0
    do j = 1, sz
       jj = nint((data(j)-datamin)/binwidth) 
       tally(1+jj) = tally(1+jj)+1
    end do
    if ( sum(tally) /= sz ) call nrerror("histogram: something amiss with tally")
    area = (sz*Nbins)*binwidth
    histo(:,2) = real(tally,sp)
    if ( Ncols > 2 ) then
       histo(:,3) = sqrt(histo(:,2))
       histo(:,3) = histo(:,3)/area
    end if
    histo(:,2) = histo(:,2)/area
  END SUBROUTINE histogram

  FUNCTION kolm(N,Delta)
    INTEGER(I4B), INTENT(IN) :: N
    REAL(SP), INTENT(IN) :: Delta
    REAL(SP) :: kolm
    INTEGER(I4B), PARAMETER :: Maxit = 100
    ! Asymptotic two-sided Kolmogorov test in the form of
    ! M.A. Stephens, J. Royal Stat. Soc. B 32 (1970) 115.
    REAL(SP) :: sqn, a, alt, cut, add
    INTEGER(I4B) :: i
    sqn = sqrt(real(N,sp))
    a = -2.0_sp*(sqn*Delta+0.12_sp*Delta+0.11*Delta/sqn)**2
    alt = 2.0_sp
    kolm = 0.0_sp
    cut = 0.0_sp
    do i = 1, Maxit
       add = alt*exp(a*i**2)
       kolm = kolm+add
       if ( abs(add) <= cut ) return
       alt = -alt
       cut = abs(add)/1000.0_sp
    end do
    call nrerror("kolm: Kolmogorov test doesn't converge")
  END FUNCTION kolm

  FUNCTION PDF_s(Fser,a,b,x)
    REAL(SP), DIMENSION(:), INTENT(IN) :: Fser
    REAL(SP), INTENT(IN) :: a, b, x
    REAL(SP) :: PDF_s
    ! Compute the PDF as derivative of the CDF
    REAL(SP), DIMENSIOn(size(Fser)) :: aux
    REAL(SP) :: range, arg
    INTEGER(I4B) :: i, sz
    sz = size(Fser)
    range = b-a
    arg = (x-a)/range
    do i = 1, sz 
       aux(i) = PI*i*Fser(i)*cos(PI*i*arg)
    end do
    PDF_s = 1.0_sp+sum(aux)
    PDF_s = PDF_s/range
  END FUNCTION PDF_s

  FUNCTION PDF_v(Fser,a,b,x)
    REAL(SP), DIMENSION(:), INTENT(IN) :: Fser
    REAL(SP), INTENT(IN) :: a, b
    REAL(SP), DIMENSION(:), INTENT(IN) :: x
    REAL(SP), DIMENSION(size(x)) :: PDF_v
    ! Compute the PDF as derivative of the CDF
    REAL(SP), DIMENSION(size(x),size(Fser)) :: aux
    REAL(SP) :: range
    REAL(SP), DIMENSION(size(x)) :: arg
    INTEGER(I4B) :: i, sz
    sz = size(Fser)
    range = b-a
    arg = (x-a)/range
    do i = 1, sz 
       aux(:,i) = PI*i*Fser(i)*cos(PI*i*arg)
    end do
    PDF_v = 1.0_sp+sum(aux,dim=2)
    PDF_v = PDF_v/range
  END FUNCTION PDF_v

  FUNCTION Percentile_s(data,frac)
    REAL(SP), DIMENSION(:), INTENT(IN) :: data
    REAL(SP), INTENT(IN) :: frac
    REAL(SP) :: percentile_s
    ! Compute the percentile corresponding to frac in a ecdf.
    ! The format of ecdf is first column datapoints, second column
    ! value of the empirical cumulative distribution function.
    INTEGER(I4B) :: sz, ifrac
    if ( frac <= 0.0 .or. frac > 1.0 ) call nrerror("Percentile: 0 < frac <= 1")
    sz = size(data)
    ifrac = nint(frac*sz)
    percentile_s = data(ifrac)
  END FUNCTION Percentile_s

  FUNCTION Percentile_v(data,frac)
    REAL(SP), DIMENSION(:), INTENT(IN) :: data
    REAL(SP), DIMENSION(:), INTENT(IN) :: frac
    REAL(SP), DIMENSION(size(frac)) :: percentile_v
    ! Compute the percentile corresponding to frac in a ecdf.
    ! The format of ecdf is first column datapoints, second column
    ! value of the empirical cumulative distribution function.
    INTEGER(I4B) :: sz
    INTEGER(I4B), DIMENSION(size(frac)) :: ifrac
    if ( any(frac <= 0.0) .or. any(frac > 1.0) ) call nrerror("Percentile: 0 < frac <= 1")
    sz = size(data)
    ifrac = nint(frac*sz)
    percentile_v = data(ifrac)
  END FUNCTION Percentile_v

  SUBROUTINE ProbDistr(data,xval,Q,j,ECDF,CDF,EPDF)
    REAL(SP), DIMENSION(:), INTENT(IN) :: data, xval
    REAL(SP), INTENT(OUT) :: Q
    INTEGER(I4B), INTENT(OUT) :: j
    REAL(SP), DIMENSION(:), INTENT(OUT) :: ECDF, CDF, EPDF
    REAL(SP), PARAMETER :: Qcut = 0.5_sp
    INTEGER(I4B), PARAMETER :: MaxJ = 100
    ! Subroutine for analysis of probability distributions. Based on Berg and Harris, arXiv:0712.3852.
    ! Given a vector of sorted data, it computes the approximated probability distribution
    ! Returns also the result of the Kolmogorov test Q and the order of the Fourier series j.
    REAL(SP) :: step, a, b, delta
    REAL(SP), DIMENSION(MaxJ) :: Fser
    REAL(SP), DIMENSION(size(data)) :: R
    INTEGER(I4B) :: i, N, dim
    N = assert_eq(size(data),size(ECDF),size(CDF),'ProbDistr: N')
    dim = assert_eq(size(xval),size(EPDF),'ProbDistr: dim')
    step = 1.0_sp/real(N,sp)
    ECDF = arth(step,step,N) ! Empirical Cumulative Distribution Function
    ! Build the remainder function
    a = data(1); b = data(N)
    CDF = F0(a,b,data)   ! initial guess for the analytical CDF
    R = ECDF-CDF         ! remainder function
    FSer = 0.0_sp
    do j = 1, MaxJ
       ! Kolmogorov test
       delta = maxval(abs(ECDF-CDF))
       Q = kolm(N,delta)
       if ( Q >= Qcut ) exit
       ! Add one term of the Fourier expansion
       Fser(j) = Fcoef(ECDF,data,j)
       CDF = CDF+Fser(j)*sin(PI*j*(data-a)/(b-a))
    end do
    EPDF = PDF(Fser,a,b,xval)
  END SUBROUTINE ProbDistr

  SUBROUTINE JKProbDistr(data,xval,JKPDF)
    REAL(SP), DIMENSION(:), INTENT(IN) :: data, xval
    REAL(SP), DIMENSION(:,:), INTENT(OUT) :: JKPDF
    ! Produce the Probability Density Function according to the recipe 
    ! given by Berg and Harris, arXiv:0712.3852.
    ! The input data should be already sorted in ascending order.
    ! Errors are evaluated using a simple JackKnife procedure.
    REAL(SP), DIMENSION(size(data),size(JKPDF,1)) :: matr
    REAL(SP), DIMENSION(size(data)-1) :: v, ECDF, CDF
    REAL(SP) :: Q
    INTEGER(I4B) :: i, N, dim, m
    N = size(data)
    dim = assert_eq(size(xval),size(JKPDF,1),'JKProbDistr: dim')
    if (size(JKPDF,2)<2) call nrerror('JKProbDistr: wrong dimension of JKPDF')
    matr = 0.0_sp
    do i = 1, N
       if ( i == 1 ) then
          v = data(2:)
       else if ( i == N ) then
          v = data(:N-1)
       else
          v(:i-1) = data(:i-1); v(i:) = data(i+1:)
       end if
       call ProbDistr(v,xval,Q,m,ECDF,CDF,matr(i,:))
    end do
    JKPDF(:,1) = sum(matr,dim=1)/real(N,sp)
    matr = matr-spread(JKPDF(:,1),dim=1,ncopies=N)
    matr = matr*matr
    JKPDF(:,2) = sum(matr,dim=1)
    JKPDF(:,2) = JKPDF(:,2)*(N-1.0_sp)/real(N,sp)
    JKPDF(:,2) = sqrt(JKPDF(:,2))
  END SUBROUTINE

END MODULE amProbDistrib
