! DataSplit.
! Llegeix les dades i amb el COMfile separa les trajectories de pulling per u i f
! guarda: cnt, y, t, fy
! Deixa 2 linies en blanc de separació quan hi ha alguna pausa (index 0, 1, 2, ... en el gnuplot)

! creat: 110829 (DataSplitPulling reutilitzat)
! modificat: 130819

! abans de fer córrer això és interessant posar al terminal:
! for file in */*COM.txt; do awk 'FNR==NR{if (NR>=115) {a[NR]=$1} else {print $0}} {if (NR>=115) {if ($1/1==0) {print a[FNR-1], $0} else {print $0}}}' $file > $file-new; mv $file-new $file; done

! for file in 060813/*COM.txt; do 
! awk 'FNR==NR{if (FNR>114) {
! if ($1/1==0) {
! print cnt[FNR-1], $0
! } else {
! print $0; cnt[FNR]=$1
! }
! } else
! {print $0}
! }' $file > $file-kk;  mv $file-kk $file; 
! done


MODULE Parametres
  USE nrtype
  IMPLICIT NONE
  INTEGER(I4B), PARAMETER :: nocol = 50
  CHARACTER(17), DIMENSION(nocol), PARAMETER :: channels = (/'*CycleCount/n____','*A_PsdX__________','*A_PsdY__________',& ! 3
         '*A_PsdSum________','*A_Iris__________','*A_LeverX________','*A_LeverY________','*A_LeverSum______',& 
         '*B_PsdX__________','*B_PsdY__________','*B_PsdSum________','*B_Iris__________','*B_LeverX________',&
         '*B_LeverY________','*B_LeverSum______','*A_Temperature___','*B_Temperature___','*X_force_________',&
         '*A_forceX________','*B_forceX________','*Y_force_________','*A_forceY________','*B_forceY________',&
         '*Z_force_________','*A_forceZ________','*B_forceZ________','*Tension_________','*A_dist-X________',&
         '*A_dist-Y________','*B_dist-X________','*B_dist-Y________','*Rel_dist_X______','*Rel_dist_Y______',&
         '*Motor_X_________','*Motor_Y_________','*Motor_Z_________','*MotorVel_X______','*MotorVel_Y______',& 
         '*MotorVel_Z______', & ! 32
	 '*BeadDis_X_______','*BeadDis_Y_______',&
	 '*A_trapPIC_Err___','*B_trapPIC_Err___','*comPIC_Err______','*A_CycleCount____',&
         '*B_CycleCount____','*Motor_CycleCount','*time(sec)_______','*Status__________','*Null____________'/) 
END MODULE Parametres

PROGRAM Data_Split
  USE nrtype; USE nrutil, ONLY : nrerror
  USE amFileUtil, ONLY : IntToText
  USE Parametres
  IMPLICIT NONE
  REAL(SP) :: yt, fy, fyA, fyB, t, t0, t00
  REAL(SP), DIMENSION(:), ALLOCATABLE :: rawdata
  INTEGER(I4B) :: pos_fy, pos_fyA, pos_fyB, pos_yA, pos_yB, pos_t
  INTEGER(I4B) :: a, counter0, counter1, counter2, counter3, cnt, cnt0, NU, NF
  INTEGER(I4B) :: OpenStatus, ReadStatus, ii, jj
  INTEGER(I4B), DIMENSION(nocol) :: key
  CHARACTER(20) :: label, fold*1
  CHARACTER(256) :: inputfile, molecula
  CHARACTER(1), DIMENSION(8), PARAMETER :: alphabet = (/'A','B','C','D','E','F','G','H'/)

!  print *, "Please enter the path of the input file"
!  print *, '(no spaces allowed, enclose between quotes "" if slashes / are required):'
!  read *, inputfile 
  inputfile = 'input_split.txt'
  open(12, file=trim(inputfile), status='old', iostat=OpenStatus)
  if (OpenStatus /= 0) call nrerror("ERR001: Inputfile not found.")
  outerloop: do
    read(12,*,iostat=ReadStatus) molecula
    if (ReadStatus < 0) exit outerloop
    if (ReadStatus > 0) call nrerror("ERR002: Error reading the Inputfile")
    if (molecula(:1)=='#') cycle
    ! Read the COMfile
    open(14, file=trim(molecula)//"COM.txt", status='old', iostat=OpenStatus)
    if ( OpenStatus /= 0 ) call nrerror("ERR003: COM file not found")
    print *, "Opening ",trim(molecula)//"COM.txt..." 
    call FindText(14,'*CycleCount/n____',ReadStatus)
    if ( ReadStatus /= 0 ) call nrerror("ERR004: Error reading COM file")
    key = ReadChannels(14) ! read from the COMfile which columns are present
    allocate(rawdata(sum(key)-1))
    call FindText(14,'Force-limit',ReadStatus) ! aquest protocol SEMPRE comença amb un unfolding
    ! Now we are ready to go and open the datafiles
    a = 1  ! index over the letters of the alphabet
    NU = 0; NF = 0  ! counters over folding and unfolding events
    open(16, file=trim(molecula)//alphabet(a)//".txt", status='old', iostat=OpenStatus)
    if ( OpenStatus /= 0 ) call nrerror('ERR005: datafile not found')
    open(15, file=trim(molecula)//"/LOG.txt", status='replace')
    read(14,*,iostat=ReadStatus) counter0, label
    counter2 = counter0-1; counter3 = counter0-1
    readCOM: do
      read(14,*,iostat=ReadStatus) counter1, label
      if ( ReadStatus == 0 ) then
	select case(trim(label))
	case('minForce-limit')
	  fold = 'f'
	case('maxForce-limit')
	  fold = 'u'
	case('ReleasePause-begin')
	  counter2 = counter1; cycle readCOM
	case('ReleasePause-end')
	  counter3 = counter1; cycle readCOM
	case('StretchPause-begin')
	  counter2 = counter1; cycle readCOM
	case('StretchPause-end')
	  counter3 = counter1; cycle readCOM
	case default
	  cycle readCOM       ! skip lines that do not interest
	end select
      else if ( ReadStatus > 0 ) then
	print *, counter1, label
	call nrerror('Usar BASH abans d''aquest programa')
      else if ( ReadStatus < 0 ) then
	  exit readCOM
      end if
      call FindCounter(16,counter0,ReadStatus)
      read(16,*) cnt, rawdata
      backspace(16)
      select case(fold)
      case('u')   
	NU = NU+1
	open(13,file=trim(molecula)//"/events_raw/u"//trim(IntToText(NU))//".txt",status='replace')
	write(15,*) "u"//trim(IntToText(NU)), counter0, counter1-1
	ii = 0 ; jj = 0
      case('f')
	NF = NF+1
	open(13,file=trim(molecula)//"/events_raw/f"//trim(IntToText(NF))//".txt",status='replace')
	if (counter2 > counter0) open(31,file=trim(molecula)//"/events_raw/p"//trim(IntToText(NF))//".txt",status='replace')
	write(15,*) "f"//trim(IntToText(NF)), counter0, counter1-1
	ii = 0; jj = 0
      case default
	call nrerror("ERR007: fold must be either 'u' or 'f'")
      end select
! Escriptura de dades a events_raw
      copydata: do
	read(16,*,iostat=ReadStatus) cnt, rawdata
	if ( ReadStatus > 0 ) cycle
	if ( ReadStatus < 0 ) then 
	  call OpenNext(trim(molecula),16,a)
	  cycle
	end if
	if ( cnt >= counter1 ) then
	  backspace(16)
	  exit copydata
	end if
	! Determination of yt (control parameter)
	if ( key(29) == 1 .and. key(31) == 1 ) then
	  pos_yA = sum(key(1:29))-1
	  pos_yB = sum(key(1:31))-1
	  yt = 0.5_sp*(rawdata(pos_yA)+rawdata(pos_yB))
	else
	  call nrerror("ERR006: data don't allow to determine y-distance")
	end if
	! Determination of the y-force 
	if ( key(21) == 1 ) then
	  pos_fy = sum(key(1:21))-1 
	  fy = -rawdata(pos_fy)
	  fyA = 0.5_sp*fy; fyB = 0.5_sp*fy
	else if ( key(22) == 1 .and. key(23) == 1 ) then
	  pos_fyA = sum(key(1:22))-1
	  pos_fyB = sum(key(1:23))-1
	  fyA = -rawdata(pos_fyA)
	  fyB = -rawdata(pos_fyB)
	  fy = fyA+fyB
	else
	  call nrerror("ERR009: data don't allow to determine y-force")
	end if 
! 	   Determination of time
	if ( key(nocol-2) == 1) then
	  pos_t = sum(key(1:nocol-2))-1
	  t = rawdata(pos_t)
	else
	  call nrerror("ERR: No hi ha temps a guardar!")
	end if
	ii = ii+1
	if (ii == 1) then
	  t0 = t
	end if
	t = t-t0
	if ((counter2 < counter0).and.(counter3 < counter0)) then
	  write(13,*) cnt, yt, t, fy
	else
	  if ((cnt < counter2).and.(cnt < counter3)) then
	    write(13,*) cnt, yt, t, fy
	  else if ((cnt > counter2).and.(cnt < counter3)) then
	    write(13,*) cnt, yt, t, fy
	    jj = jj+1
	    if (jj==1) t00 = t
	    write(31,*) cnt, yt, t-t00, fy
	  else if (cnt > counter3) then
	    write(13,*) cnt, yt, t, fy
	    counter2 = counter0-1; counter3 = counter0-1
	  end if
	end if
      end do copydata
      close(13)
      counter0 = counter1
    end do readCOM
    close(16)
    close(14)
    deallocate(rawdata)
  end do outerloop
  close(12)

CONTAINS

  SUBROUTINE FindCounter(unit,n,ReadStatus)
    INTEGER(I4B), INTENT(IN) :: unit, n
    INTEGER(I4B), INTENT(OUT) :: ReadStatus
    ! find in unit the line beginning with the counter n
    INTEGER(I4B) :: cntr
    do
       read(unit,*,iostat=ReadStatus) cntr
       if ( ReadStatus < 0 ) exit    ! end-of-file
       if ( ReadStatus > 0 ) cycle   ! skip the lines with comments
       if ( cntr >= n ) then
          backspace(unit)
          exit
       end if
    end do
  END SUBROUTINE FindCounter

  SUBROUTINE FindText(unit, text, ReadStatus)
    INTEGER(I4B), INTENT(IN) :: unit
    CHARACTER(*), INTENT(IN) :: text
    INTEGER(I4B), INTENT(OUT) :: ReadStatus
    ! find in the file connected to unit the text specified
    INTEGER(I4B) :: sz, i
    CHARACTER(1) :: char
    sz = len(text)
    i = 1
    do 
       read(unit,'(A1)',advance='no',iostat=ReadStatus) char 
       if ( ReadStatus == -1 ) exit   ! end-of-file
       if ( ReadStatus == -2 ) then   ! end-of-line
          i = 1
          cycle 
       end if
       if ( ReadStatus > 0 ) call nrerror('ERR009: FindText')
       if ( char /= text(i:i) ) then
          i = 1
          cycle
       end if
       if ( i == sz ) then
          backspace(unit) 
          exit
       else
          i = i+1
       end if
    end do
  END SUBROUTINE FindText

  SUBROUTINE OpenNext(address,unit,progr)
    CHARACTER(*), INTENT(IN) :: address
    INTEGER(I4B), INTENT(IN) :: unit
    INTEGER(I4B), INTENT(INOUT) :: progr
    ! Open next datafile
    INTEGER(I4B) :: OpenStatus
    close(unit)
    progr = progr+1
    open(unit, file=address//alphabet(progr)//".txt", status='old', iostat=OpenStatus)
print *, address//alphabet(progr)//".txt"
    if ( OpenStatus /= 0 ) call nrerror('ERR008: datafile not found')
  END SUBROUTINE OpenNext

  FUNCTION ReadChannels(unit)
    INTEGER(I4B), INTENT(IN) :: unit
    INTEGER(I4B), DIMENSION(nocol) :: ReadChannels
!     CHARACTER(17), DIMENSION(50), PARAMETER :: channels = (/'*CycleCount/n____','*A_PsdX__________','*A_PsdY__________',&
!          '*A_PsdSum________','*A_Iris__________','*A_LeverX________','*A_LeverY________','*A_LeverSum______',&
!          '*B_PsdX__________','*B_PsdY__________','*B_PsdSum________','*B_Iris__________','*B_LeverX________',&
!          '*B_LeverY________','*B_LeverSum______','*A_Temperature___','*B_Temperature___','*X_force_________',&
!          '*A_forceX________','*B_forceX________','*Y_force_________','*A_forceY________','*B_forceY________',&
!          '*Z_force_________','*A_forceZ________','*B_forceZ________','*Tension_________','*A_dist-X________',&
!          '*A_dist-Y________','*B_dist-X________','*B_dist-Y________','*Rel_dist_X______','*Rel_dist_Y______',&
!          '*Motor_X_________','*Motor_Y_________','*Motor_Z_________','*MotorVel_X______','*MotorVel_Y______',&
!          '*MotorVel_Z______', &
! 	 '*BeadDis_X_______','*BeadDis_Y_______',&
! 	 '*A_trapPIC_Err___','*B_trapPIC_Err___','*comPIC_Err______','*A_CycleCount____',&
!          '*B_CycleCount____','*Motor_CycleCount','*time(sec)_______','*Status__________','*Null____________'/) 
    ! For reading the COMfile
    INTEGER(I4B) :: ReadStatus, i, n
    CHARACTER(17) :: anything
    do i = 1, nocol
       read(unit,'(A17,t19,I1)',iostat=ReadStatus) anything, n
       if ( ReadStatus /= 0 ) call nrerror("ReadValue")
       if ( anything /= channels(i) ) call nrerror("ReadValue")
       ReadChannels(i) = n
    end do
  END FUNCTION ReadChannels

END PROGRAM Data_Split
