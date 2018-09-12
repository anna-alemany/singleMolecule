! Representa gràficament les trajectòries del LOG que no comencen per #.

PROGRAM PlotLog
  IMPLICIT NONE
  CHARACTER(len=5) :: traj0, traj1
  CHARACTER(len=100) :: molecula
  INTEGER :: ReadStatus, OpenStatus


do
  print *, 'Molecula a representar?'; read *, molecula
  if ((molecula(1:1)=='e').or.(molecula(1:1)=='s').or.(molecula(1:1)=='q')) stop
  open(unit=12, file=trim(molecula)//'/LOG.txt', status='old', iostat=OpenStatus)
  if (OpenStatus/=0) then
    print *, 'Problemes en obrir el '//trim(molecula)//'/LOG.txt'
    stop
  end if
  open(unit=10, file=trim(molecula)//'/events_raw/PlotLog.pl', status='replace')

  WritePlot: do
    write(10,*) 'pause-1'
    read(12,*,iostat=ReadStatus) traj0
    if (ReadStatus > 0) then
      print *, 'Problemes llegint el LOG.txt (per segon cop!)'
      stop
    else if (ReadStatus < 0) then
      exit WritePlot
    end if
!     if (traj0(1:1)=='#') cycle
    if (traj0(1:1) /= 'f') cycle
    write(10,*) 'pl '''//trim(molecula)//'/events_raw/'//trim(traj0)//'.txt'' us 2:4 w l title '''//trim(traj0)//''''

    read(12,*,iostat=ReadStatus) traj1
    if (ReadStatus > 0) then
      print *, 'Problemes llegint el LOG.txt (per segon cop!)'
      stop
    else if (ReadStatus < 0) then
      exit WritePlot
    end if
    if (traj1(1:1)=='#') cycle
    write(10,*) 'repl '''//trim(molecula)//'/events_raw/'//trim(traj1)//'.txt'' us 2:4 w l title '''//trim(traj1)//''''
    
  end do WritePlot
  close(10); close(12)
end do

END PROGRAM PlotLog


!   INTEGER :: refF, refU, it
!   INTEGER :: OpenStatus, ReadStatus
!   CHARACTER(len=5) :: traj, referenciaF, referenciaU
!   CHARACTER(len=1) :: label
!   CHARACTER(len=256) :: molecula
!   it = 0; refF = 0; refU = 0 
!   escrivintShift: do
!     it = it+1
!     read(12,*,iostat=ReadStatus) traj!, N0, sh0, N1, sh1, n
!     if (ReadStatus > 0) then
!       print *, 'Problemes llegint el LOG.txt (per segon cop!)'
!       stop
!     else if (ReadStatus < 0) then
!       exit escrivintShift
!     end if
!     label = traj(1:1)
!     select case (label)
!     case('u')
!       refU = refU+1; referenciaU = trim(traj);
!     case('f')
!       refF = refF+1; referenciaF = trim(traj);
!     case('#')
!       cycle
!     case default
!       print *, 'Problema llegint algun nom de trajectoria: '//trim(traj); stop
!     end select
!     if ((refF == 1).and.(refU == 1)) then
!       write(10,*)'plot '''//trim(molecula)//'/events_raw/'//trim(referenciaF)//'.txt'' us 2:4 w l'
!       write(10,*)'replot '''//trim(molecula)//'/events_raw/'//trim(referenciaU)//'.txt'' us 2:4 w l; pause-1'
!       refF = 0; refU = 0
!     else if ((refF == 2).and.(refU == 0)) then
!       write(10,*)'plot '''//trim(molecula)//'/events_raw/'//trim(referenciaF)//'.txt'' us 2:4 w l; pause-1'
!       refF = 0; refU = 0
!     else if ((refF == 0).and.(refU == 2)) then
!       write(10,*)'plot '''//trim(molecula)//'/events_raw/'//trim(referenciaU)//'.txt'' us 2:4 w l; pause-1'
!       refF = 0; refU = 0
!     else if ((refF > 2).or.(refU > 2)) then
!       print *, 'Algun error molt estrany! refF = ',refF,' i refU = ',refU
!       stop
!     end if
!   end do escrivintShift
