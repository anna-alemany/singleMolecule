! Utilities for manipulating text files.
! Written by Alessandro Mossa.

! CountRows, FindText, IntToText, ReadData
! CountColumns

MODULE amFileUtil
  USE nrtype; USE nrutil, ONLY : nrerror
  IMPLICIT NONE

  INTERFACE ReadData
     MODULE PROCEDURE ReadData_s1, ReadData_d1, ReadData_s2, ReadData_d2
  END INTERFACE

CONTAINS

  FUNCTION CountRows(unit,tag)
    INTEGER(I4B), INTENT(IN) :: unit
    CHARACTER(1), INTENT(IN), OPTIONAL :: tag
    INTEGER(I4B) :: CountRows
    ! Count the number of rows in the file specified by unit.
    ! Optionally, one can specify the tag 'n' to count only the 
    ! lines that begin with a number, or the tag '#' to count 
    ! only the lines that do not begin with '#'. If tag is absent, 
    ! then all the lines which contains at least a character will be 
    ! counted (but not blank lines). Rewind the unit when done.
    REAL(SP) :: anynumber
    INTEGER(I4B) :: ReadStatus
    CHARACTER(1) :: anything
    CountRows = 0
    if ( present(tag) .and. tag == 'n' ) then
       do 
          read(unit,*,iostat=ReadStatus) anynumber
          if ( ReadStatus == 0 ) then
             CountRows = CountRows+1
          else if ( ReadStatus < 0 ) then
             exit
          else if ( ReadStatus > 0 ) then
             cycle
          end if
       end do
    else if ( present(tag) .and.( (tag == '#').or.(tag=='u').or.(tag=='f')) ) then
       do 
          read(unit,*,iostat=ReadStatus) anything
          if ( ReadStatus == 0 ) then
             if ( anything /= tag) then
                CountRows = CountRows+1
             else
                cycle
             end if
          else if ( ReadStatus > 0 ) then
             call nrerror("CountRows")
          else if ( ReadStatus < 0 ) then
             exit
          end if
       end do
    else
       do
          read(unit,*,iostat=ReadStatus) anything
          if ( ReadStatus == 0 ) then
             CountRows = CountRows+1
          else if ( ReadStatus < 0 ) then
             exit
          else if ( ReadStatus > 0 ) then
             call nrerror("CountRows")
          end if
       end do
    end if
    rewind(unit)
  END FUNCTION CountRows

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
       if ( ReadStatus > 0 ) call nrerror('FindText')
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

  FUNCTION IntToText(n)
    INTEGER(I4B), INTENT(IN) :: n
    CHARACTER(11) :: IntToText
    ! Convert a non-negative integer number into a string
    CHARACTER(11) :: text
    if ( n >= 0 ) then
       WRITE(text,'(I11)') n
       IntToText = adjustl(text)
    else
       call nrerror('IntToText requires non-negative argument')
    end if
  END FUNCTION IntToText

  FUNCTION ParentFolder(address)
    CHARACTER(*), INTENT(IN) :: address
    CHARACTER(len(address)) :: ParentFolder
    ! Find the Parent Folder for the file specified by address.
    INTEGER(I4B) :: i, L
    L = len(address)
    do i = L, 1, -1
       if ( address(i:i) == '/' ) exit 
    end do
    if ( i == 0 ) then
       ParentFolder = "./"
    else
       ParentFolder = address(1:i)
    end if
  END FUNCTION ParentFolder

  SUBROUTINE ReadData_s1(unit,x)
    INTEGER(I4B), INTENT(IN) :: unit
    REAL(SP), DIMENSION(:), INTENT(OUT) :: x
    ! Read the data from the file specified by unit into the vector x.
    REAL(SP) :: data
    INTEGER(I4B) :: ReadStatus, i
    i = 0
    do
       read(unit,*,iostat=ReadStatus) data
       if ( ReadStatus < 0 ) exit
       if ( ReadStatus > 0 ) call nrerror('ReadData_s1: error reading datafile')
       i = i+1
       x(i) = data
    end do
    if ( i /= size(x) ) call nrerror('ReadData_s1: dimension mismatch')
  END SUBROUTINE ReadData_s1

  SUBROUTINE ReadData_d1(unit,x)
    INTEGER(I4B), INTENT(IN) :: unit
    REAL(DP), DIMENSION(:), INTENT(OUT) :: x
    ! Read the data from the file specified by unit into the vector x.
    REAL(DP) :: data
    INTEGER(I4B) :: ReadStatus, i
    i = 0
    do
       read(unit,*,iostat=ReadStatus) data
       if ( ReadStatus < 0 ) exit
       if ( ReadStatus > 0 ) call nrerror('ReadData_d: error reading datafile')
       i = i+1
       x(i) = data
    end do
    if ( i /= size(x) ) call nrerror('ReadData_d: dimension mismatch')
  END SUBROUTINE ReadData_d1

  SUBROUTINE ReadData_s2(unit,x)
    INTEGER(I4B), INTENT(IN) :: unit
    REAL(SP), DIMENSION(:,:), INTENT(OUT) :: x
    ! Read the data from the file specified by unit into the matrix x.
    REAL(SP), DIMENSION(size(x,2)) :: data
    INTEGER(I4B) :: ReadStatus, i
    i = 0
    do
       read(unit,*,iostat=ReadStatus) data
       if ( ReadStatus < 0 ) exit
       if ( ReadStatus > 0 ) call nrerror('ReadData_s: error reading datafile')
       i = i+1
       x(i,:) = data
    end do
    if ( i /= size(x,1) ) call nrerror('ReadData_s: dimension mismatch')
  END SUBROUTINE ReadData_s2

  SUBROUTINE ReadData_d2(unit,x)
    INTEGER(I4B), INTENT(IN) :: unit
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: x
    ! Read the data from the file specified by unit into the matrix x.
    REAL(DP), DIMENSION(size(x,2)) :: data
    INTEGER(I4B) :: ReadStatus, i
    i = 0
    do
       read(unit,*,iostat=ReadStatus) data
       if ( ReadStatus < 0 ) exit
       if ( ReadStatus > 0 ) call nrerror('ReadData_d: error reading datafile')
       i = i+1
       x(i,:) = data
    end do
    if ( i /= size(x,1) ) call nrerror('ReadData_d: dimension mismatch')
  END SUBROUTINE ReadData_d2

  FUNCTION CountColumns(unit)
    INTEGER(I4B), INTENT(IN) :: unit
    ! Troba el nombre de columnes que hi ha a l'arxiu en qüestió.
    ! La primera és per defecte un caràcter, i les altres són reals.
    INTEGER(I4B) :: i, CountColumns, k, ReadStatus
    REAL(SP) :: num
    CHARACTER(100) :: paraula
    rewind(unit)
    i = 0
    do 
      read(unit,*,iostat=ReadStatus) paraula, (num, k=1,i+1)
      backspace(unit)
      if (ReadStatus == 0) then
	i = i+1
      else if (Readstatus /= 0) then
	exit
      end if
    end do
    CountColumns = i; rewind(unit)
  END FUNCTION CountColumns

END MODULE amFileUtil
