







!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : functionalsubs
! OBJECT: TYPE(funcsubs)
! USED  :
! DATE  : 2017-12-12
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!               Do some operations.
!
! STANDARD:
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  []
!
!
! avalable gets:
!                  [fun] get_string_upper(strIn)
!                        input string strIn and return the upper string.
!
!                  [fun] GetIntegerToCondenseSting(int4)
!                        return character::(16)
!
!                  [fun] ReportTime(t)   ! sees there are still some problem
!                        real*8::t
!                        return character::(64)
!
!                  [sub] LinuxMkdir(path)
!                        make directory. If it is existed, delete it (recersively) first.
! avalable Is:
!
! others:
!
!
!
!
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


module functionalsubs
    implicit none



    type::funcsubs


    contains
        procedure,nopass::get_string_upper=>to_upper
        procedure,nopass::GetIntegerToCondenseSting
        procedure,nopass::ReportTime
        procedure,nopass::LinuxMkdir
    end type

    private::to_upper
    private::LinuxMkdir

    contains



function to_upper(strIn) result(strOut)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
! Original author: Clive Page
! I found it from (20170628)
! https://stackoverflow.com/questions/10759375/how-can-i-write-a-to-upper-or-to-lower-function-in-f90
     implicit none
     character(len=*), intent(in) :: strIn
     character(len=len(strIn)) :: strOut
     integer :: i,j
     do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do
end function to_upper









character(16) function GetIntegerToCondenseSting(i)
               implicit none
               integer,intent(in)::i
               !----------------------
               write(GetIntegerToCondenseSting,*)i
               GetIntegerToCondenseSting = trim(adjustl(GetIntegerToCondenseSting))
endfunction


character(64) function ReportTime(t)
              implicit none
              real*8::t
              !-----------------------------
              integer::Time,Day,Hour,Min,Sec
              ! character(16)::temp
              Time = int(t)                ! ;write(*,*)time,t
              Day  = Time / 86400          ;  Time = Time - Day  * 86400
              Hour = Time / 3600           ;  Time = Time - Hour * 3600
              Min  = Time / 60             ;  Time = Time - Min  * 60
              Sec  = Time                 !;write(*,*)Day,Hour,Min,Sec,666

              ReportTime = ""
              !if (  Day.ne.0   )then
                ReportTime = trim(getint(Day))//" Day, "
              !endif
              !if (  Hour.ne.0   )then
                ReportTime = trim(adjustl(ReportTime))//trim(getint(Hour))//" Hou, "
              !endif
              !if (  Min.ne.0   )then
                ReportTime = trim(adjustl(ReportTime))//trim(getint(Min))//" Min, "
              !endif
              !if (  Min.ge.0   )then
                ReportTime = trim(adjustl(ReportTime))//trim(getint(Sec))//" Sec. "
              !endif
                      !WRITE(*,*)trim(getint(Sec))//" Sec. "
                      !WRITE(*,*)777,ReportTime
            contains
              character(16) function getint(i)
                implicit none
                integer,intent(in)::i
                !--------------
                write(getint,*)i
                getint = trim(adjustl(getint))
              endfunction

            endfunction






  subroutine LinuxMkdir(path)
    implicit none
    character(len=*),intent(in)::path
    !--------------------------------------
    logical::lexist,CheckDirExisted

    !---------------------------------------------
    !            make it and test it
    call system(("[[ ! -e ")//trim(adjustl(path))//" ]] && mkdir "//trim(adjustl(path)))
    !---------------------------------------------
    !            delete it
    call system(  "rm -r -f "//trim(adjustl(path)) )
    !---------------------------------------------
    !            make it
    call system("mkdir "//trim(adjustl(path)))
  endsubroutine





end module
