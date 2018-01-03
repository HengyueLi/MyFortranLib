

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : CodeObject
! OBJECT: TYPE(Object)
! USED  :
! DATE  : 2017-12-24
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source :
!------------------
! DESCRIPTION:
!            A basis class used for coding.
!
! STANDARD:
!
!
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                   [sub] SetPrint(print)
!                         integer::print
!
!                   [sub] SetShow(show)
!                         integer::show
!
!                   [sub] SetInitiated(TF)
!                         logical::TF
!
! avalable gets:
!                   [fun] GetPrint()
!
!                   [fun] GetShow()
!
!
!
! avalable is :
!                   [fun] IsInitiated()
!
!
! others      :
!                   [sub] CheckInitiatedOrStop()
!                         check if this object is initiated, if not, stop and report.
!
!
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


module CodeObject
  implicit none

  type::Object
    private
    logical :: initiated = .false.


    !------------------------------
    integer :: print = 6
    integer :: show  = 0

  contains




    procedure,pass::CheckInitiatedOrStop

    procedure,pass::IsInitiated
    procedure,pass::SetPrint
    procedure,pass::SetShow
    procedure,pass::SetInitiated
    procedure,pass::GetPrint
    procedure,pass::GetShow

  endtype



  private::CheckInitiatedOrStop
  private::IsInitiated
  private::SetPrint
  private::SetInitiated
  private::GetPrint
  private::GetShow

contains


  logical function IsInitiated(self)
    implicit none
    class(Object),intent(in)::self
    !----------------------------------
    IsInitiated = self%initiated
  endfunction




  subroutine SetPrint(self,print_)
    implicit none
    class(Object),intent(inout)::self
    integer,intent(in)::print_
    !----------------------------------
    self%print = print_
  endsubroutine

  subroutine SetShow(self,show_)
    implicit none
    class(Object),intent(inout)::self
    integer,intent(in)::show_
    !----------------------------------
    self%show = show_
  endsubroutine

  subroutine SetInitiated(self,TF)
    implicit none
    class(Object),intent(inout)::self
    logical,intent(in)::TF
    !----------------------------------
    self%initiated = TF
  endsubroutine

  integer function GetPrint(self)
    implicit none
    class(Object),intent(in)::self
    !----------------------------------
    GetPrint = self%print
  endfunction

  integer function GetShow(self)
    implicit none
    class(Object),intent(in)::self
    !----------------------------------
    GetShow = self%show
  endfunction





  subroutine CheckInitiatedOrStop(self)
    implicit none
    class(Object),intent(in)::self
    !----------------------------------
    if (.not.self%initiated)then
       write(self%print,*)"ERROR: object is used before Initialization. Program exit."
       stop
    endif
  endsubroutine




End Module CodeObject
