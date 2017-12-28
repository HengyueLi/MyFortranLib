


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : LaPrimaryH
! OBJECT: TYPE(PH)
! USED  : CodeObject,functionalsubs
! DATE  : 2017-12-23
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            Fully discribe a Hamitonian of a Lattice system by choosing a primary cell.
!
!            When appending the interaction, only half of it will be inputted. Pos > 0 is kept.
!
! STANDARD:
!            *CALL Initialization(  )
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                   [sub] Initialization(  )
!
!                   [sub] StartAppending()
!
!                   [sub] Append(Disc,Pos,Itype,Ipara,Ivalu)
!                         character(DiscLen)  :: Disc          ( DiscLen = 32?  )
!                         integer             :: pos(3)         denote cluster position
!                         character(ITypeLen) :: Itype         ( ITypeLen = 16? )
!                         integer             :: Ipara(8)       1-basis for site index.  0,1 for spin up and down.
!                         complex*16          :: Ivalu
!                         Itype can be checked in module FermionHamiltonian
!                         Only for the terms that the first two elements in para(8) represent sites can be used.

!                   [sub] EndAppending()
!
!
! avalable gets:
!                   [fun] GetLen()
!                         get the total number of interacting
!
!                   [fun] GetIdata(i)
!                         type(idata)::GetIdata
!                         get the i-th interacting
!
!                   [fun] GetState()
!
!
!
! avalable is :
!                  ![fun] i
! others      :
!                  ![sub] p
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$





module LaprimaryHusedDatatype
  implicit none

  integer,parameter::DiscLen  = 32
  integer,parameter::ITypeLen = 16

  !-----------------
  !  contains interacting
  type::idata
    character(DiscLen)  :: Disc    ! discription of the interacting
    character(ITypeLen) :: Itype   !
    integer             :: Para(8) !
    complex*16          :: V       !
    integer             :: P(3)    ! position of connectedx cluster (PC).

  contains
    procedure,pass::SetIntoMatix
  endtype


  private::SetIntoMatix


contains

  !  add one signle particle term into matrix
  !  NOTICE: Matrix = Matrix + data
  !  recognized nteracting:
  !     [
  !     "SpinOnSite"  ,  "OnSite"  ,  "SpinHopping"  ,  "Hopping"
  !     ]
  !  logcal::IsLocal  ->  put on the conjugate term.
  subroutine SetIntoMatix(data,ns,Matrix,spini,spinj,IsLocal)
    use functionalsubs
    implicit none
    class(idata),intent(in) :: data
    integer,intent(in)      :: Ns
    complex*16,intent(inout):: Matrix(Ns,Ns)
    integer,intent(in)      :: spini,spinj
    logical,intent(in)      :: IsLocal
    !-----------------------------------
    TYPE(funcsubs)::f                         

    select case(adjustl(trim(f%get_string_upper(data%Itype))))
    case("SPINONSITE")               !     "SpinOnSite"
      if ( (spini==spinj) .and.  (spini==data%Para(3)) ) then
         Matrix(data%Para(1),data%Para(1)) = Matrix(data%Para(1),data%Para(1)) + data%v
      endif
    case("ONSITE")                   !     "OnSite"
      if ( (spini==spinj)                              ) then
         Matrix(data%Para(1),data%Para(1)) = Matrix(data%Para(1),data%Para(1)) + data%v
      endif
    case("SPINHOPPING")             !   "SpinHopping"
      if ( (spini==spinj).and.  (spini==data%Para(3))  ) then
        Matrix(data%Para(1),data%Para(2)) = Matrix(data%Para(1),data%Para(2)) + data%v
        if ( IsLocal ) then!------local term
          Matrix(data%Para(2),data%Para(1)) = Matrix(data%Para(2),data%Para(1)) + conjg(data%v)
        endif
      endif
    case("HOPPING")                !      "Hopping"
      if ( (spini==spinj)                             ) then
        Matrix(data%Para(1),data%Para(2)) = Matrix(data%Para(1),data%Para(2)) + data%v
        if ( IsLocal )then !------local term
          Matrix(data%Para(2),data%Para(1)) = Matrix(data%Para(2),data%Para(1)) + conjg(data%v)
        endif
      endif
    end select

  endsubroutine

endmodule




module LaPrimaryHUsedList
  use LaprimaryHusedDatatype , only : Data => idata
  INCLUDE "../../ListStructure/ListStructure.finc"
endmodule






module  LaPrimaryH
  use LaPrimaryHUsedList    , only : Ilist => ListStru
  use LaprimaryHusedDatatype
  use CodeObject
  implicit none


  type,extends(object)::PH
    private
    integer    ::state = 0
    type(Ilist)::L


  contains
    procedure,pass::Initialization


    procedure,pass::StartAppending
    procedure,pass::Append
    procedure,pass::EndAppending

    procedure,pass::GetLen
    procedure,pass::GetIdata
    procedure,pass::GetState
  endtype


  private::Initialization

  private::Allowed_S2,Allowed_S3,StartAppending,Append,EndAppending

  private::GetLen,GetIdata,GetState

contains

  subroutine Initialization(self)
    implicit none
    class(PH),intent(inout)::self
    !--------------------------------------------------
    call self%SetInitiated(.true.)

    self%state = 1

  endsubroutine


  SUBROUTINE StartAppending(self)
    implicit none
    class(PH),intent(inout)::self
    !--------------------------------------------------
    call self%L%Initialization(self%getprint())
    self%state = 2
  endsubroutine

  subroutine Append(self,Disc,Pos,Itype,Ipara,Ivalu)
    implicit none
    class(PH),intent(inout)::self
    character(DiscLen),intent(in)::Disc
    integer,intent(in)::Pos(3)
    character(ITypeLen),intent(in)::Itype
    integer,intent(in)::Ipara(8)
    complex*16,intent(in)::Ivalu
    !--------------------------------------------------
    class(idata),pointer::p
    call Allowed_S2(self)
    allocate(p)
    p%Disc  = Disc
    p%Itype = Itype
    p%Para  = Ipara
    p%P     = pos
    p%V     = Ivalu
    Call self%L%append(p)
  endsubroutine


  subroutine EndAppending(self)
    implicit none
    class(PH),intent(inout)::self
    !--------------------------------
    class(idata),pointer::p
    integer::jc


    call Allowed_S2(self)
    self%state = 3
    !-------------
    do jc = 1 , self%L%getlen()
       p => self%L%GetDataPointer(jc)
       if (  (p%P(1).ge.0) .and. (p%P(2).ge.0) .and. (p%P(3).ge.0) )then
       else
         write(self%getprint(),*)"ERROR: Only positiv-connected primary cell need to be input."
         write(self%getprint(),*)"Unknow part: i=",jc,"Pos:",p%p
         stop
       endif
    enddo
  endsubroutine



  subroutine Allowed_S2(self)
    implicit none
    class(PH),intent(inout)::self
    !--------------------------------
    call self%CheckInitiatedOrStop()
    if (self%state.eq.2)then
    else
      write(self%getprint(),*)"ERROR: PH is not allowed to Append interacting"
      write(self%getprint(),*)"call StartAppending first"
      stop
    endif
  endsubroutine

  subroutine Allowed_S3(self)
    implicit none
    class(PH),intent(inout)::self
    !--------------------------------
    call self%CheckInitiatedOrStop()
    if (self%state.eq.3)then
    else
      write(self%getprint(),*)"ERROR: In PH, EndAppending is not called yet."
      stop
    endif
  endsubroutine

  integer function GetLen(self)
    implicit none
    class(PH),intent(inout)::self
    !--------------------------------
    call Allowed_S3(self)
    GetLen = self%l%getlen()
  endfunction

  type(idata) function GetIdata(self,i)
    implicit none
    class(PH),intent(inout):: self
    integer,intent(in)     :: i
    !--------------------------------
    class(idata),pointer::p
    call Allowed_S3(self)                !;
    p => self%l%GetDataPointer(i)        !;write(*,*)i,associated(p)
    GetIdata = p
  endfunction


  integer function GetState(self)
    implicit none
    class(PH),intent(inout):: self
    !---------------------------------
    GetState = self%state
  endfunction







endmodule
