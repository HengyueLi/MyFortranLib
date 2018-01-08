


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : LaPrimaryH
! OBJECT: TYPE(PH)
! USED  : CodeObject,functionalsubs,CPTInterType
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
!                  [sub] report(wtp)
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$







module LaPrimaryHUsedList
  use CPTInterType , only : Data => idata
  INCLUDE "../../ListStructure/ListStructure.finc"
endmodule






module  LaPrimaryH
  use LaPrimaryHUsedList    , only : Ilist => ListStru
  use CPTInterType
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

    procedure,pass::report



    ! procedure :: Copy
    ! generic :: assignment(=) => Copy
  endtype


  private::Initialization

  private::Allowed_S2,Allowed_S3,StartAppending,Append,EndAppending

  private::GetLen,GetIdata,GetState
  private::Copy
  private::report

contains

  subroutine Copy(a,b)
    implicit none
    class(PH),intent(out)::a
    class(PH),intent(in )::b
    !--------------------------------------------------
    a%state = b%state
    a%L     = b%L
  endsubroutine

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


  ! subroutine EndAppending(self)
  !   implicit none
  !   class(PH),intent(inout)::self
  !   !--------------------------------
  !   class(idata),pointer::p
  !   integer::jc
  !
  !
  !   call Allowed_S2(self)
  !   self%state = 3
  !   !-------------
  !   do jc = 1 , self%L%getlen()
  !      p => self%L%GetDataPointer(jc)
  !      if (  (p%P(1).ge.0) .and. (p%P(2).ge.0) .and. (p%P(3).ge.0) )then
  !      else
  !        write(self%getprint(),*)"ERROR: Only positiv-connected primary cell need to be input."
  !        write(self%getprint(),*)"Unknow part: i=",jc,"Pos:",p%p
  !        stop
  !      endif
  !   enddo
  ! endsubroutine
  subroutine EndAppending(self)
    implicit none
    class(PH),intent(inout)::self
    !--------------------------------
    class(idata),pointer::p
    integer::jc

    call Allowed_S2(self)
    self%state = 3
    !-------------
    call self%L%SetMark(1)
    do jc = 1 , self%L%getlen()
       p => self%L%GetMarkedPointerAndNext()
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


  subroutine report(self,wtp)
    implicit none
    class(PH),intent(inout):: self
    integer,intent(in)::wtp
    !---------------------------------
    character(32),allocatable::Disc(:)
    complex*16,allocatable::V(:)
    integer::jc,jc1
    class(idata),pointer::p
    Call self%L%SetMark(1)

    Allocate(   Disc(self%L%getlen()) ,   V(self%L%getlen())   )

    do jc = 1 , self%L%getlen()
      p => self%l%GetMarkedPointerAndNext()
      Disc(jc) = p%Disc
      V(jc)    = p%v      !;write(*,*)p%Disc,p%v
    enddo
    write(wtp,"(A67)",advance='no')"    +=============================================================+"
    write(wtp,*)" "
    write(wtp,'(A68)')"    |                 Hamitonian parameters:                      | "
    do jc = 1 , self%L%getlen()
       do jc1 = 1 , jc - 1
          if ( Disc(jc1)== Disc(jc) ) goto 100
       enddo
       call reportone(wtp,Disc(jc),v(jc))
   100 continue
    enddo
    deallocate(disc,v)
  contains

    subroutine reportone(wtp,disc,v)
      implicit none
      integer::wtp
      character(32)::disc
      complex*16::v
      !------------------------------
      write(wtp,100,advance='no')"    |",trim(adjustl(disc)) ," = " ,real(v), " , " ,imag(v),"         | "
      write(wtp,*)" "
  100 format( A5 , A32 ,A3 ,  f7.4  ,A3,  f7.4 ,A11)
    endsubroutine

  endsubroutine







endmodule








!
!
! #      Append all the Hamitonial terms in a primary cell.
! #The format looks like:
! #
! #        (Append)
! #         Disc                ! Character(32)
! #         type                ! Character(16)
! #         rp                  ! real*8
! #         ip                  ! real*8
! #         pos11 , pos12 , pos13  , par11 , par12 , par13
! #         pos21 , pos22 , pos23  , par21 , par22 , par23
! #         ...
! # where:
! # pos1,pos2,pos3  represent the position of Interaction
! # par1mpar2,par3  represent the parameters of Hamitonian terms
! #
! #
!
!
!
! (Append)
!  $t_c$
!  Hopping
!  1
!  0
!  1 , 0 , 0     ,  1 , 1 , 0
!  0 , 1 , 0     ,  1 , 1 , 0
!
!
!
! (Append)
!  $t_v$
!  Hopping
!  1
!  0
!  1 , 0 , 0     ,  2 , 2 , 0
!  0 , 1 , 0     ,  2 , 2 , 0
!
!
!
! (Append)
!  $D$
!  DiffOnSite
!  1
!  0
!  0 , 0 , 0     ,  1 , 2 , 0
!
!
!
! (Append)
!  $U$
!  OnSiteU
!  10.4
!  0
!  0 , 0 , 0     ,  1 , 1 , 0
!  0 , 0 , 0     ,  2 , 2 , 0
!
!
! (Append)
!  $V_1$
!  InterV
!  6
!  0
!  0 , 0 , 0     ,  1 , 2 , 0
!
!
!
! (Append)
!  $J$
!  Hund
!  0
!  0
!  0 , 0 , 0     ,  1 , 2 , 0
!
! (Append)
!  $I$
!  PairHopping
!  0
!  0
!  0 , 0 , 0     ,  1 , 2 , 0
!
!
! (Append)
!  $\mu$
!  "OnSite"
!  - (     )
!  0
!  0 , 0 , 0     ,  1 , 2 , 0
!
