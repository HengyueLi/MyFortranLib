


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : CPTInterType
! OBJECT: TYPE(idata)
! USED  : FermionHamiltonian
! DATE  : 2018-01-02
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!           A interacting type used for the connection between type(Ham) and cluster method
!
! STANDARD:
!
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  [sub] SetIntoMatix(ns,Matrix,spini,spinj,IsLocal)
!                        integer,intent(in)      :: Ns
!                        complex*16,intent(inout):: Matrix(Ns,Ns)
!                        integer,intent(in)      :: spini,spinj
!                        logical,intent(in)      :: IsLocal
!
!                       set one term into matrix
!
! avalable gets:
!
!
! avalable is :
!                  ![fun] i
! others      :
!                  [sub] AppendDataToHam(H)
!                        1-basis index will transfor to 0-basis index. 
!                        for a input Type(Ham)::H, append all the interacting (type idata) into it.
!                        H maybe contains some other terms already, but here it does not check that.
!
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




module CPTInterType
  use FermionHamiltonian, only : Ham
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
    procedure,pass::AppendToHam
  endtype


  private::SetIntoMatix
  private::AppendToHam


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


  impure elemental subroutine AppendToHam(self,H)
    implicit none
    class(idata),intent(inout) :: self
    class(Ham),intent(inout)::H
    !---------------------------------
    Integer::para(8)
    para = self%Para
    para(1:2) = para(1:2) - 1
    call H%AppendingInteraction(  self%Itype  , Para , self%v)
  endsubroutine


endmodule
