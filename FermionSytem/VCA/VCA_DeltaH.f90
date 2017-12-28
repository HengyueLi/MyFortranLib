

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : VCA_DeltaH
! OBJECT: TYPE(VCAdH)
! USED  : CodeObject,LatticeConfig,LaPrimaryH
! DATE  : 2017-12-28
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            Saving VCA meanfiled.
!
! STANDARD:
!            *CALL Initialization( PrH,LaC,print_,show_ )
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                   [sub] I
!
!                   [sub] StartAppending()
!
!                   [sub] EndAppending()
!
!                   [sub] Append(data)
!                         type(idata)::data
!                         append data into VCAdH.
!                         In data, all index are 1-basied and represent the index in LC.
!
!                   [sub] AppendDeltaAF(B,Qm,discription)
!                         real*8,intent(in)::B,Qm(3)
!                         character(DiscLen),intent(in)::discription   DiscLen=32?
!
!                         Delta =   sum_{r,spin} (-1)^spin Exp(i Qm * r )
!
!                   [sub] AppendOrbitalCDW(E,Qm,oi,oj,discription)
!                         real*8,intent(in)::E,Qm(3)
!                         integer,intent(in)::Oi,Oj
!                         character(DiscLen),intent(in)::discription   DiscLen=32?
!
!                         Delta = sum_{r} Exp(i Qm * r ) * (n_oi - n _oj)
!                             notice the order of oi,oj
!
!
!
! avalable gets:
!                   [fun] G
!
! avalable is :
!                  ![fun] i
! others      :
!                  ![sub] p
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module VCA_DeltaH
  use CodeObject
  use LatticeConfig
  use LaprimaryHusedDatatype
  implicit none


  integer,parameter::Nmax = 200

  type,extends(Object)::VCAdH
    private
    integer::state
    class(LaCon),pointer :: Lconf => null()    ! lattice configruatioin


    !---------------------------------
    integer::Nli
    type(idata)::LocalInter(Nmax)
    !---------------------------------
  contains
   procedure,pass::Initialization
   procedure,pass::StartAppending
   procedure,pass::EndAppending
   procedure,pass::Append
   procedure,pass::AppendOrbitalCDW
  endtype

  private::Initialization,StartAppending,EndAppending,Append
  private::AppendDeltaAF,AppendOrbitalCDW

contains


  subroutine Initialization(self,Lconf)
    implicit none
    class(VCAdH),intent(inout)::self
    class(LaCon),target       ::Lconf
    !---------------------------------------
    call self%SetInitiated(.True.)
    self%Lconf => Lconf
    self%state = 1

  endsubroutine


  subroutine StartAppending(self)
    implicit none
    class(VCAdH),intent(inout)::self
    !---------------------------------------
    self%state = 2
    self%Nli   = 0
  endsubroutine

  subroutine EndAppending(self)
    implicit none
    class(VCAdH),intent(inout)::self
    !---------------------------------------
    if (self%state==2)then
      self%state = 3
    else
      write(self%getprint(),*)"ERROR: StartAppending is not called while trying to call EndAppending"
      stop
    endif
  endsubroutine

  subroutine Append(self,data)
    implicit none
    class(VCAdH),intent(inout)::self
    class(idata),intent(in)::data
    !---------------------------------------
    if (self%state==2)then
      self%Nli = self%Nli + 1
      if ( self%Nli .gt. Nmax )then
        write(self%getprint(),*)"ERROR: Nmax in VCAdH is too small";stop
      endif
      self%LocalInter(self%Nli) = data
    else
      write(self%getprint(),*)"ERROR: Append is not allowed, StartAppending is not calleld."
    endif
  endsubroutine


  subroutine AppendDeltaAF(self,B,Qm,discription)
    implicit none
    class(VCAdH),intent(inout)::self
    real*8,intent(in)::B,Qm(3)
    character(DiscLen),intent(in)::discription
    !---------------------------------------
    type(idata)::I
    integer::jc,spin
    complex*16::expiqr
    Do jc = 1 , self%Lconf%GetNs()
       expiqr = Zexp( sum( self%Lconf%GetSiteRealP(jc) * Qm ) * (0._8,1._8) )
       I%Disc=discription
       I%Itype="SpinOnSite"
       I%Para(1:2) = jc
       do spin = 0 , 1
          i%para(3) = spin
          i%V = expiqr * (-1)**spin * B
          call self%append(I)
       enddo
    enddo
  endsubroutine


  Subroutine AppendOrbitalCDW(self,E,Qm,oi,oj,discription)
    implicit none
    class(VCAdH),intent(inout)::self
    real*8,intent(in)::E,Qm(3)
    integer,intent(in)::Oi,Oj
    character(DiscLen),intent(in)::discription
    !---------------------------------------
    integer::jc,o
    type(idata)::I
    complex*16::expiqr
    do jc = 1 , self%Lconf%GetNs()
       expiqr = Zexp( sum( self%Lconf%GetSiteRealP(jc) * Qm ) * (0._8,1._8) )
       I%Disc=discription
       I%Itype="OnSite"
       I%Para(1:2) = jc
       o =  self%Lconf%GetOrbitIndex(jc)
       if ( o == oi )then
         I%V =   E * expiqr
         CALL SELF%Append(I)
       elseif (o == oj)then
         I%v = - E * expiqr
         CALL SELF%Append(I)
       endif
    enddo
  endsubroutine






endmodule
