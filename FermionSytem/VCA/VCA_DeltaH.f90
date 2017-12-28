

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : VCA_DeltaH
! OBJECT: TYPE(VCAdH)
! USED  : CodeObject,LatticeConfig,LaPrimaryH,functionalsubs
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
!                   [sub] Initialization(Lconf)
!                         class(LaCon),target::Lconf
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
!                         Delta =   B * sum_{r,spin} (-1)^spin Exp(i Qm * r )
!
!                   [sub] AppendOrbitalCDW(E,Qm,oi,oj,discription)
!                         real*8,intent(in)::E,Qm(3)
!                         integer,intent(in)::Oi,Oj
!                         character(DiscLen),intent(in)::discription   DiscLen=32?
!
!                         Delta = E * sum_{r} Exp(i Qm * r ) * (n_oi - n _oj)
!                             notice the order of oi,oj
!
!                   [sub] AppendOrbitalESDW(E,Qm,oi,oj,discription)
!                         real*8,intent(in)::E,Qm(3)
!                         integer,intent(in)::Oi,Oj
!                         character(DiscLen),intent(in)::discription   DiscLen=32?
!
!                         Delta = E * sum_r Exp(i Qm * r ) * (-1)^spin * (  c^+_{r,oi} c_{r,oj} +h.c.  )
!
!                   [sub] SetValueByDiscription(discription,v)
!                         character(DiscLen),intent(in)::discription   DiscLen=32?
!                         complex*16,intent(in)::v
!
! avalable gets:
!                   [fun] GetDelatMatrix(spini,spinj)
!                         complex*16::GetDelatMatrix(ns,ns)
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
    integer::ns
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
   procedure,pass::AppendDeltaAF
   procedure,pass::AppendOrbitalCDW
   procedure,pass::AppendOrbitalESDW
   procedure,pass::GetDelatMatrix
   procedure,pass::SetValueByDiscription
  endtype

  private::Initialization,StartAppending,EndAppending,Append
  private::AppendDeltaAF,AppendOrbitalCDW,AppendOrbitalESDW
  private::GetDelatMatrix
  private::SetValueByDiscription

contains


  subroutine Initialization(self,Lconf)
    implicit none
    class(VCAdH),intent(inout)::self
    class(LaCon),target       ::Lconf
    !---------------------------------------
    call self%SetInitiated(.True.)
    self%Lconf => Lconf
    self%state = 1
    self%ns    = self%Lconf%GetNs()
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


  subroutine AppendOrbitalESDW(self,E,Qm,oi,oj,discription)
    implicit none
    class(VCAdH),intent(inout)::self
    real*8,intent(in)::E,Qm(3)
    integer,intent(in)::Oi,Oj
    character(DiscLen),intent(in)::discription
    !---------------------------------------
    integer::jci,jcj,o1,o2,pi(3),pj(3),spin
    type(idata)::I
    complex*16::expiqr
    do jci = 1 , self%Lconf%GetNs()
       do jcj = jci + 1 , self%Lconf%GetNs()
          pi = self%Lconf%GetSiteP(jci)
          pj = self%Lconf%GetSiteP(jcj)
          o1 = self%Lconf%GetOrbitIndex(jci)
          o2 = self%Lconf%GetOrbitIndex(jcj)
          if ( sum( abs(pi-pj) )==0 )then
             if (  ( (o1==oi) .and. (o2==oj)  ) .or. ( (o2==oi) .and. (o1==oj)  ) )then
               !-----------------------------
               !  site jci,jcj
               expiqr = Zexp( sum( self%Lconf%GetSiteRealP(jci) * Qm ) * (0._8,1._8) )
               I%Disc=discription
               I%Itype="SpinHopping"
               I%Para(1) = jci   ;  i%Para(2) = jcj
               do spin = 0 , 1
                  I%para(3) = spin
                  I%v = (-1)**spin * expiqr *  E
                  CALL SELF%Append(I)
               enddo
             endif
          endif
       enddo
    enddo
  endsubroutine


  function GetDelatMatrix(self,spini,spinj) result(r)
    implicit none
    class(VCAdH),intent(inout)::self
    integer,intent(in)::spini,spinj
    complex*16::r(self%ns,self%ns)
    !---------------------------------------
    integer::jc
    call self%CheckInitiatedorStop()
    if (self%state==3)then
      r = (0._8,0._8)
      do jc = 1 , self%Nli
         call self%LocalInter(jc)%SetIntoMatix(self%ns,r,spini,spinj,.true.)   !;write(*,*)r;pause
      enddo
    else
      write(self%getprint(),*)"ERROR: EndAppend should be called before GetDelatMatrix";stop
    endif
  endfunction

  subroutine SetValueByDiscription(self,discription,v)
    use functionalsubs
    implicit none
    class(VCAdH),intent(inout)::self
    character(DiscLen),intent(in)::discription
    complex*16,intent(in)::v
    !---------------------------------------
    integer::jc
    character(DiscLen)::DisIn,DisCheck
    TYPE(funcsubs)::f

    call self%CheckInitiatedOrStop()


    DisIn = adjustl(trim(f%get_string_upper(discription)))
    do jc = 1 , self%Nli
      DisCheck = self%LocalInter(jc)%Disc
      DisCheck = adjustl(trim(f%get_string_upper(DisCheck)))
      if (DisIn.eq.DisCheck)then
         self%LocalInter(jc)%v = v
      endif
    enddo
  endsubroutine






endmodule
