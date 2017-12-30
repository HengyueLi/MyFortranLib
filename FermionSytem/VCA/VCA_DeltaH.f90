

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : VCA_DeltaH
! OBJECT: TYPE(VCAdH)
! USED  : CodeObject,LatticeConfig,LaPrimaryH,functionalsubs,basic_math_functions
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

  !----------------------------------
  !  for idata array
  integer,parameter::Nmax = 200
  !----------------------------------
  !  max number of variational terms
  integer,parameter:: NVmax = 20

  type,private::Vterm
    character(32) :: disc   ! discription
    character(16) :: Vtype  ! check type
    real*8        :: V
    real*8        :: rpara(5)
    integer       :: ipara(5)
    character(3)  :: how    ! = "MIN", "MAX"
  endtype




  type,extends(Object)::VCAdH
    private
    integer::EigenID
    integer::ns
    integer::state
    class(LaCon),pointer :: Lconf => null()    ! lattice configruatioin

    !--------------------------------
    integer    ::Nvar
    type(Vterm)::Var(NVmax)
    !---------------------------------
  contains
   procedure,pass::Initialization
   procedure,pass::StartAppending
   procedure,pass::EndAppending
  !  procedure,pass::Append
   procedure,pass::AppendDeltaAF
   procedure,pass::AppendOrbitalCDW
   procedure,pass::AppendOrbitalESDW
   procedure,pass::GetDelatMatrix
   procedure,pass::SetValueByDiscription
  endtype

  private::Initialization,StartAppending,EndAppending!,Append
  private::AppendDeltaAF,AppendOrbitalCDW,AppendOrbitalESDW
  private::GetDelatMatrix
  private::SetValueByDiscription
  private::CheckNarOverflow,SetToUpperCase
  private::checkNmaxSamll,GetTotalIdataArray
  private::GetIdataArrayFromVarType
  private::GetVarIdByDiscription
  private::SetRadomEigenId

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
    ! self%Nli   = 0

    self%Nvar  = 0
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

  ! subroutine Append(self,data)
  !   implicit none
  !   class(VCAdH),intent(inout)::self
  !   class(idata),intent(in)::data
  !   !---------------------------------------
  !   if (self%state==2)then
  !     self%Nli = self%Nli + 1
  !     if ( self%Nli .gt. Nmax )then
  !       write(self%getprint(),*)"ERROR: Nmax in VCAdH is too small";stop
  !     endif
  !     self%LocalInter(self%Nli) = data
  !   else
  !     write(self%getprint(),*)"ERROR: Append is not allowed, StartAppending is not calleld.";stop
  !   endif
  ! endsubroutine




 subroutine CheckNarOverflow(self)
   implicit none
   class(VCAdH),intent(inout)::self
    !---------------------------------------
    if (self%Nvar.gt.NVmax)then
      write(self%getprint(),*)"ERROR: NVmax=",NVmax," is too small.";stop
    endif
 endsubroutine

 subroutine SetToUpperCase(self,varid)
   use functionalsubs
   implicit none
   class(VCAdH),intent(inout)::self
   integer,intent(in)::varid
   !---------------------------------
   TYPE(funcsubs)::f
   self%var(varid)%Disc = trim(adjustl(f%get_string_upper(self%var(varid)%Disc)))
 endsubroutine

 SUBROUTINE CheckAllowedAppendOrStop(self )
   implicit none
   class(VCAdH),intent(inout)::self
   !-------------------------------------
   if (.not.self%state==2)then
     write(self%getprint(),*)"ERROR: append is not allowed in VCAdH, consider to call startappending"
     stop
   endif
 endsubroutine


  subroutine AppendDeltaAF(self,B,Qm,discription)
    implicit none
    class(VCAdH),intent(inout)::self
    real*8,intent(in)::B,Qm(3)
    character(DiscLen),intent(in)::discription
    !---------------------------------------
    call CheckAllowedAppendOrStop(self )
    self%Nvar = self%Nvar + 1
    call CheckNarOverflow(self)
    self%var(self%Nvar)%vtype      = "AF"
    self%var(self%Nvar)%disc       = discription
    self%var(self%Nvar)%v          = B
    self%var(self%Nvar)%rpara(1:3) = Qm
    self%var(self%Nvar)%how        = "MIN"
    call SetToUpperCase(self,self%Nvar)
  endsubroutine

  Subroutine AppendOrbitalCDW(self,E,Qm,oi,oj,discription)
    implicit none
    class(VCAdH),intent(inout)::self
    real*8,intent(in)::E,Qm(3)
    integer,intent(in)::Oi,Oj
    character(DiscLen),intent(in)::discription
    !---------------------------------------
    call CheckAllowedAppendOrStop(self )
    self%Nvar = self%Nvar + 1
    call CheckNarOverflow(self)
    self%var(self%Nvar)%vtype      = "ORBITALCDW"
    self%var(self%Nvar)%disc       = discription
    self%var(self%Nvar)%v          = E
    self%var(self%Nvar)%rpara(1:3) = Qm
    self%var(self%Nvar)%ipara(1)   = oi
    self%var(self%Nvar)%ipara(2)   = oj
    self%var(self%Nvar)%how        = "MIN"
    call SetToUpperCase(self,self%Nvar)
  endsubroutine

! Append_Va_DeltaAF,Append_Va_OrbitalCDW,Append_Va_OrbitalESDW


subroutine AppendOrbitalESDW(self,E,Qm,oi,oj,discription)
  implicit none
  class(VCAdH),intent(inout)::self
  real*8,intent(in)::E,Qm(3)
  integer,intent(in)::Oi,Oj
  character(DiscLen),intent(in)::discription
  !---------------------------------------
  call CheckAllowedAppendOrStop(self )
  self%Nvar = self%Nvar + 1
  call CheckNarOverflow(self)
  self%var(self%Nvar)%vtype      = "ORBITALESDW"
  self%var(self%Nvar)%disc       = discription
  self%var(self%Nvar)%v          = E
  self%var(self%Nvar)%rpara(1:3) = Qm
  self%var(self%Nvar)%ipara(1)   = oi
  self%var(self%Nvar)%ipara(2)   = oj
  self%var(self%Nvar)%how        = "MIN"
  call SetToUpperCase(self,self%Nvar)
endsubroutine





subroutine checkNmaxSamll(self,i)
  implicit none
  class(VCAdH),intent(inout):: self
  integer,intent(in)::i
  !----------------------------------
  if (i.gt.nmax)then
    write(self%getprint(),*)"ERROR: Nmax in VCAdH is too small";stop
  endif
endsubroutine



SUBROUTINE GetIdataArrayFromVarType(self,Var ,idataarry,n)
  implicit none
  class(VCAdH),intent(inout):: self
  type(Vterm),intent(in)    :: Var
  type(idata),intent(inout) :: idataarry(Nmax)
  integer,intent(out)       :: n
  !----------------------------------------
  complex*16::expiqr
  integer::jc,spin,o,o1,o2,pi(3),pj(3),oj,oi,jci,jcj
  type(idata)::i
  n = 0
  I%Disc= var%disc
  select case(trim(adjustl(Var%Vtype)))
  case("AF")
    Do jc = 1 , self%Lconf%GetNs()
       expiqr = Zexp( sum( self%Lconf%GetSiteRealP(jc) * Var%rpara(1:3) ) * (0._8,1._8) )
       I%Itype="SpinOnSite"
       I%Para(1:2) = jc
       do spin = 0 , 1
          i%para(3) = spin
          i%V = expiqr * (-1)**spin * Var%v
          n = n + 1
          call checkNmaxSamll(self,n)
          idataarry(n) = i
       enddo
    enddo
  case("ORBITALCDW")
    do jc = 1 , self%Lconf%GetNs()
       expiqr = Zexp( sum( self%Lconf%GetSiteRealP(jc) * Var%rpara(1:3) ) * (0._8,1._8) )
       I%Itype="OnSite"
       I%Para(1:2) = jc
       o =  self%Lconf%GetOrbitIndex(jc)
       if ( o == Var%ipara(1) )then
         I%V =   Var%v * expiqr
         n = n + 1
         call checkNmaxSamll(self,n)
         idataarry(n) = i
       elseif (o == Var%ipara(2))then
         I%v = - Var%v * expiqr
         n = n + 1
         call checkNmaxSamll(self,n)
         idataarry(n) = i
       endif
    enddo
  case("ORBITALESDW")
    oi = var%ipara(1)
    oj = var%ipara(2)
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
               expiqr = Zexp( sum( self%Lconf%GetSiteRealP(jci) * Var%rpara(1:3) ) * (0._8,1._8) )
               I%Itype="SpinHopping"
               I%Para(1) = jci   ;  i%Para(2) = jcj
               do spin = 0 , 1
                  I%para(3) = spin
                  I%v = (-1)**spin * expiqr *  var%v
                  n = n + 1
                  call checkNmaxSamll(self,n)
                  idataarry(n) = i
               enddo
             endif
          endif
       enddo
    enddo
  endselect

endsubroutine


subroutine GetTotalIdataArray(self,idataarray,n)
  implicit none
  class(VCAdH),intent(inout):: self
  type(idata),intent(out)::idataarray(Nmax)
  integer,intent(out)::n
  !---------------------------------------------
  integer::jc , NN  ,Ntry
  type(idata)::Inn(Nmax)

  n = 0
  do jc = 1 , self%Nvar
    call GetIdataArrayFromVarType(self, self%Var(jc),Inn,NN )
    Ntry = n + NN
    call checkNmaxSamll(self, Ntry )
    idataarray(n+1:Ntry) = Inn(1:nn)
    n = Ntry
  enddo
endsubroutine


function GetDelatMatrix(self,spini,spinj) result(r)
  implicit none
  class(VCAdH),intent(inout)::self
  integer,intent(in)::spini,spinj
  complex*16::r(self%ns,self%ns)
  !---------------------------------------
  type(idata)::Totidata(nmax)
  integer::jc,n
  call self%CheckInitiatedorStop()
  if (self%state==3)then
     call GetTotalIdataArray(self,Totidata,n)
     r = (0._8,0._8)
     do jc = 1 , n
        call Totidata(jc)%SetIntoMatix(self%ns,r,spini,spinj,.true.)   !;write(*,*)r;pause
     enddo
  else
    write(self%getprint(),*)"ERROR: EndAppend should be called before GetDelatMatrix";stop
  endif
endfunction



  integer function GetVarIdByDiscription(self,Disc)
    use functionalsubs
    implicit none
    class(VCAdH),intent(inout)::self
    character(32),intent(in)::Disc
    !-----------------------------------------------
    TYPE(funcsubs)::f
    character(32)::DiscIn
    integer::jc
    DiscIn = trim(adjustl(f%get_string_upper(Disc)))
    GetVarIdByDiscription = -1
     do jc = 1 , self%Nvar
       if (DiscIn == trim(adjustl(f%get_string_upper(self%var(jc)%Disc))))then
          GetVarIdByDiscription = jc
          goto 999
       endif
      enddo
999 continue
   endfunction


   subroutine SetRadomEigenId(self)
     use basic_math_functions
     implicit none
     class(VCAdH),intent(inout)::self
     !----------------------------------
     TYPE(bmathf)::f
     self%EigenId = f%get_random_int(int(0,4),int(2147483648,4))
   endsubroutine

   subroutine SetValueByDiscription(self,Disc,rV)
     implicit none
     class(VCAdH),intent(inout)::self
     character(32),intent(in)::Disc
     real*8,intent(in)::rV
     !-----------------------------------------------
     integer::jc
     call SetRadomEigenId(self)
     jc = GetVarIdByDiscription(self,Disc)
     if (jc==-1)then
       write(self%getprint(),*)"WARNNING: IN SetValueByDiscription, no matching Disciption is found.&
                               Program Stop";stop
     endif
     self%var(jc)%v = rv
   endsubroutine








endmodule
