

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  :  MODULE
! NAME  :  CEsolver
! OBJECT:  TYPE(CES)
! USED  :  CodeObject,LA_GCESpace,ED_WholeSpace,functionalsubs
! DATE  :  2017-12-29
! AUTHOR:  hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            Canonical Emsenbel Sytem.   Using ED and LA.
!
! STANDARD:
!            *CALL Initialization
!            *call SynchronizeWithHamiltonian()
!
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  [sub] Initialization(Svpara,Ta,Ha,print_,show_)
!                        class(SolverPara),intent(in) :: Svpara
!                        class(table),target          :: Ta
!                        class(Ham),target            :: Ha
!                        integer,intent(in),optional::print_,show_
!
!
! avalable gets:
!                   [fun] GetNs()
!                         integer::GetNs system site.
!
!                   [fun] GetSolverType()
!                         character(*)::GetSolverType
!
!                   [fun] GetTemperature()
!
!
!                   [fun] GerEDpara()
!                         return  type(SolverPara)
!
!                   [fun] GetSymmetry()
!
!                   [fun] GetEdPointer()
!                         return TYPE(ED_GCE) pointer
!
!                   [fun] GetLaPointer()
!                         return TYPE(LA_GCE) pointer
!
!                   [fun] GetRZ()
!                         "ED" only. The reduced partition function.
!
!                   [fun] GetDe()
!                         "LA" only. Degeneracy
!
!                   [fun] GetEg()
!                         Ground state energy
!
!                   [fun] GetOperatorProduct(A)
!                         TYPE(FermOper)::A
!                         COMPLEX*16::GetOperatorProduct
!                         return
!                               <A> = Tr( e^-{\beta*H} A ) / Tr( e^-{\beta*H} )
!
!                   [fun] GetGrandPotential()
!                         real*8
!
! avalable is :
!                  ![fun] i
! others      :
!                  [sub] SynchronizeWithHamiltonian()
!                        diagonalization (if needed)
!
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


Module CEsolver
  use LA_GCESpace
  use ED_WholeSpace
  use CodeObject
  implicit none



  type::SolverPara
    character(2):: SvType      = "NO"    ! = "ED" / "LA" /
    integer     :: job         =  0      !
    logical     :: IsReal      = .false. !
    integer     :: symmetry    =  0      !
    real*8      :: Temperature =  0      !
    !--------------------------
    ! ED use
    !
    !--------------------------
    !  Lanczos use
    real*8  :: pre      = 1.e-13
    real*8  :: bzero    = 0.000001
    integer :: M        = 30
    logical :: oth      = .true.
    real*8  :: DegPre   = 1.e-6
    !***************************
    ! Dictionary for SvType:
    character(2):: EDGCE  = "ED"
    character(2):: LAGCE  = "LA"
    !----------------
  endtype



  TYPE,extends(Object)::CES

    class(table),pointer :: ta
    type(SolverPara)     :: Spara
    TYPE(LA_GCE)         :: LA
    TYPE(ED_GCE)         :: ED

  contains
    procedure,pass::Initialization
    procedure,pass::GetEdPointer
    procedure,pass::GetLaPointer
    procedure,pass::SynchronizeWithHamiltonian

    !-----------------------------------------
    !   ED FUNCTION
    procedure,pass::GetRZ

    !-----------------------------------------
    !   LA function
    procedure,pass::GetDe
    !-----------------------------------------
    !  common
    procedure,pass::GetSolverType
    procedure,pass::GetEg
    procedure,pass::GetOperatorProduct
    procedure,pass::GetNs
    procedure,pass::GetGrandPotential
    procedure,pass::GetSymmetry
    procedure,pass::GerEDpara
    procedure,pass::GetTemperature
  endtype



  private::Initialization
  private::GetEdPointer
  private::GetLaPointer
  private::InterGetEDPointer,InterGetLAPointer
  private::SynchronizeWithHamiltonian
  private::CheckTypeOrStop

  private::GetRZ,GetDe
  private::GetOperatorProduct
  private::GetSolverType
  private::GetNs
  private::GetGrandPotential
  private::GetSymmetry
  private::GerEDpara
  private::GetTemperature
contains





    subroutine Initialization(self,Svpara,Ta,Ha,print_,show_)
      use functionalsubs
      implicit none
      class(CES),intent(inout)     :: self
      class(SolverPara),intent(in) :: Svpara
      class(table),target          :: Ta
      class(Ham),target            :: Ha

      integer,intent(in),optional::print_,show_
      !-------------------------------------------
      TYPE(funcsubs)::f
      call self%SetInitiated(.true.)


      if (present(PRINT_))  call self%setprint(print_)
      if (present(SHOW_))   call self%setshow(show_)!self%show   = SHOW_

      self%Spara = svpara
      !-------------------------------------------------

      self%Spara%SvType = f%get_string_upper(self%Spara%SvType)
      self%Ta => ta
      ! select case(self%Spara%SvType)
      ! case("ED")
      !   call self%ed%Initialization(T=self%Spara%Temperature,Ta=Ta,CH=Ha,IsReal=self%spara%IsReal &
      !      ,print_=self%getprint())
      ! case("LA")
      !   call self%la%Initialization(Ta=Ta,H=Ha,IsReal_=self%spara%IsReal,PRINT_=self%getprint(),&
      !        SHOW_=self%getshow(),Pre_=self%spara%pre,DegPre_=self%spara%DegPre,Bzero_=self%spara%bzero,&
      !        M_=self%spara%m,oth_=self%spara%oth)
      ! case default
      !   write(self%getprint(),*)"ERROR: Unknow solver type:",self%Spara%SvType;stop
      ! endselect
      if (self%Spara%SvType==self%Spara%EDGCE)then
         call self%ed%Initialization(T=self%Spara%Temperature,Ta=Ta,CH=Ha,IsReal=self%spara%IsReal &
           ,print_=self%getprint())
         goto 101 ! out case
      endif
      if (self%Spara%SvType==self%Spara%LAGCE)then
        call self%la%Initialization(Ta=Ta,H=Ha,IsReal_=self%spara%IsReal,PRINT_=self%getprint(),&
             SHOW_=self%getshow(),Pre_=self%spara%pre,DegPre_=self%spara%DegPre,Bzero_=self%spara%bzero,&
             M_=self%spara%m,oth_=self%spara%oth)
        GOTO 101
      endif

      write(self%getprint(),*)"ERROR: Unknow solver type:",self%Spara%SvType;stop
  101 continue

    endsubroutine


  function InterGetEDPointer(ED)  result(p)
    implicit none
    class(ED_GCE),target::ED
    class(ED_GCE),POINTER::P
    !--------------------------
    p => ed
  endfunction
  function InterGetLAPointer(LA)  result(p)
    implicit none
    class(LA_GCE),target::LA
    class(LA_GCE),POINTER::P
    !--------------------------
    p => LA
  endfunction


  function GetEdPointer(self) result(r)
    implicit none
    class(CES),intent(inout)   :: self
    class(ED_GCE),pointer::r
    !-------------------------------------
    call SELF%CheckInitiatedOrStop()
    if (self%Spara%SvType==self%Spara%EDGCE)then
      r => InterGetEDPointer(self%ed)!self%ed
    else
      write(self%getprint(),*)"ERROR: Rectent solver is not ED type";stop
    endif
  endfunction

  function GetLaPointer(self) result(r)
    implicit none
    class(CES),intent(inout)   :: self
    class(LA_GCE),pointer::r
    !-------------------------------------
    call SELF%CheckInitiatedOrStop()
    if (self%Spara%SvType==self%Spara%laGCE)then
      r => InterGetLAPointer(SELF%LA) !self%la
    else
      write(self%getprint(),*)"ERROR: Rectent solver is not LA type";stop
    endif
  endfunction


  subroutine SynchronizeWithHamiltonian(self)
    implicit none
    class(CES),intent(inout)   :: self
    !-------------------------------------
    ! select case(self%Spara%SvType)
    ! case("ED")
    !   call self%ed%SynchronizeWithHamiltonian()
    ! case("LA")
    !   call self%la%SynchronizeWithHamiltonian()
    ! end select

    if(self%Spara%SvType==self%Spara%EDGCE)then
      call self%ed%SynchronizeWithHamiltonian()
      GOTO 101
    endif
    if(self%Spara%SvType==self%Spara%LAGCE)then
      call self%la%SynchronizeWithHamiltonian()
      GOTO 101
    endif

    WRITE(SELF%GETPRINT(),*)"Type is not defined in SynchronizeWithHamiltonian@CES"
    stop
101 continue
  endsubroutine



 !--------------------------
 subroutine CheckTypeOrStop(self,Tp)
   implicit none
   class(CES),intent(inout)   :: self
   character(2),intent(in)    :: Tp
   !-------------------------------------
   call self%CheckInitiatedOrStop()
   if (self%Spara%SvType==Tp)then
   else
      write(self%getprint(),*)"ERROR: Solver type is not ",TP;stop
   endif
 endsubroutine
 !--------------------------


  real*8 function GetRZ(self)
    implicit none
    class(CES),intent(inout)   :: self
    !-------------------------------------
    ! call CheckTypeOrStop(self,"ED")
    ! GetRZ = self%ed%GetRz()
    !-----------------------------------------------
    if (self%Spara%SvType==self%Spara%EDGCE)then
        GetRZ = self%ed%GetRz()
        goto 101
    endif
    !-----------------------------------------------
    write(self%getprint(),*)"Unknow ED type in GetRZ@CESolver";stop
101 continue
  endfunction

  integer FUNCTION GetDe(self)
    implicit none
    class(CES),intent(inout)   :: self
    !-------------------------------------
    ! call CheckTypeOrStop(self,"LA")
    ! GetDe = SELF%LA%GetDe()
    if (self%Spara%SvType==self%Spara%LAGCE)then
        GetDe = SELF%LA%GetDe()
        goto 101
    endif
    !-----------------------------------------------
    write(self%getprint(),*)"Unknow ED type in GetDe@CESolver";stop
101 continue
  endfunction

  real*8 function GetEg(self)
    implicit none
    class(CES),intent(inout)   :: self
    !-------------------------------------
    ! select case(self%Spara%SvType)
    ! case("ED")
    !   GetEg = self%ed%get_Eg()
    ! case("LA")
    !   GetEg = self%LA%GetEg()
    ! endselect
    if (self%Spara%SvType==self%Spara%EDGCE)then
        GetEg = self%ed%get_Eg()
        goto 101
    endif
    !--------------------------------------------
    if (self%Spara%SvType==self%Spara%LAGCE)then
        GetEg = self%LA%GetEg()
        goto 101
    endif
    !--------------------------------------------
    write(self%getprint(),*)"Unknow ED type in GetEg@CESolver";stop
101 continue
  endfunction


  complex*16 function GetOperatorProduct(self,opt)
    implicit none
    class(CES),intent(inout)::self
    class(FermOper),intent(inout)::opt
    !-------------------------------------
    call self%CheckInitiatedOrStop()
    call self%SynchronizeWithHamiltonian()
    ! select case(self%Spara%SvType)
    ! case("ED")
    !   GetOperatorProduct = self%ed%get_trace_value(opt)
    ! case("LA")
    !   GetOperatorProduct = self%la%GetOperateProduct(opt)
    ! endselect
    if (self%Spara%SvType==self%Spara%EDGCE)then
        GetOperatorProduct = self%ed%get_trace_value(opt)
        goto 101
    endif
    !--------------------------------------------
    if (self%Spara%SvType==self%Spara%LAGCE)then
        GetOperatorProduct = self%la%GetOperateProduct(opt)
        goto 101
    endif
    !--------------------------------------------
    write(self%getprint(),*)"Unknow ED type in GetOperatorProduct@CESolver";stop
101 continue
  endfunction



  function GetSolverType(self) result(r)
    implicit none
    class(CES),intent(inout)::self
    character(len=len(self%Spara%SvType))::r
    !-------------------------------------
    call self%CheckInitiatedOrStop()
    r = self%Spara%SvType
  endfunction


  integer function GetSymmetry(self)
    implicit none
    class(CES),intent(inout)::self
    !-------------------------------------
    GetSymmetry = self%Spara%symmetry
  endfunction

  integer function GetNs(self)
    implicit none
    class(CES),intent(inout)::self
    !-------------------------------------
    GetNs = self%ta%get_ns()
  endfunction

  ! real*8 function GetGrandPotential(self)
  !   implicit none
  !   class(CES),intent(inout)::self
  !   !-------------------------------------
  !   select case(self%Spara%SvType)
  !   case("ED")
  !     GetGrandPotential = self%ed%get_Eg() &
  !                       - self%Spara%Temperature*dlog( self%ed%GetRz() ) !;write(*,*)self%ed%get_Eg(),self%Spara%Temperature,self%ed%GetRz(),999;stop
  !   case("LA")
  !     GetGrandPotential = self%la%GetEg()
  !   endselect
  ! endfunction
  real*8 function GetGrandPotential(self)
    implicit none
    class(CES),intent(inout)::self
    !-------------------------------------
    if (self%Spara%SvType==self%Spara%EDGCE)then
        GetGrandPotential = self%ed%get_Eg() &
                            - self%Spara%Temperature*dlog( self%ed%GetRz() )
        goto 101
    endif
    !--------------------------------------------
    if (self%Spara%SvType==self%Spara%LAGCE)then
        GetGrandPotential = self%la%GetEg()
        goto 101
    endif
    !--------------------------------------------
    write(self%getprint(),*)"Unknow ED type in GetGrandPotential@CESolver";stop
101 continue
  endfunction



  type(SolverPara) function GerEDpara(self)
    implicit none
    class(CES),intent(inout)     :: self
    !-------------------------------------
    GerEDpara = self%Spara
  endfunction



  real*8 function GetTemperature(self)
    use functionalsubs
    implicit none
    class(CES),intent(in) :: self
    !--------------------------------------
    GetTemperature = self%Spara%Temperature
  endfunction





endmodule
