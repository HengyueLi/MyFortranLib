

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
!                  [sub] Initialization(Tp,Te,Ta,Ha,isreal_,print_,show_,Pre_,DegPre_,Bzero_,M_,oth_)
!                        Character(2),intent(in)       :: Tp           ! = "ed" or "LA"
!                        Real*8,intent(in)             :: Te           ! temperature
!                        class(table),target           :: Ta           ! table
!                        class(Ham),target             :: Ha           ! Hamiltonian
!                        logical,intent(in),optional   :: IsReal_      ! real problem or not
!                        integer,intent(in),optional   :: print_,show_ !
!                        !----------lanczos options---------------------------
!                        real*8,intent(in),optional    :: Pre_,DegPre_,Bzero_
!                        integer,intent(in),optional   :: M_
!                        logical,intent(in),optional   :: oth_
!-------------------------------------------
!
!
!
!
! avalable gets:
!                   [fun] GetNs()
!                         integer::GetNs system site.
!
!                   [fun] GetSolverType()
!                         character(2)::GetSolverType
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



  TYPE,extends(Object)::CES

    Character(2)::Tp   ! = "LA" or "ED" to choose solver

    TYPE(LA_GCE)::LA
    TYPE(ED_GCE)::ED
    class(table),pointer ::ta
    !-------------------------------------
    !  default value for solvers
    logical::Isreal
    ! la
    real*8 ::pre      = 1.e-14
    real*8 ::bzero    = 0.000001
    integer::M        = 30
    logical::oth      = .true.
    real*8 ::DegPre   = 1.e-6


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
contains


  subroutine Initialization(self,Tp,Te,Ta,Ha,isreal_,print_,show_    ,Pre_,DegPre_,Bzero_,M_,oth_)
    use functionalsubs
    implicit none
    class(CES),intent(inout)   :: self
    Character(2),intent(in)    :: Tp
    Real*8,intent(in)          :: Te
    class(table),target        :: Ta
    class(Ham),target          :: Ha
    logical,intent(in),optional:: IsReal_
    integer,intent(in),optional::print_,show_
    real*8,intent(in),optional     :: Pre_,DegPre_,Bzero_
    integer,intent(in),optional    :: M_
    logical,intent(in),optional    :: oth_
    !-------------------------------------------
    TYPE(funcsubs)::f
    call self%SetInitiated(.true.)

    if(present(IsReal_))  self%isreal = IsReal_
    if (present(PRINT_))  call self%setprint(print_)
    if (present(SHOW_))   call self%setshow(show_)!self%show   = SHOW_
    if (present(DegPre_)) self%DegPre = DegPre_
    if (present(Pre_))    self%Pre    = Pre_
    if (present(Bzero_))  self%bzero  = Bzero_
    if (present(M_))      self%M      = M_
    if (present(oth_))    self%oth    = oth_
    !-------------------------------------------------

    self%Tp = f%get_string_upper(Tp)
    self%Ta => ta
    select case(self%Tp)
    case("ED")
      call self%ed%Initialization(T=Te,Ta=Ta,CH=Ha,IsReal=self%isreal,print_=self%getprint())
    case("LA")
      call self%la%Initialization(Ta=Ta,H=Ha,IsReal_=self%isreal,PRINT_=self%getprint(),&
           SHOW_=self%getshow(),Pre_=self%pre,DegPre_=self%DegPre,Bzero_=self%bzero,M_=self%m,oth_=self%oth)
    case default
      write(self%getprint(),*)"ERROR: Unknow solver type:",Tp;stop
    endselect

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
    if (self%Tp=="ED")then
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
    if (self%Tp=="LA")then
      r => InterGetLAPointer(SELF%LA) !self%la
    else
      write(self%getprint(),*)"ERROR: Rectent solver is not LA type";stop
    endif
  endfunction


  subroutine SynchronizeWithHamiltonian(self)
    implicit none
    class(CES),intent(inout)   :: self
    !-------------------------------------
    select case(self%tp)
    case("ED")
      call self%ed%SynchronizeWithHamiltonian()
    case("LA")
      call self%la%SynchronizeWithHamiltonian()
    end select
  endsubroutine



 !--------------------------
 subroutine CheckTypeOrStop(self,Tp)
   implicit none
   class(CES),intent(inout)   :: self
   character(2),intent(in)    :: Tp
   !-------------------------------------
   call self%CheckInitiatedOrStop()
   if (self%Tp==Tp)then
   else
      write(self%getprint(),*)"ERROR: Solver type is not ",TP;stop
   endif
 endsubroutine
 !--------------------------


  real*8 function GetRZ(self)
    implicit none
    class(CES),intent(inout)   :: self
    !-------------------------------------
    call CheckTypeOrStop(self,"ED")
    GetRZ = self%ed%GetRz()
  endfunction

  integer FUNCTION GetDe(self)
    implicit none
    class(CES),intent(inout)   :: self
    !-------------------------------------
    call CheckTypeOrStop(self,"LA")
    GetDe = SELF%LA%GetDe()
  endfunction

  real*8 function GetEg(self)
    implicit none
    class(CES),intent(inout)   :: self
    !-------------------------------------
    select case(self%tp)
    case("ED")
      GetEg = self%ed%get_Eg()
    case("LA")
      GetEg = self%LA%GetEg()
    endselect
  endfunction


  complex*16 function GetOperatorProduct(self,opt)
    implicit none
    class(CES),intent(inout)::self
    class(FermOper),intent(inout)::opt
    !-------------------------------------
    call self%CheckInitiatedOrStop()
    call self%SynchronizeWithHamiltonian()
    select case(self%tp)
    case("ED")
      GetOperatorProduct = self%ed%get_trace_value(opt)
    case("LA")
      GetOperatorProduct = self%la%GetOperateProduct(opt)
    endselect
  endfunction



  character(2) function GetSolverType(self)
    implicit none
    class(CES),intent(inout)::self
    !-------------------------------------
    call self%CheckInitiatedOrStop()
    GetSolverType = self%tp
  endfunction



  integer function GetNs(self)
    implicit none
    class(CES),intent(inout)::self
    !-------------------------------------
    GetNs = self%ta%get_ns()
  endfunction






endmodule
