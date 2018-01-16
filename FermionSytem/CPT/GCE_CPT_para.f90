


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : GCE_CPT_para
! OBJECT: TYPE(GCECPTpara)
! USED  : CodeObject , fermion_table , LaLatticeH , CEsolver , FermionHamiltonian , CE_Green
! DATE  : 2017-12-27
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            Used for CPT calcuation. This type used to save a CPT problem.
!          When we do some CPT calculation, this type would be one of the input.
!
! STANDARD:
!            *CALL Initialization( CPTH , EDpara )
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                   [sub] Initialization(edpara,Gpara,cpth)
!                          class(SolverPara),intent(in)::edpara
!                          class(GreenPara),intent(in) :: Gpara
!                          class(lh),intent(in) ::cpth
!
!
!
! avalable gets:
!
!                   [fun] GetSolverPointer()
!                         Class(CES),pointer::GetSolverPointer
!
!                   [fun] GetGp()
!                         return Type(GreenPara)
!
!                   [fun] GetCPTHpointer()
!                         class(lh)::CPTH
!
!                   [fun] GetLatticePointer()
!                         return class(LaCon)
!
!                   [fun] GetTemperature()
!
!
!
! avalable is :
!                  ![fun] i
! others      :
!                   [sub] ReportPara(wtp,mode)
!                         report parameters
!
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module GCE_CPT_para
  use CodeObject
  use fermion_table
  use LaLatticeH
  use CEsolver
  use CE_Green
  use FermionHamiltonian
  implicit none


  type,extends(object)::GCECPTpara
    private
    TYPE(LH)        :: CPTH
    TYPE(table)     :: Ta
    TYPE(Ham)       :: SolverH
    TYPE(CES)       :: solver
    Type(GreenPara) :: GP

  contains
    procedure,pass::Initialization
    procedure,pass::GetSolverPointer
    procedure,pass::ReportPara
    procedure,pass::GetCPTHpointer
    procedure,pass::GetLatticePointer
    procedure,pass::GetTemperature
    procedure,pass::GetGp

  endtype

  private::Initialization
  private::GetSolverPointer
  private::ReportPara
  private::GetCPTHpointer
  private::GetLatticePointer
  private::GetTemperature
  private::GetGp



contains


  subroutine Initialization(self,edpara,Gpara,cpth)
    implicit none
    class(GCECPTpara),intent(inout) :: self
    class(SolverPara),intent(in)    :: edpara
    class(GreenPara),intent(in)     :: Gpara
    class(lh),intent(in)            :: cpth
    !----------------------------------------
    call self%SetInitiated(.True.)

    self%gp = Gpara
    !--------------------------------------------------------------------
    ! set CPTH
    self%CPTH    = cpth
    !--------------------------------------------------------------------
    ! set H
    call self%SolverH%Initialization(  self%CPTH%GetNs() )
    call self%SolverH%StartAppendingInteraction()
    call self%CPTH%AppendLocalDataToHam(self%SolverH)
    call self%SolverH%EndAppendingInteraction()
    !--------------------------------------------------------------------
    ! set Ta
    Call self%Ta%Initialization( self%CPTH%GetNs(),  EDpara%symmetry  )
    !--------------------------------------------------------------------
    ! set solver
    Call self%solver%Initialization(Svpara=EDpara,Ta=self%Ta,Ha=self%SolverH,&
                                    print_=self%getprint(),show_=self%getshow())
    !--------------------------------------------------------------------

  endsubroutine


  function GetSolverPointer(self) result(r)
    implicit none
    class(GCECPTpara),intent(inout)::self
    class(CES),pointer::r
    !--------------------------------------
    Call self%CheckInitiatedOrStop()
    r => maketypetopointer(self%solver)
  contains
    function maketypetopointer(ed) result(r)
      class(CES),target::ed
      class(CES),pointer::r
      !-------------------------
      r => ed
    endfunction
  endfunction


  subroutine ReportPara(self,wtp,mode)
    implicit none
    class(GCECPTpara),intent(inout)::self
    integer,intent(in)::wtp,mode
    !--------------------------------------
    call self%CPTH%Report(wtp,mode)
  endsubroutine

  function GetCPTHpointer(self) result(r)
    implicit none
    class(GCECPTpara),intent(in)::self
    class(lh),pointer::r
    !---------------------------------------
    call self%CheckInitiatedOrStop()
    r => getp(self%CPTH)
  contains
    function getp(x) result(y)
      implicit none
      class(lh),target::x
      class(lh),pointer::y
      !--------------------------
      y => x
    endfunction
  endfunction



    function GetLatticePointer(self) result(r)
      implicit none
      class(GCECPTpara),intent(in)::self
      class(LaCon),pointer::r
      !----------------------------------------
      call self%CheckInitiatedOrStop()
      r => self%CPTH%GetLatticeConfigPointer()
    endfunction


    real*8 function GetTemperature(self)
      implicit none
      class(GCECPTpara),intent(in)::self
      !---------------------------------------
      GetTemperature = self%solver%GetTemperature()
    endfunction

    type(GreenPara) FUNCTION GetGp(self)
      implicit none
      class(GCECPTpara),intent(inout) :: self
      !----------------------------------------
      GetGp = self%GP
    endfunction



endmodule
