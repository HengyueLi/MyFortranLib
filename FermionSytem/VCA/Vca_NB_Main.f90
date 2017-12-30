
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : VCA_NB_MAIN
! OBJECT: TYPE(VCANB)
! USED  : CodeObject , LaPrimaryH  , LatticeConfig  ,VCA_DeltaH , LaLatticeH ,fermion_table
! DATE  : 2017-12-30
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            Vational Cluster Approach without bath sites.
!
!            The Major input is CPT_H. The output is DeltaH. We fixed DeltaH by variational process.
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
!                   [sub] I
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



module VCA_NB_MAIN
  use CodeObject
  use LaPrimaryH
  use LatticeConfig
  use VCA_DeltaH
  use LaLatticeH
  use fermion_table
  implicit none


  type,extends(object)::VCANB
    private

    !--------------------------------------------------------------------
    ! temperature
      real*8      :: Tem
    ! symmetry
      integer     :: Symm
    ! H is real or not
      logical     :: IsReal
    !--------------------------------------------------------------------
    ! saving a Hamitonian in a primary cell. No geometry is contained.
      TYPE(PH)    :: PrimaryH
    !--------------------------------------------------------------------
    ! the lattice configuration of the system.
      TYPE(LaCon) :: LattConf
    !--------------------------------------------------------------------
    ! variational meanfiled
      TYPE(VCAdH) :: DeltaH
    !--------------------------------------------------------------------



    !====================================================================
    !  inner usage
    !====================================================================
    ! Frequently be used by solver. Thus we save the value.
      Type(table) :: SolverTa
    ! Cpt H
      TYPE(LH) :: CPTH


  contains
    !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    !  override

    !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&




    procedure,pass::Initialization

  endtype


  private::Initialization




contains





  subroutine Initialization(self,T,Symm,IsReal,PrimaryH,LattConf,DeltaH)
    implicit none
    class(VCANB),intent(out):: self
    real*8,intent(in)       :: T
    integer,intent(in)      :: Symm
    logical,intent(in)      :: IsReal
    class(PH),intent(in)    :: PrimaryH
    class(LaCon),intent(in) :: LattConf
    class(VCAdH),intent(in) :: DeltaH
  !--------------------------------------------------------------------

    !------------------------------------
    call self%SetInitiated(.true.)
    !------------------------------------
    self%Tem    = T
    self%symm   = symm
    self%IsReal = IsReal
    ! Check Initiated
    if (.not.PrimaryH%IsInitiated() )then
       write(self%getprint(),*)"ERROR: In vca input, PrimaryH is not initiated yet";stop
    endif
    if (.not.LattConf%IsInitiated() )then
       write(self%getprint(),*)"ERROR: In vca input, LattConf is not initiated yet";stop
    endif
    if (.not.DeltaH%IsInitiated() )then
       write(self%getprint(),*)"ERROR: In vca input, DeltaH is not initiated yet";stop
    endif
    self%PrimaryH = PrimaryH
    self%LattConf = LattConf
    self%DeltaH   = DeltaH
    !------------------------------------
    call self%CPTH%Initialization( self%PrimaryH,self%LattConf,self%getprint(),self%getshow())
    !-----------------------------------
    CALL SELF%SolverTa%Initialization( SELF%LattConf%GetNs() , self%symm ,SELF%GETPRINT()) 



  endsubroutine











endmodule















!
! subroutine SetPrimaryCellHamiltonian(self,PrimaryH)
!   implicit none
!   class(VCANB),intent(inout)::self
!   class(PH),intent(inout)   ::PrimaryH
!   !------------------------------------
!   write(self%getprint(),*)"ERROR: SetPrimaryCellHamiltonian is not defined yet."
!   write(self%getprint(),*)"-------------------------------------------------------------------"
!   write(self%getprint(),*)"use Subroutine PrimaryH%Append(Disc,Pos,Itype,Ipara,Ivalu) where"
!   write(self%getprint(),*)"character(DiscLen)  :: Disc          ( DiscLen = 32?  )"
!   write(self%getprint(),*)"integer             :: pos(3)         denote cluster position"
!   write(self%getprint(),*)"character(ITypeLen) :: Itype         ( ITypeLen = 16? )"
!   write(self%getprint(),*)"integer             :: Ipara(8)       1-basis for site index.  0,1 for spin up and down."
!   write(self%getprint(),*)"complex*16          :: Ivalu"
!   write(self%getprint(),*)"Itype can be checked in module FermionHamiltonian."
!   write(self%getprint(),*)"And only for the terms that the first two elements in para(8) represent sites can be used."
!   stop
! endsubroutine
!
