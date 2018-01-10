


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : CPT_ARPES
! OBJECT: TYPE(ARPES)
! USED  : CodeObject , GCE_CPT_para , CE_Green , CreateKspace , LatticeConfig
! DATE  : 2017-12-27
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            For input CPTpara, calculate the ARPES spectral.
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
!                   [sub] Initialization(
!
!
!
! avalable gets:
!
!                   [fun] G
! avalable is :
!                  ![fun] i
! others      :
!                   [sub] R
!
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


module CPT_ARPES
  USE CodeObject
  USE GCE_CPT_para
  use CE_Green
  use CreateKspace
  use LatticeConfig
  implicit none


  type,extends(object)::ARPES
    private

    class(GCECPTpara),pointer   :: CPT      => null()
    class(LaCon),pointer        :: lattice  => null()
    type(CEG)                   :: GreenFun
    ! class(GreenPara)            :: Gp
    type(Kspace)                :: pcKs           ! choose prH basis.


  contains
    procedure,pass :: Initialization

  endtype


  private::Initialization


contains


  subroutine Initialization(self  ,CPTpara,GP,nk)
    implicit none
    class(ARPES),intent(out)    :: self
    class(GCECPTpara),target    :: CPTpara
    class(GreenPara),intent(in) :: GP
    integer,intent(in)          :: nk(3)
    !-----------------------------------------
    call self%SetInitiated(.true.)
    self%CPT     => CPTpara
    self%lattice => self%CPT%GetLatticePointer()

    !---------------------------------------------------------------------------
    ! Green function
    call self%GreenFun%Initialization( solver = self%CPT%GetSolverPointer() ,&
                                                GP = GP,PRINT_=self%getprint() )
    !---------------------------------------------------------------------------
    ! PC k-space
    call self%pcKs%Initialization( a= self%lattice%GetVp(), n=nk ,meshtype=1,print_=self%getprint() )


  endsubroutine




endmodule
