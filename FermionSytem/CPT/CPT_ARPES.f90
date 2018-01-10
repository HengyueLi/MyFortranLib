


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : CPT_ARPES
! OBJECT: TYPE(ARPES)
! USED  : CodeObject , GCE_CPT_para , CE_Green , CreateKspace , LatticeConfig , LaLatticeH
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
  use LaLatticeH
  implicit none


  type,extends(object)::ARPES
    private

    class(GCECPTpara),pointer   :: CPT      => null()
    class(lh),pointer           :: CPTh     => null()
    class(LaCon),pointer        :: lattice  => null()
    type(CEG)                   :: GreenFun
    ! class(GreenPara)            :: Gp
    type(Kspace)                :: pcKs           ! choose prH basis.


    !-----------------------
    integer::ns
  contains
    procedure,pass :: Initialization

  endtype


  private::Initialization
  private::GetAkMatrixFile,GetGkMatrixbyG,GetContractedGk


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
    self%lattice => Self%CPT%GetLatticePointer()

    !---------------------------------------------------------------------------
    ! Green function
    call self%GreenFun%Initialization( solver = self%CPT%GetSolverPointer() ,&
                                                GP = GP,PRINT_=self%getprint() )
    !---------------------------------------------------------------------------
    ! PC k-space
    call self%pcKs%Initialization( a= self%lattice%GetVp(), n=nk ,meshtype=1,print_=self%getprint() )
    !---------------------------------------------------------------------------
    ! cpth
    self%CPTh => self%cpt%GetCPTHpointer()
    !-----------------------
    ! some other integer.
    self%ns = self%lattice%GetNs()

  endsubroutine

  Subroutine GetAkMatrixFile(self,SpinI,SpinJ,o1,o2,No,eta,NKstep,Klist)
    implicit none
    class(ARPES),intent(inout) :: self
    integer,intent(in):: SpinI,SpinJ
    real*8,intent(in) :: o1,o2,eta
    integer,intent(in):: No,NKstep
    real*8,intent(in) :: Klist(:,:) ! Klist(3,N)
    !-------------------------------------------------------------
    integer::jc
    complex*16::Omega(No),G(No,self%ns,self%ns)


    !-------------------------------------------------------------
    ! Get Omega list
     do jc = 1 , No
      Omega(jc) = o1 + (o2-o1)/No * (jc-1) + ( 0._8 , 1.0_8 ) * eta
     enddo
    !-------------------------------------------------------------
    ! Get G matrix
    call self%GreenFun%GetGreenMatrix(spini,spinj,NO,Omega,G)
    !-------------------------------------------------------------

  endsubroutine


  Subroutine GetGkMatrixbyG(self,Gmatrix,spini,spinj,k,Gk)
    implicit none
    class(ARPES),intent(inout) :: Self
    complex*16,intent(in)      :: Gmatrix(self%ns,self%ns)
    integer,intent(in)         :: spini,spinj
    real*8,intent(in)          :: k(3)
    complex*16,intent(out)     :: Gk(self%ns,self%ns)
    !-------------------------------------------------------------
    Gk = Gmatrix - self%CPTh%GetTqMatrix(k,spini,spinj)
  endsubroutine


  !  calculate gk for orbital list OrbitalList. Only count the site belongs to OrbitalList
  ! GKm is GK matrix .
  ! gk = 1/L * sum_{ij} G_{ij}(k,\omega)
  function GetContractedGk(self,GKm,k,OrbitalList) result(gk)
    implicit none
    class(ARPES),intent(inout) :: Self
    complex*16,intent(in)::Gkm(self%ns,self%ns)
    real*8,intent(in)    :: k(3)
    integer,intent(in)::OrbitalList(:)
    complex*16::gk
    !-------------------------------------------------------------
    integer::jc1,jc2,Nllist
    complex*16::expkx
    !------------------------
    gk = (0._8,0._8)
    !------------------------
    Nllist = size(OrbitalList)
    do jc1 = 1 , self%ns   ;  do jc2 = 1 , self%ns
      if (belongs(jc1,Nllist,OrbitalList) .and. belongs(jc2,Nllist,OrbitalList))then
        expkx=zexp(sum((self%lattice%GetSiteRealP(jc2)-self%lattice%GetSiteRealP(jc1))* k)*(0._8,1._8))
        gk = gk  +  expkx * GKm(jc1,jc2)
      endif
    enddo                  ;  enddo
    gk = gk / self%Ns
  contains

    logical function belongs(i,N,list)
      implicit none
      integer::i,n,list(n)
      !---------------------
      integer::jc
      belongs = .false.
      do jc = 1 , n
        if (i==list(jc))then
          belongs = .true.
          goto 99
        endif
      enddo
  99  continue
    endfunction
  endfunction











endmodule
