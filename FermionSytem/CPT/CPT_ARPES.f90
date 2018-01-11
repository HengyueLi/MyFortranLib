


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
!                   [sub] Initialization(CPTpara,GP)
!                          class(GCECPTpara),target    :: CPTpara
!                          class(GreenPara),intent(in) :: GP
!
!
!
! avalable gets:
!
!                   [sub] GetReducedGkOmegaMatrix(spini,spinj,GpriMatrix,OrbitalList,Klist,Kstep,g)
!                          integer,intent(in)::spini,spinj
!                          complex*16,intent(in)::GpriMatrix(:,:,:)!(Nomega , self%ns,self%ns)
!                          integer,intent(in)::OrbitalList(:)
!                          real*8,intent(in)::Klist(:,:)           !klist(3,l)
!                          integer,intent(in)::Kstep
!                          complex*16,intent(out)::g(:,:)        ! g(Omega,nk= (len(klist)-1)*kstep   )
!                          !----------------------------------------------------------------------------
!                          input:
!                               GpriMatrix(No,ns,ns) to fix y data, and input Klist to fix x data
!                               OrbitalList to check orbital index
!                          ouput: g (matrix)
!                                          ┌────┐
!                                   omega  │    │
!                                          └────┘         (Notice the order of Omega)
!                                             k
!
!                   [sub] GetReducedGkOmegaMatrixbyGpara(spini,spinj,No,o1,o2,eta,OrbitalList,Klist,Kstep,g)
!                         see the details of 'GetReducedGkOmegaMatrix'
!                         the different is to Caluculate Gpri automatically.
!
!                   [sub] GetARPES_spectral_to_file(spini,spinj,No,o1,o2,eta,OrbitalList,Klist,Kstep,file)
!                          integer,intent(in)::spini,spinj,No
!                          real*8,intent(in)::o1,o2,eta
!                          integer,intent(in)::OrbitalList(:)
!                          real*8,intent(in)::Klist(:,:)           !klist(3,l)
!                          integer,intent(in)::Kstep
!                          character(*),intent(in)::file
!
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

    procedure,pass::GetReducedGkOmegaMatrix
    procedure,pass::GetReducedGkOmegaMatrixbyGpara
    procedure,pass::GetARPES_spectral_to_file

  endtype


  private::Initialization
  private::GetAkMatrixFile,GetGkMatrixbyG,GetContractedGk,Get_ARPES_K_Points
  private::GetReducedGkOmegaMatrixbyGpara,GetReducedGkOmegaMatrix
  private::GetARPES_spectral_to_file


contains


  subroutine Initialization(self  ,CPTpara,GP)
    implicit none
    class(ARPES),intent(out)    :: self
    class(GCECPTpara),target    :: CPTpara
    class(GreenPara),intent(in) :: GP
    ! integer,intent(in)          :: nk(3)
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
    ! call self%pcKs%Initialization( a= self%lattice%GetVp(), n=nk ,meshtype=1,print_=self%getprint() )
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
    logical::ISbelo1,iSbelo2
    !------------------------
    gk = (0._8,0._8)
    !------------------------
    Nllist = size(OrbitalList)
    do jc1 = 1 , self%ns   ;  do jc2 = 1 , self%ns
      ISbelo1 = belongs(self%lattice%GetOrbitIndex(jc1),Nllist,OrbitalList)
      ISbelo2 = belongs(self%lattice%GetOrbitIndex(jc2),Nllist,OrbitalList)
      if ( ISbelo1 .and. ISbelo2 )then
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


  subroutine Get_ARPES_K_Points(self,Klist,NKstep,KPointsOut)
    implicit none
    class(ARPES),intent(inout) :: Self
    real*8,intent(in)::Klist(:,:) !(3,:)
    integer,intent(in)::NKstep
    real*8,intent(out)::KPointsOut(:,:)
    !---------------------------------------------------------------------------
    integer::NK,Klistlen,jc1
    real*8::tempK(3,NKstep)
    !--------------
    ! check interface
    if (size(Klist,1).ne.3)then
      write(self%getprint(),*)"size of first script in Klist=",size(Klist,1),", which can only be 3"
      stop
    endif
    Klistlen = size(Klist     ,2)
    NK       = size(KPointsOut,2)
    IF (NK.ne.(Klistlen-1)*NKstep)then
      write(self%getprint(),*)"Should have: NKstep * (len(klist) - 1) = len(KPointsOut) in "
      write(self%getprint(),*)"Get_ARPES_K_Points@CPT_ARPES"  ;stop
    endif
    !--------------
    do jc1 = 1 , Klistlen-1
      Call TwoKline(Klist(:,jc1),Klist(:,jc1+1),NKstep,tempK)
      KPointsOut(:,  (jc1-1)*NKstep + 1 : jc1*NKstep     ) = tempK
    enddo
  contains
    subroutine TwoKline(k1,k2,nk,kl)
      implicit none
      integer,intent(in)::nk
      real*8::k1(3),k2(3),kl(3,nk)
      !-------------------------------
      integer:: jc
      real*8::x,dx
      dx = 1._8/nk
      do jc = 1 , nk
         x = (jc - 1) * dx + dx/2
         kl(:,jc) = k1 + (k2 - k1) * x
      enddo
    endsubroutine
  endsubroutine

  subroutine GetReducedGkOmegaMatrix(self,spini,spinj,GpriMatrix,OrbitalList,Klist,Kstep,g)
    implicit none
    class(ARPES),intent(inout) :: Self
    integer,intent(in)::spini,spinj
    complex*16,intent(in)::GpriMatrix(:,:,:)!(Nomega , self%ns,self%ns)
    integer,intent(in)::OrbitalList(:)
    real*8,intent(in)::Klist(:,:)           !klist(3,l)
    integer,intent(in)::Kstep
    complex*16,intent(out)::g(:,:)        ! g(Omega,nk= (len(klist)-1)*kstep   )
    !--------------------------------------------------------------------------------------
    integer::Nk,jc,jc1,No
    real*8,allocatable::AllKpoints(:,:)
    complex*16::tempg( size( GpriMatrix,1  )   ),Gt1(Self%ns,self%ns),Gk1(Self%ns,self%ns)
    real*8::k(3)

    NK = ( size(klist,2) - 1 ) * Kstep
    No = size( GpriMatrix,1  )

    Allocate( AllKpoints(3,nk) )

    Call Get_ARPES_K_Points(self,Klist,Kstep,AllKpoints)

    do jc = 1 , nk
      k = AllKpoints(:,jc)
      do jc1 = 1 , No
        Gt1 = GpriMatrix(jc1,:,:)
        call GetGkMatrixbyG(self,Gmatrix=Gt1,spini=spini,spinj=spinj,k=k,Gk=gk1)
        g(jc1,jc) = GetContractedGk(self,GKm=gk1,k=k,OrbitalList=OrbitalList)
      enddo
    enddo

    deallocate(AllKpoints)
  endsubroutine


  Subroutine GetReducedGkOmegaMatrixbyGpara(self,spini,spinj,No,o1,o2,eta,OrbitalList,Klist,Kstep,g)
    implicit none
    class(ARPES),intent(inout) :: Self
    integer,intent(in)::spini,spinj,No
    real*8,intent(in)::o1,o2,eta
    integer,intent(in)::OrbitalList(:)
    real*8,intent(in)::Klist(:,:)           !klist(3,l)
    integer,intent(in)::Kstep
    complex*16,intent(out)::g(No,   Kstep* (size( Klist,2 ) - 1) )        ! g(Omega,nk= (len(klist)-1)*kstep   )
    !--------------------------------------------------------------------------------------
    integer::jc
    complex*16::Omega(No),Gpri(No,self%ns,self%ns)

    Call self%CheckInitiatedOrStop()
    !-------------------------------------------------------------
    ! Get Omega list
     do jc = 1 , No
      Omega(jc) = o1 + (o2-o1)/No * (jc-1) + ( 0._8 , 1.0_8 ) * eta
    enddo
    !-------------------------------------------------------------
    ! Get G matrix
    call self%GreenFun%GetGreenMatrix(spini,spinj,NO,Omega,Gpri)
    !-------------------------------------------------------------
    Call GetReducedGkOmegaMatrix(self,spini,spinj,Gpri,OrbitalList,Klist,Kstep,g)
  endsubroutine


  subroutine GetARPES_spectral_to_file(self,spini,spinj,No,o1,o2,eta,OrbitalList,Klist,Kstep,file)
    implicit none
    class(ARPES),intent(inout) :: Self
    integer,intent(in)::spini,spinj,No
    real*8,intent(in)::o1,o2,eta
    integer,intent(in)::OrbitalList(:)
    real*8,intent(in)::Klist(:,:)           !klist(3,l)
    integer,intent(in)::Kstep
    character(*),intent(in)::file
    !--------------------------------------------------------------------------------------
    real*8,parameter::pi = 3.141592653589793238462643383279502884_8
    complex*16::g(No,   Kstep* (size( Klist,2 ) - 1) )
    integer::jc1,Nk
    character(32)::Nkc
    character(128)::form

    call GetReducedGkOmegaMatrixbyGpara(self,spini,spinj,No,o1,o2,eta,OrbitalList,Klist,Kstep,g)
    Nk = Kstep* (size( Klist,2 ) - 1)
    write(Nkc,*)nk
    form = "("//trim(adjustl(Nkc))//"ES22.14)"
    open(9541,file=trim(adjustl(file)))
    do jc1 = No , 1 , -1
       write(9541,trim(adjustl(form)),advance='no') -imag(g(jc1,:))/pi  ;write(9541,*)
    enddo
    close(9541)
  endsubroutine












endmodule
