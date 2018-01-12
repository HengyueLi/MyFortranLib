


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : CPT_ARPES
! OBJECT: TYPE(ARPES)
! USED  : CodeObject , GCE_CPT_para , CE_Green , CreateKspace , LatticeConfig , LaLatticeH , class_numerical_method
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
!                   [sub]
!
!
!                   [sub] GetReducedGkOmegaMatrix(spini,spinj,GpriMatrix,OrbitalList,Karray,g)
!                          integer,intent(in)::spini,spinj
!                          complex*16,intent(in)::GpriMatrix(:,:,:)!(Nomega , self%ns,self%ns)
!                          integer,intent(in)::OrbitalList(:)
!                          real*8,intent(in)::Karray(:,:)           !klist(3,l)   saving all k points
!                          complex*16,intent(out)::g(:,:)        ! g(Omega,size(Karray,2)   )
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
!
!                   [sub] GetReducedGkbyOmegaArrayAndKArray(spini,spinj,OmegaArray,KArray,OrbitalList,g)
!                         see the detals of GetReducedGkOmegaMatrix, input OmegaArray directly.
!
!
!
!
!                   [sub] GetReducedGkOmegaMatrixbyGpara(spini,spinj,No,o1,o2,eta,OrbitalList,Karray,g)
!                         see the details of 'GetReducedGkOmegaMatrix'
!                         the different is to Caluculate Gpri automatically.
!
!                   [sub] GetARPES_spectral_to_file(spini,spinj,No,o1,o2,eta,OrbitalList,Klist,Kstep,file)
!                          integer,intent(in)::spini,spinj,No
!                          real*8,intent(in)::o1,o2,eta
!                          integer,intent(in)::OrbitalList(:)  ! if one do not care orbital, input all orbitals.
!                          real*8,intent(in)::Klist(:,:)           !klist(3,l)
!                          integer,intent(in)::Kstep
!                          character(*),intent(in)::file
!
!                   [sub] GetIntergrationOfGkInKpsace(spini,spinj,OmegaArray,nk3,OrbitalList,g)
!                          implicit none
!                          integer,intent(in)::spini,spinj
!                          complex*16,intent(in)::OmegaArray(:)
!                          integer,intent(in)::nk3(3)
!                          integer,intent(in)::OrbitalList(:)
!                          complex*16,intent(out)::g(Size(OmegaArray))
!                          !--------------------------------------------------------------------
!                          used to calculate total Spectral
!
!                   [sub]  Get_Lattic_Dos_to_file(self,spini,spinj,No,o1,o2,eta,OrbitalList,nk3,file)
!                           implicit none
!                            class(ARPES),intent(inout) :: Self
!                            integer,intent(in)::spini,spinj,No
!                            real*8,intent(in)::o1,o2,eta
!                            integer,intent(in)::OrbitalList(:)
!                            integer,intent(in)::nk3(3)
!                            character(*),intent(in)::file
!                          !--------------------------------------------------------------------
!                           calculate lattice dos and output to file.
!
!                   [sub] Get_K_DOS_to_file(self,spini,spinj,No,o1,o2,eta,OrbitalList,k,file)
!                         calculate the DOS at one K point.
!
!
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
  use class_numerical_method
  implicit none


  type,extends(object)::ARPES
    private

    class(GCECPTpara),pointer   :: CPT      => null()
    class(lh),pointer           :: CPTh     => null()
    class(LaCon),pointer        :: lattice  => null()
    type(CEG)                   :: GreenFun

    !-----------------------
    integer::ns
  contains
    procedure,pass :: Initialization

    procedure,pass::GetReducedGkOmegaMatrix
    procedure,pass::GetReducedGkOmegaMatrixbyGpara
    procedure,pass::GetARPES_spectral_to_file
    procedure,pass::Get_Lattic_Dos_to_file
    procedure,pass::Get_K_DOS_to_file

  endtype


  private::Initialization
  private::GetGkMatrixbyGinvs,GetContractedGk,Get_ARPES_K_Points
  private::GetReducedGkOmegaMatrix
  private::GetReducedGkOmegaMatrixbyGpara
  private::GetARPES_spectral_to_file
  private::Get_Lattic_Dos_to_file
  private::Get_K_DOS_to_file

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


  Subroutine GetGkMatrixbyGinvs(self,GInvsMatrix,spini,spinj,k,Gk)
    implicit none
    class(ARPES),intent(inout) :: Self
    complex*16,intent(in)      :: GInvsMatrix(self%ns,self%ns)
    integer,intent(in)         :: spini,spinj
    real*8,intent(in)          :: k(3)
    complex*16,intent(out)     :: Gk(self%ns,self%ns)
    !-------------------------------------------------------------
    TYPE(nummethod)::f                                            ! ;integer::jc1,jc2;complex*16::m(4,4)

    Gk = GInvsMatrix - self%CPTh%GetTqMatrix(k,spini,spinj)
    if (self%CPTh%GetMeanFieldState())then
       Gk = Gk + self%CPTh%GetMeanField(spini,spinj)
                                                                ! m = self%CPTh%GetMeanField(spini,spinj)
                                                                ! do jc1=1,4;do jc2=1,4
                                                                !    write(*,*)jc1,jc2,real(m(jc1,jc2))
                                                                !  enddo;enddo
                                                                !  stop





    endif
    call f%MatrixInverse(self%ns,gk)
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
        gk = gk  +  GKm(jc1,jc2) * expkx
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
         x = (jc - 1) * dx! + dx/2
         kl(:,jc) = k1 + (k2 - k1) * x
      enddo
    endsubroutine
  endsubroutine



  subroutine GetReducedGkOmegaMatrix(self,spini,spinj,GpriInvsMatrix,OrbitalList,Karray,g)
    implicit none
    class(ARPES),intent(inout) :: Self
    integer,intent(in)::spini,spinj
    complex*16,intent(in)::GpriInvsMatrix(:,:,:)!(Nomega , self%ns,self%ns)
    integer,intent(in)::OrbitalList(:)
    real*8,intent(in)::Karray(:,:)           !klist(3,l)
    complex*16,intent(out)::g(:,:)        ! g(Omega,nk= (len(klist)-1)*kstep   )
    !--------------------------------------------------------------------------------------
    integer::jc,jc1
    complex*16::tempg( size( GpriInvsMatrix,1  )   ),Gt1(Self%ns,self%ns),Gk1(Self%ns,self%ns)
    real*8::k(3)

    do jc = 1 , size(Karray,2)
      k = Karray(:,jc)
      do jc1 = 1 , size( GpriInvsMatrix,1  )
        Gt1 = GpriInvsMatrix(jc1,:,:)
        call GetGkMatrixbyGinvs(self,GInvsMatrix=Gt1,spini=spini,spinj=spinj,k=k,Gk=gk1)
        g(jc1,jc) = GetContractedGk(self,GKm=gk1,k=k,OrbitalList=OrbitalList)
      enddo
    enddo

  endsubroutine



  subroutine GetReducedGkbyOmegaArrayAndKArray(self,spini,spinj,OmegaArray,KArray,OrbitalList,g)
    implicit none
    class(ARPES),intent(inout) :: Self
    integer,intent(in)::spini,spinj
    complex*16,intent(in)::OmegaARRAY(:)!(Nomega)
    real*8,intent(in)::Karray(:,:)     !Karray(3,l)
    integer,intent(in)::OrbitalList(:)
    complex*16,intent(out)::g(:,:)        ! g(Omega,nk= (len(klist)-1)*kstep   )
    !--------------------------------------------------------------------------------------
    integer::jc
    complex*16::Gt(Self%ns,self%ns)
    complex*16::GpriInvs(size(OmegaArray),self%ns,self%ns)
    TYPE(nummethod)::f

    call self%GreenFun%GetGreenMatrix(spini,spinj,Size(OmegaArray),OmegaArray,GpriInvs)
    do jc = 1 , size(OmegaArray)
      Gt = GpriInvs(jc,:,:)
      call f%MatrixInverse(self%ns,Gt)
      GpriInvs(jc,:,:)  = gt
    enddo

   call GetReducedGkOmegaMatrix(self,spini,spinj,GpriInvs,OrbitalList,Karray,g)
  endsubroutine



  Subroutine GetReducedGkOmegaMatrixbyGpara(self,spini,spinj,No,o1,o2,eta,OrbitalList,Karray,g)
    implicit none
    class(ARPES),intent(inout) :: Self
    integer,intent(in)::spini,spinj,No
    real*8,intent(in)::o1,o2,eta
    integer,intent(in)::OrbitalList(:)
    real*8,intent(in)::Karray(:,:)           !klist(3,l)
    complex*16,intent(out)::g(No,   size(karray,2) )        ! g(Omega,nk= (len(klist)-1)*kstep   )
    !--------------------------------------------------------------------------------------
    integer::jc
    complex*16::Omega(No),GpriInvs(No,self%ns,self%ns),Gt(self%ns,self%ns)
    TYPE(nummethod)::f

    Call self%CheckInitiatedOrStop()
    !-------------------------------------------------------------
    ! Get Omega list
     do jc = 1 , No
      Omega(jc) = o1 + (o2-o1)/No * (jc-1) + ( 0._8 , 1.0_8 ) * eta
    enddo
    !-------------------------------------------------------------
   call GetReducedGkbyOmegaArrayAndKArray(self,spini,spinj,Omega,KArray,OrbitalList,g)
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
    real*8::KArray(3,  Kstep* (size( Klist,2 ) - 1)   )
    integer::jc1
    character(32)::Nkc
    character(128)::form

    !---------------------------------------------------------------------------
    ! get All K points from Klist
    Call Get_ARPES_K_Points(self,Klist,Kstep,KArray)               !;write(*,*)KArray

    call GetReducedGkOmegaMatrixbyGpara(self,spini,spinj,No,o1,o2,eta,OrbitalList,KArray,g)
    write(Nkc,*)Kstep* (size( Klist,2 ) - 1)
    form = "("//trim(adjustl(Nkc))//"ES22.14)"
    open(9541,file=trim(adjustl(file)))
    do jc1 = No , 1 , -1
       write(9541,trim(adjustl(form)),advance='no') -imag(g(jc1,:))/pi  ;write(9541,*)
    enddo
    close(9541)
  endsubroutine

  subroutine GetIntergrationOfGkInKpsace(self,spini,spinj,OmegaArray,nk3,OrbitalList,g)
    implicit none
    class(ARPES),intent(inout) :: Self
    integer,intent(in)::spini,spinj
    complex*16,intent(in)::OmegaArray(:)
    integer,intent(in)::nk3(3)
    integer,intent(in)::OrbitalList(:)
    complex*16,intent(out)::g(Size(OmegaArray))
    !--------------------------------------------------------------------
    type(Kspace) :: Ks
    real*8,allocatable::KArray(:,:),Wk(:)
    integer::jc
    complex*16,allocatable::gmatrix(:,:)

    call Ks%Initialization( a= self%lattice%GetVp(), n=nk3 ,meshtype=1,print_=self%getprint() )
    allocate( KArray( 3,ks%getnk() )  ,  Wk(ks%getnk())  ,  gmatrix(Size(OmegaArray),ks%getnk()) )
    do jc = 1 , ks%getnk()
      Karray(:,jc) = ks%getk(jc)
      wk(jc)       = ks%getw(jc)
    enddo

    call GetReducedGkbyOmegaArrayAndKArray(self,spini,spinj,OmegaArray,KArray,OrbitalList,gmatrix)
    g = (0._8,0._8)
    do jc = 1 , ks%getnk()
       g = g + gmatrix(:,jc) * wk(jc)
    enddo
    deallocate( KArray  ,  Wk,  gmatrix )
  endsubroutine

  subroutine Get_Lattic_Dos_to_file(self,spini,spinj,No,o1,o2,eta,OrbitalList,nk3,file)
    implicit none
    class(ARPES),intent(inout) :: Self
    integer,intent(in)::spini,spinj,No
    real*8,intent(in)::o1,o2,eta
    integer,intent(in)::OrbitalList(:)
    integer,intent(in)::nk3(3)
    character(*),intent(in)::file
    !--------------------------------------------------------------------------------------
    real*8,parameter::pi = 3.141592653589793238462643383279502884_8
    integer::jc
    complex*16::g(NO),omega(no)

    !-------------------------------------------------------------
    ! Get Omega list
     do jc = 1 , No
      Omega(jc) = o1 + (o2-o1)/No * (jc-1) + ( 0._8 , 1.0_8 ) * eta
    enddo
    !-------------------------------------------------------------
    call GetIntergrationOfGkInKpsace(self,spini,spinj,Omega,nk3,OrbitalList,g)

    open(9541,file=trim(adjustl(file)))
    do jc = 1,no
       write(9541,*)real(omega(jc)) , -imag(g(jc))/pi
    enddo
    close(9541)
  endsubroutine


  subroutine Get_K_DOS_to_file(self,spini,spinj,No,o1,o2,eta,OrbitalList,k,file)
    implicit none
    class(ARPES),intent(inout) :: Self
    integer,intent(in)::spini,spinj,No
    real*8,intent(in)::o1,o2,eta
    integer,intent(in)::OrbitalList(:)
    real*8,intent(in)::k(3)
    character(*),intent(in)::file
    !--------------------------------------------------------------------------------------
    real*8,parameter::pi=3.141592653589793238462643383279_8
    TYPE(nummethod)::f
    complex*16::g(NO),omega(no)
    complex*16::GpriInvs(No,self%ns,self%ns),Gt(self%ns,self%ns),gt2(self%ns,self%ns)
    integer::jc
    !-------------------------------------------------------------
    ! Get Omega list
     do jc = 1 , No
      Omega(jc) = o1 + (o2-o1)/No * (jc-1) + ( 0._8 , 1.0_8 ) * eta
    enddo
    !-------------------------------------------------------------
    call self%GreenFun%GetGreenMatrix(spini,spinj,No,Omega,GpriInvs)
    open(200,file=trim(adjustl(file)))
    do jc = 1 , no
      Gt = GpriInvs(jc,:,:)
      call f%MatrixInverse(self%ns,Gt)
      call GetGkMatrixbyGinvs(self,Gt,spini,spinj,k,Gt2)
      g(jc) = GetContractedGk(self,Gt2,k,OrbitalList)
      write(200,*)real(omega(jc)),-imag(g(jc))/pi
    enddo
    close(200)
  endsubroutine




endmodule
