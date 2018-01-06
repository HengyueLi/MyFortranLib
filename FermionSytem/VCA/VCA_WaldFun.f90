


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : VCA_WaldFun
! OBJECT: TYPE(waldf)
! USED  : CodeObject , LaLatticeH , VCA_DeltaH , CEsolver , CE_Green , CreateKspace , fermion_table , class_numerical_method
! DATE  : 2018-01-05
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!
!             Evaluate the Wald functional of the system. That is :
!                               \Omega = \Omega' + Trln[-G] - Trln[-G'] = \Omega' - I
!             where  I = T * \sum_{omega_n}\sum_k lndet[1-Vkd * G']
!
!             basically, we seperate the Integration into two part: I = I_i + I_s
!             where
!                  I_i is an integration part and I_s is an summation part.
!             The different of I_i and I_s can be concludes as the only difference of integration weight.
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
!                   [sub] Initialization(Ta,SP,GP,CPTH,dH,jobi,jobr,print_,show_)
!                          class(table),target,intent(inout):: Ta
!                          class(SolverPara),intent(in)     :: SP
!                          class(GreenPara),intent(in)      :: GP
!                          class(LH),target,intent(inout)   :: CPTh
!                          class(VCAdH),target,intent(inout):: dH
!                          integer,intent(in)               :: jobi(:)
!                          real*8 ,intent(in)               :: jobr(:)
!                          integer,intent(in),optional      :: print_
!                          integer,intent(in),optional      :: show_
!
!                          -----------------------
!                          jobi(1)   = method
!                          jobi(2:4) = n1,n2,n3 the samoling points on each direction in Kspace.
!
! avalable gets:
!
!                   [fun] GetLatticeOmegaPerSite()
!
!
! avalable is :
!                  ![fun] i
! others      :
!                   [sub] A
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module CreatArcPath
  use class_integration
  implicit none


  type::VCAWandArc
     private
     logical :: Initiated = .false.
     real*8  :: R
     real*8  :: theta   ! a small arc
     integer :: nr
     integer :: ntheta  ! correcponding to the small arc

    !---------
    integer::print = 6

  contains
    procedure,pass::Initialization
    final::finalization
    procedure,pass::GetOmegaWeigth
  endtype

  private::Initialization,UnInitialization,Finalization

contains


  subroutine Initialization(self,R,theta,nr,ntheta)
    implicit none
    class(VCAWandArc),intent(inout)::self
    real*8,intent(in)::R,theta
    integer,intent(in)::nr,ntheta
    !-----------------------------------------------
    call UnInitialization(self) ; self%Initiated = .true.

    self%r = r
    self%theta = theta
    self%nr = nr
    self%ntheta = ntheta


  endsubroutine


  subroutine UnInitialization(self)
    implicit none
    class(VCAWandArc),intent(inout)::self
    !-----------------------------------------------
    if (self%Initiated)then  ; self%Initiated = .false.

    endif
  endsubroutine

  impure elemental subroutine finalization(self)
    implicit none
    type(VCAWandArc),intent(inout)::self
    !-----------------------------------------------
    call UnInitialization(self)
  endsubroutine



  subroutine GetOmegaWeigth(self,Omega,Wei)
    implicit none
    class(VCAWandArc),intent(inout)::self
    complex*16,intent(out)::Omega(self%ntheta*2+self%nr*2),Wei(self%ntheta*2+self%nr*2)
    real*8,parameter::pi=3.141592653589793238462643383279_8
    !-----------------------------------------------
    TYPE(Integration)::I1,I2
    real*8::x1(self%nr),w1(self%nr),x2(self%ntheta),w2(self%ntheta)
    real*8::x(self%nr+self%ntheta),w(self%nr+self%ntheta)
    real*8::t1,t2,t3
    integer::jc1,N



    call I1%Initialization(self%nr     ,1)
    call I2%Initialization(self%ntheta ,1)

    N = self%nr + self%ntheta
    t1 = 0._8    ;  t2 = (90._8 - self%theta)/180._8*pi  ;  t3 = 0.5_8*pi

    call I1%get_xw(xs=t1,xe= t2 ,x=x1,w=w1)
    call I2%get_xw(xs=t2,xe= t3 ,x=x2,w=w2)

    do jc1 = 1 , self%nr
      x(jc1) = x1(jc1)          ; w(jc1) = w1(jc1)
    enddo
    do jc1 = 1 , self%ntheta
      x(jc1+self%nr ) = x2(jc1) ; w(jc1+self%nr ) = w2(jc1)
    enddo
    w = w * self%r

    !write(*,*)x;stop

    do jc1 = 1 , N
       call GetXW_left(self%r,x(jc1),omega(jc1),wei(jc1))
       wei(jc1) = wei(jc1) * w(jc1)
    enddo
    !write(*,*)omega(1:n);stop


    do jc1 = 1 , N
      omega(2*N-jc1+1) = -real(omega(jc1)) + imag(omega(jc1)) * (0._8,1._8)
      wei(  2*N-jc1+1) =  real(wei(jc1)  ) - imag(wei(jc1)  ) * (0._8,1._8)
    enddo



  contains
    subroutine GetXW_left(r,phi,x,w)
      implicit none
      real*8,intent(in)::r,phi
      complex*16,intent(out)::x,w
      !---------------------------
      real*8::a,b

      a = 3*pi/4 - phi
      b = a - pi/2
      x = dcos(a) * r + dsin(a) * r * (0._8,1._8)  + (-r/dsqrt(2.0_8)*(1._8,1._8))
      w = dcos(b)     + dsin(b) * (0._8,1._8)

    endsubroutine
  endsubroutine

endmodule









module VCA_WaldFun
  use CodeObject
  use LaLatticeH
  USE VCA_DeltaH
  use CEsolver
  use CE_Green
  use CreatArcPath
  use phyfunc
  use class_integration
  use CreateKspace
  use fermion_table
  use class_numerical_method
  implicit none


  integer,parameter::NparaMax = 50

  type,extends(object)::waldf
    private
    integer::jobi(NparaMax)
    real*8 ::jobr(NparaMax)
    class(table),pointer:: Ta => null()
    type(SolverPara) :: SP
    TYPE(GreenPara)  :: Gp

    class(LH),pointer    :: cptH => NULL()
    class(VCAdH),pointer :: dH   => null()
    !------------------------------------------
    real*8::temperature  ! we always need the temperature. It can be set to be very small.
    integer:: IntNomega
    complex*16,allocatable:: IntOmega(:)
    complex*16,allocatable:: IntOmegaWeight(:)

    integer::ns
    complex*16,allocatable :: SavedG(:,:,:)     ! G( Nomega , 2NS , 2NS   )
    real*8                 :: OmegaPri
    TYPE(Kspace)::k

  contains
    procedure,pass::Initialization
    final::Finalization

    procedure,pass::GetLatticeOmegaPerSite
  endtype


  private::Initialization,Finalization
  private::AllocateIntergrationPath,AllocateIntergrationPath_method2
  private::GetSavingGandGrandPotetial
  private::FuncF,GetVkMatrix,GetI,GetLatticeOmegaPerSite

contains

  subroutine Initialization(self,Ta,SP,GP,CPTH,dH,jobi,jobr,print_,show_)
    implicit none
    class(waldf),intent(inout)       :: self
    class(table),target,intent(inout):: Ta
    class(SolverPara),intent(in)     :: SP
    class(GreenPara),intent(in)      :: GP
    class(LH),target,intent(inout)   :: CPTh
    class(VCAdH),target,intent(inout):: dH
    integer,intent(in)               :: jobi(:)
    real*8 ,intent(in)               :: jobr(:)
    integer,intent(in),optional      :: print_
    integer,intent(in),optional      :: show_
    !-------------------------------------------------------------
    real*8 ::ba(3,3)
    integer::kn(3)

    call Finalization(self)
    call self%SetInitiated(.true.)
    if (present(print_)) call self%SetPrint(print_)
    if (present(show_))  call self%SetShow(show_)

    self%Ta   => Ta
    self%SP   =  SP
    self%Gp   =  gp
    self%cptH => cptH
    self%dh   => dH

    if (  max( size(jobi) , size(jobr) ).gt.NparaMax )then
       write(self%getprint(),*)"ERROR: NparaMax is too small in Waldf";stop
    endif
    self%jobi(1:size(jobi)) = jobi
    self%jobr(1:size(jobr)) = jobr
    self%temperature = SP%temperature
    self%ns          = ta%get_ns()


    ba = self%cptH%GetVbasis()
    kn = self%jobi(2:4)
    call self%k%Initialization( a = ba,n=kn,meshtype=1,&
                                    print_=self%getprint(),show_=self%getshow() )


    call AllocateIntergrationPath(self)

  endsubroutine


  impure elemental subroutine Finalization(self)
    implicit none
    type(waldf),intent(inout) :: self
    !---------------------------------------
    if (  self%IsInitiated()  )then
       deallocate(  self%IntOmega , self%IntOmegaWeight , self%SavedG )
    endif
  endsubroutine

  subroutine AllocateIntergrationPath(self)
    implicit none
    class(waldf),intent(inout)  :: self
    !---------------------------------------
    select case(self%jobi(1))
    case(2)
      call AllocateIntergrationPath_method2(self)
    case default
      write(self%getprint(),*)"ERROR: Unknow method in waldf@AllocateIntergrationPath"
      write(self%getprint(),*)"input method = jobi(1)=",self%jobi(1)
    endselect

  endsubroutine


  ! Arc intergration
  ! jobr(1) = arc R. The scale of circle
  ! jobr(2) = theta. A small angle used to treat the near-zero problem.
  ! jobi(5) = sampling points in R part(large part).
  ! jobi(6) = sampling points in theta regime.
  subroutine AllocateIntergrationPath_method2(self)
    implicit none
    class(waldf),intent(inout)  :: self
    !----------------------------------------------
    real*8,parameter::Pi=3.141592653589793238462643383279_8
    !==============================================
    TYPE(phyfun)::f
    integer::jc
    INTEGER::I,jc1,jc2,con   ;REAL(8)::T,XI,WI,I0,I1,I2
    real*8,allocatable::x(:),w(:)
    TYPE(VCAWandArc)::VcAW

    call VcAW%Initialization(R=self%jobr(1),theta=self%jobr(2)&
      ,nr=self%jobi(5),ntheta=self%jobi(6) )


    self%IntNomega = sum(self%jobi(5:6))  * 2
    allocate(  self%IntOmega(self%IntNomega)                                     )
    allocate(  self%IntOmegaWeight(self%IntNomega)                               )
    allocate(  self%SavedG(self%IntNomega,self%ta%get_ns()*2,self%ta%get_ns()*2) )

    call VcAW%GetOmegaWeigth(self%IntOmega,self%IntOmegaWeight)

     !------------------------
     ! absorb f(omega) into weight
     do jc2 = 1 , self%IntNomega
        self%IntOmegaWeight(jc2) = (0._8,1._8)/2._8/pi*self%IntOmegaWeight(jc2) &
                        * f%FermionFunc(self%IntOmega(jc2),self%temperature)
     enddo
  endsubroutine



  subroutine GetSavingGandGrandPotetial(self)
    implicit none
    class(waldf),intent(inout)  :: self
    !----------------------------------------------
    type(CES)::solver
    TYPE(CEG)::G
    TYPE(Ham)::H
    integer::ns,jc1,jc2                                     !,jc,para(8),wtp
    complex*16::tempG(self%IntNomega,self%ns,self%ns)       !;class(FermOper),pointer::p

    ns = self%ta%get_ns()


    !----------------------------------------------
    !  get H
    call H%Initialization(ns,self%getprint())
    call H%StartAppendingInteraction()
    Call self%cptH%AppendLocalDataToHam(H)
    Call self%dH%AppendDataToHam(H)
    call h%EndAppendingInteraction()



    ! open(99,file="/home/hengyue/Desktop/eivca.dat")
    ! wtp = 99
    ! do jc = 1 , H%GetOptN()
    !   p => h%GetOptact(jc)
    !   write(wtp,*)jc,"/",H%GetOptN()
    !   write(wtp,*)"optid=",p%get_optid()
    !   write(wtp,*)"V=",H%GetOptV(jc)
    !   para = p%get_para()
    !   write(wtp,*)"para=",para(1:6)
    !
    ! enddo   ;stop


    !-----------------------------------------------

    Call solver%Initialization( self%SP , self%Ta, H , self%getprint() , self%getshow() )
    Call solver%SynchronizeWithHamiltonian()

    self%OmegaPri = solver%GetGrandPotential()


    Call G%Initialization(solver,self%Gp,self%getprint() , self%getshow()+1)

    call G%GetGreenMatrix(spini=0,spinj=0,NOmega=self%IntNomega,Omega=self%IntOmega&
                    , GM=tempG)     ;  self%SavedG(:,1:ns,1:ns) = tempG
    call G%GetGreenMatrix(spini=1,spinj=1,NOmega=self%IntNomega,Omega=self%IntOmega&
                    , GM=tempG)     ;  self%SavedG(:,ns+1:2*ns,ns+1:2*ns)= tempG
    call G%GetGreenMatrix(spini=0,spinj=1,NOmega=self%IntNomega,Omega=self%IntOmega&
                    , GM=tempG)     ;  self%SavedG(:,1:ns,ns+1:2*ns)  = tempG

    forall(jc1=1:self%IntNomega)
       self%SavedG(jc1,ns+1:2*ns,1:ns) = transpose( conjg( self%SavedG(jc1,1:ns,ns+1:2*ns) ) )
    endforall

  endsubroutine

  ! V(2ns,2ns)
  function GetVkMatrix(self,k,alpha) result(vk)
    implicit none
    class(waldf),intent(inout)  :: self
    real*8,intent(in)           :: k(3)
    real*8,intent(in)           :: alpha
    complex*16::vk(self%ns*2,self%ns*2)
    !----------------------------------------------
    vk = self%cptH%GetSpinSuppresedTq(k) * alpha - self%dH%GetDeltaMatrixSpinSupp()
  endfunction



  ! F = \sum_k \ln\det[1-V_{k\Delta}G']
  complex*16 function FuncF(self,omegaid)
    IMPLICIT NONE
    class(waldf),intent(inout)  :: self
    integer,intent(in)::omegaid
    !------------------------------
    TYPE(nummethod)::f
    integer::jck ,jc , symm ,spin ,i1,i2
    complex*16::temp(self%ns*2,self%ns*2),Vk(self%ns*2,self%ns*2)
    complex*16::G(self%ns*2,self%ns*2) ,temp1 ,TempS2(self%ns,self%ns)
    real*8::k(3),w
    !--------- initiate ------------
    FuncF = (0._8,0._8)
    !----------get G----------------
    G = self%SavedG(omegaid,:,:)
    !----------sum k---------------
    symm = self%Ta%get_symmetry()
    do jck = 1 , self%k%getnk()
       k =  self%k%getk(jck)
       w =  self%k%getw(jck)
       temp = (0._8,0._8)
       do jc = 1 , self%ns * 2
         temp(jc,jc) = (1._8,0._8)
       enddo
       vk = GetVkMatrix(self,k, self%dh%GetAlpha()  )
       temp = temp - matmul( Vk , G )

       !------------------------------------------------------------------------
       ! if symmtry = 2 ,  there is no cross terms. The off diagonal terms =0
       if (symm==2)then
         temp1 = (0._8,0._8)
         do spin = 0 , 1
           i1 = spin * self%ns + 1  ; i2 = i1 + self%ns - 1
           TempS2 = temp(i1:i2,i1:i2)
           temp1 = temp1 + zlog( f%get_determinant_of_Matrix(self%ns,TempS2) )
         enddo
       else
         temp1 = zlog( f%get_determinant_of_Matrix(self%ns*2,temp) )
       endif
       !------------------------------------------------------------------------

       FuncF = FuncF +  temp1 * w
     enddo
  endfunction


  real*8 function GetI(self)
    IMPLICIT NONE
    class(waldf),intent(inout)  :: self
    !--------------------------------------
    integer::jc
    complex*16::ret

    ret = (0._8,0._8)
    do jc = 1 , self%IntNOmega
       ret = ret + FuncF(self,jc) * self%IntOmegaWeight(jc)
    enddo
    GetI = 2._8 * real(ret)
  endfunction

  real*8 function GetLatticeOmegaPerSite(self)
    IMPLICIT NONE
    class(waldf),intent(inout)  :: self
    !--------------------------------------
    call self%CheckInitiatedOrStop()
    call GetSavingGandGrandPotetial(self)
    GetLatticeOmegaPerSite = self%OmegaPri - GetI(self)             !;write(*,*)self%OmegaPri, GetI(self),666
    GetLatticeOmegaPerSite = GetLatticeOmegaPerSite / self%ns
  endfunction






endmodule
