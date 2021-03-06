


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : VCA_WaldFun
! OBJECT: TYPE(waldf)
! USED  : CodeObject , LaLatticeH , VCA_DeltaH , CEsolver , CE_Green , CreateKspace , fermion_table , class_numerical_method
! DATE  : 2018-01-27
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
!                          select case :
!                          ===========================
!                          method = 2
!                          Arc intergration
!                            jobr(1) = arc R. The scale of circle
!                            jobr(2) = theta. A small angle used to treat the near-zero problem.
!                            jobi(5) = sampling points in R part(large part).
!                            jobi(6) = sampling points in theta regime.
!                          ===========================
!                          method = 3
!                          square intergration   (select summation number on image axis)
!                            jobi(5):Nx  (from -R to R)
!                            jobi(6):Ny  in each side
!                            jobr(1):R   (range to enclose regime)
!                            jobi(7) is summation number on image axis
!                          ===========================
!                          method = 4
!                          square intergration   (select summation length on image axis)
!                            jobi(5):Nx  (from -R to R)
!                            jobi(6):Ny  in each side
!                            jobr(1):R   (range to enclose regime)
!                            jobr(2) is summation length on image axis
!
!
!
!
! avalable gets:
!
!                   [fun] GetLatticeOmegaPerSite()
!                         Omega
!
!                   [fun] GetClusterOmegaPriPerSite()
!                         Omega'
!
!                   [fun] GetLatticeIPerSite()
!                         Omega = Omega' - I
!                         return I
!
!                   [fun] GetOmegaValue(typeid)
!                         return real*8
!                         integer::typeid
!                         This is a inteface for all the get above.
!                           typeid = 1 : GetLatticeOmegaPerSite()
!                           typeid = 2 : GetClusterOmegaPriPerSite()
!                           typeid = 3 : GetLatticeIPerSite()
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
    !------------------------------------------
    integer::ns
    complex*16,allocatable :: SavedG(:,:,:)     ! G( Nomega , 2NS , 2NS   )
    real*8                 :: OmegaPri
    TYPE(Kspace)::k
    !-----------------------------------------

  contains
    procedure,pass::Initialization
    final::Finalization

    procedure,pass::GetLatticeOmegaPerSite
    procedure,pass::GetClusterOmegaPriPerSite
    procedure,pass::GetLatticeIPerSite
    procedure,pass::GetOmegaValue
  endtype


  private::Initialization,Finalization
  private::AllocateIntergrationPath
  private::AllocateIntergrationPath_method2
  private::AllocateIntergrationPath_method34
  private::GetSavingGandGrandPotetial
  private::FuncF,GetVkMatrix,GetI,GetLatticeOmegaPerSite
  private::ReccorectI
  private::GetClusterOmegaPriPerSite
  private::GetLatticeIPerSite
  private::GetOmegaValue

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
    case(3,4)
      call AllocateIntergrationPath_method34(self)
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
    allocate(  self%IntOmega(-1:self%IntNomega)                                     )
    allocate(  self%IntOmegaWeight(-1:self%IntNomega)                               )
    allocate(  self%SavedG(-1:self%IntNomega,self%ta%get_ns()*2,self%ta%get_ns()*2) )

    call VcAW%GetOmegaWeigth(self%IntOmega(1:self%IntNomega),&
                                          self%IntOmegaWeight(1:self%IntNomega) )

     !------------------------
     ! absorb f(omega) into weight
     do jc2 = 1 , self%IntNomega
        self%IntOmegaWeight(jc2) = (0._8,1._8)/2._8/pi*self%IntOmegaWeight(jc2) &
                        * f%FermionFunc(self%IntOmega(jc2),self%temperature)
     enddo

     !----------------
     self%IntOmega(-1) =  0.0002 *  (1.0_8,1._8)/dsqrt(2.0_8)
     self%IntOmega( 0) =  0.0001 *  (1.0_8,1._8)/dsqrt(2.0_8)
  endsubroutine



  !jobi(5):Nx  (from -R to R)
  !jobi(6):Ny  in each side
  !jobr(1):R   (range to enclose regime)
  !
  ! for method3 : jobi(7) is summation number on image axis
  ! for method4 : jobr(2) is summation length on image axis
  subroutine AllocateIntergrationPath_method34(self)
    implicit none
    class(waldf),intent(inout)  :: self
    !----------------------------------------------
    real*8,parameter::Pi=3.141592653589793238462643383279_8
    complex*16,parameter::CI = (0._8,1._8)
    !==============================================
    TYPE(Integration)::intex,intey
    TYPE(phyfun)::f
    real*8::xx(self%jobi(5)),wx(self%jobi(5)),xy(self%jobi(6)),wy(self%jobi(6)),LenY
    integer::SumN,Seg1,Seg2,Seg3,Seg4,jc
    real*8::T

    T = self%temperature

    call intex%Initialization(  self%jobi(5) , 1)
    call intey%Initialization(  self%jobi(6) , 1)


    if ( self%jobi(1) == 3 )then
      SumN = self%jobi(7)
    else
      SumN = int( ( self%jobr(2)/pi/T + 1.0_8 ) / 2 )
    endif

    !---------------------------
    ! at least 1 is required
    if ( SumN < 1 ) SumN = 1
    !---------------------------
    ! integration range on Y
    LenY = 2 * SumN * Pi * T


    self%IntNomega = SumN + 2 * self%jobi(6) + self%jobi(5)

    Seg1 = self%jobi(5)
    Seg2 = Seg1 + self%jobi(6)
    Seg3 = Seg2 + self%jobi(6)
    Seg4 = self%IntNomega

    call intex%get_xw( -self%jobr(1) , self%jobr(1) ,xx ,wx)
    call intey%get_xw( 0.0_8         , LenY         ,xy ,wy)

    allocate(  self%IntOmega(self%IntNomega)                                     )
    allocate(  self%IntOmegaWeight(self%IntNomega)                               )
    allocate(  self%SavedG(self%IntNomega,self%ta%get_ns()*2,self%ta%get_ns()*2) )
    !       1      ->   Seg1
    !     .............................
    !     . Seg2      . IntNomega     . Seg3
    !     .  ^        .  ^            .  ^
    !     .  |        .  |            .  |
    !     .Seg1+1     .Seg3+1         . Seg2+1
    !
    self%IntOmega( 1        :  Seg1 ) = xx            + LenY * CI
    self%IntOmega( Seg1 + 1 :  Seg2 ) = -self%jobr(1) + xy   * CI
    self%IntOmega( Seg2 + 1 :  Seg3 ) =  self%jobr(1) + xy   * CI
    do jc = 1 , SumN
       self%IntOmega( Seg3 + jc ) = ( 2 * jc - 1 ) * pi * T  * CI
    enddo
    !======================================================================
    do jc = 1 , Seg1
        self%IntOmegaWeight( jc ) = wx(jc) * CI/2._8/pi * f%FermionFunc( real(self%IntOmega(jc)),T)
    enddo
    !-----------------------------------------------------------
    do jc = 1 , self%jobi(6)
       self%IntOmegaWeight( Seg1 + jc ) =   wy(jc) * CI / 2._8 / pi &
                                           * f%FermionFunc(  self%IntOmega(Seg1 + jc) ,T)
       self%IntOmegaWeight( Seg2 + jc ) = - wy(jc) * CI / 2._8 / pi &
                                           * f%FermionFunc(  self%IntOmega(Seg2 + jc) ,T)
    enddo
    !-----------------------------------------------------------
    self%IntOmegaWeight( Seg3 + 1 : Seg4 ) = T
  endsubroutine



  subroutine GetSavingGandGrandPotetial(self)
    implicit none
    class(waldf),intent(inout)  :: self
    !----------------------------------------------
    type(CES)::solver
    TYPE(CEG)::G
    TYPE(Ham)::H
    integer::ns,jc1,jc2 ,NOmegaG  ,bg1,bg2                              !,jc,para(8),wtp
    ! complex*16::tempG(size(self%SavedG(:,1,1)),self%ns,self%ns)       !;class(FermOper),pointer::p
    complex*16,allocatable::tempG(:,:,:)

    bg1 = LBOUND(self%SavedG,1)
    bg2 = UBOUND(self%SavedG,1)
    NOmegaG = bg2 - bg1 + 1

    allocate(tempG(bg1:bg2,self%ns,self%ns))


    ns = self%ta%get_ns()
    !----------------------------------------------
    !  get H
    call H%Initialization(ns,self%getprint())
    call H%StartAppendingInteraction()
    Call self%cptH%AppendLocalDataToHam(H)
    Call self%dH%AppendDataToHam(H)
    call h%EndAppendingInteraction()

    !-----------------------------------------------

    Call solver%Initialization( self%SP , self%Ta, H , self%getprint() , show_ = 0)
    Call solver%SynchronizeWithHamiltonian()

    self%OmegaPri = solver%GetGrandPotential()


    Call G%Initialization(solver,self%Gp,self%getprint() , show_ = 0)

    call G%GetGreenMatrix(spini=0,spinj=0,NOmega=NOmegaG,Omega=self%IntOmega(bg1:bg2)&
                    , GM=tempG)     ;  self%SavedG(:,1:ns,1:ns) = tempG
    call G%GetGreenMatrix(spini=1,spinj=1,NOmega=NOmegaG,Omega=self%IntOmega(bg1:bg2)&
                    , GM=tempG)     ;  self%SavedG(:,ns+1:2*ns,ns+1:2*ns)= tempG
    call G%GetGreenMatrix(spini=0,spinj=1,NOmega=NOmegaG,Omega=self%IntOmega(bg1:bg2)&
                    , GM=tempG)     ;  self%SavedG(:,1:ns,ns+1:2*ns)  = tempG

    forall(jc1=bg1:bg2)
       self%SavedG(jc1,ns+1:2*ns,1:ns) = transpose( conjg( self%SavedG(jc1,1:ns,ns+1:2*ns) ) )
    endforall

   deallocate(tempG)
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
  complex*16 function FuncF(self,Gpri)
    IMPLICIT NONE
    class(waldf),intent(inout)  :: self
    complex*16,intent(in)::Gpri(self%ns*2,self%ns*2)
    !------------------------------
    TYPE(nummethod)::f
    integer::jck ,jc , symm ,spin ,i1,i2
    complex*16::temp(self%ns*2,self%ns*2),Vk(self%ns*2,self%ns*2)
    complex*16::temp1 ,TempS2(self%ns,self%ns)
    real*8::k(3),w
    !--------- initiate ------------
    FuncF = (0._8,0._8)
    !----------get G----------------

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
       temp = temp - matmul( Vk , Gpri )
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

  real*8 function ReccorectI(self)
    IMPLICIT NONE
    class(waldf),intent(inout)  :: self
    !--------------------------------------
    real*8,parameter::pi = 3.141592653589793238462643383279_8
    real*8::y1,y2,x1,x2
    complex*16::Gtemp(  size(self%SavedG(1,:,1))  , size(self%SavedG(1,1,:))    )

    ReccorectI = 0._8
    select case(self%jobi(1))
    case(2)
       x1 = abs(self%IntOmega(-1))
       x2 = abs(self%IntOmega( 0))
       Gtemp = self%SavedG(-1,:,:)
       y1 = real(FuncF(self,Gtemp) * self%IntOmega(-1))
       Gtemp = self%SavedG( 0,:,:)
       y2 = real(FuncF(self,Gtemp) * self%IntOmega( 0))
       ReccorectI = (x1*y2-y1*x2)/(x1-x2)
       ReccorectI = ReccorectI / 4  !  degree = 90.   ->  d/2Pi
       ReccorectI = ReccorectI * 2  !  upper and lower plane
    endselect                         !;WRITE(*,*)ReccorectI,y1,y2
  endfunction

  real*8 function GetI(self)
    IMPLICIT NONE
    class(waldf),intent(inout)  :: self
    !--------------------------------------
    integer::jc
    complex*16::ret
    complex*16::Gtemp(  size(self%SavedG(1,:,1))  , size(self%SavedG(1,1,:))    )

    ret = (0._8,0._8)
    do jc = 1 , self%IntNOmega
       Gtemp = self%SavedG(jc,:,:)
       ret = ret + FuncF(self,Gtemp ) * self%IntOmegaWeight(jc)
    enddo
    GetI = 2._8 * real(ret)  + ReccorectI(self)
  endfunction

  real*8 function GetLatticeOmegaPerSite(self)
    IMPLICIT NONE
    class(waldf),intent(inout)  :: self
    !--------------------------------------
    call self%CheckInitiatedOrStop()
    call GetSavingGandGrandPotetial(self)
    GetLatticeOmegaPerSite = self%OmegaPri - GetI(self)            ! ;write(*,*)self%OmegaPri, GetI(self),666
    GetLatticeOmegaPerSite = GetLatticeOmegaPerSite / self%ns
  endfunction

  real*8 function GetClusterOmegaPriPerSite(self)
    IMPLICIT NONE
    class(waldf),intent(inout)  :: self
    !--------------------------------------
    call self%CheckInitiatedOrStop()
    call GetSavingGandGrandPotetial(self)
    GetClusterOmegaPriPerSite = self%OmegaPri / self%ns
  endfunction


  real*8 function GetLatticeIPerSite(self)
    IMPLICIT NONE
    class(waldf),intent(inout)  :: self
    !--------------------------------------
    call self%CheckInitiatedOrStop()
    call GetSavingGandGrandPotetial(self)
    GetLatticeIPerSite = GetI(self)/ self%ns
  endfunction


  real*8 function GetOmegaValue(self,typeid)
    IMPLICIT NONE
    class(waldf),intent(inout)  :: self
    integer                     :: typeid
    !--------------------------------------
    select case(typeid)
    case(1)
      GetOmegaValue = self%GetLatticeOmegaPerSite()
    case(2)
      GetOmegaValue = self%GetClusterOmegaPriPerSite()
    case(3)
      GetOmegaValue = self%GetLatticeIPerSite()
    case default
      write(self%getprint(),*)"Unknow typeid=",typeid,'in GetOmegaValue@VCA_Waldfun'
      stop
    endselect
  endfunction



endmodule
