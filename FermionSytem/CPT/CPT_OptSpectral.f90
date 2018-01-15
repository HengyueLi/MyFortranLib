


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : CPT_OptSpectral
! OBJECT: TYPE(cptopt1)
! USED  : CodeObject , class_integration , GCE_CPT_para , CE_Green , LaLatticeH , CreateKspace , phyfunc
! DATE  : 2018-01-05
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            Calculate optical spectral. First we define the optial conductivity.
!            The regular part of optical conductivity is
!                          e^2 * pi
!            sigma(w) = 	───────── P(w)          ( here we can chose e = 1 )
!                           N  * w                                                      1
!                                 where P(w) is spectral of JJ correlation B: P(w) = - ─── * Im( B(w+i0^+)  )
!                                                                                       π
!            B( w ) = << j | j >>_w
!
!            Definition of Current operator j = i \sum_{ij} t_{ij} ( r_i - r_j ) c^+_i c_j
!
!            In cluster representation we can use the velocity operator. (iv_q^{ab})
!
!            j = i \sum_{q,a,b} v^{ab}_q \psi^{q+}_a \psi^{q}_b
!
!           v_q = \sum_R v_q(R) here R represent the nearest connected cluster (R=0 is also included).
!
!                 v_q(R) = - t_{0a,b} exp(iqR) (R - (ra-rb) )
!
!
!           (  notation : s -> spin   )
!           Consider the spin, we can write   j = j (s1,s2) = t_ij c^+_{1s1} c_{2s2}
!
!           B ( w ) = sum_{s1,s2,s3,s4}  B(s1,s2,s3,s4)     ,  s1,s2,s3,s4 will be the input
!
!
! STANDARD:
!
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                   [sub] Initialization(CPTpara,Gpara,spin4,NOmegaP,OmegaP1,OmegaP2,EtaP,nk3,print_,show_)
!                          class(GCECPTpara),target,intent(in)   :: CPTpara
!                          class(GreenPara),intent(in)           :: Gpara
!                          integer,intent(in)                    :: spin4(4),NomegaP
!                          real*8,intent(in)                     :: OmegaP1,OmegaP2,EtaP
!                          integer,intent(in)                    :: nk3(3)
!                          integer,intent(in),optional           :: print_,show_
!
!
! avalable gets:
!                   [sub]
!
!
!                   [sub] GetB0XY(self,Omega_cpx,x,y,B0xy)
!                          complex*16,intent(in)        :: Omega_cpx(:)
!                          integer,intent(in)           :: x,y
!                          complex*16,intent(out)        :: B0xy(:)
!                          calculate the bubble approximation of jj correlation.
!
!                   [sub] GetOpticalSigmaXYTofile(x,y,Nomega,o1,o2,eta,file)
!                          integer,intent(in)::x,y,Nomega
!                          real*8,intent(in)::o1,o2,eta
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



module CPT_OptSpectral
  use CodeObject
  use class_integration
  use GCE_CPT_para
  use CE_Green
  use LaLatticeH
  use CreateKspace
  use phyfunc
  implicit none


  type,extends(object)::cptopt1
    private
    integer::spin(4)


    integer:: NomegaP
    real*8 :: etaP
    complex*16,allocatable :: Apri1(:,:,:)  !(ns,ns,Nop)     A41
    complex*16,allocatable :: Apri2(:,:,:)  !(ns,ns,Nop)     A23
    real*8,allocatable     :: OmegaP(:)
    real*8,allocatable     :: OW(:)
    TYPE(Kspace)           :: k



    integer                   :: ns
    class(GCECPTpara),pointer :: CPTPARA => null()
    Type(GreenPara)           :: Gpara
    class(LH),pointer         :: cpt     => null()


  contains
    procedure,pass::Initialization
    final::Finalization

    procedure,pass::GetB0XY
    procedure,pass::GetOpticalSigmaXYTofile
  endtype


  private::Initialization,Finalization

  private::SetSpectral,SetA
  private::GetB0XY,GetOpticalSigmaXYTofile


contains


  subroutine Initialization(self,CPTpara,Gpara,spin4,NOmegaP,OmegaP1,OmegaP2,EtaP,nk3,print_,show_)
    implicit none
    class(cptopt1),intent(inout)          :: self
    class(GCECPTpara),target,intent(in)   :: CPTpara
    class(GreenPara),intent(in)           :: Gpara
    integer,intent(in)                    :: spin4(4),NomegaP
    real*8,intent(in)                     :: OmegaP1,OmegaP2,EtaP
    integer,intent(in)                    :: nk3(3)
    integer,intent(in),optional           :: print_,show_
    !---------------------------------------------------------------------------
    TYPE(Integration)::itg


    Call Finalization(self) ; call self%SetInitiated(.true.)

    if (present(print_)) call self%setprint(print_)
    if (present(show_ )) call self%setshow(show_ )

    self%spin    = spin4
    self%NomegaP = NomegaP

    self%CPTPARA => CPTpara
    self%Gpara   =  Gpara
    self%cpt     => self%CPTPARA%GetCPTHpointer()
    self%ns      =  self%cpt%GetNs()
    self%etaP    =  EtaP

    allocate(  self%Apri1(self%ns,self%ns,self%NomegaP) )
    allocate(  self%Apri2(self%ns,self%ns,self%NomegaP) )
    allocate(  self%OmegaP(self%NomegaP)                )
    allocate(  self%Ow(self%NomegaP)                    )

    call itg%Initialization( N = self%NomegaP,typecase = 1 )
    call itg%get_xw(xs = OmegaP1,xe = OmegaP2, x = self%OmegaP ,w = self%ow )
    call self%k%Initialization(a=self%cpt%GetVbasis(),n=nk3,meshtype=1,print_=self%getprint(),show_=0)

    call SetSpectral(self)

  endsubroutine


  subroutine Finalization(self)
    implicit none
    type(cptopt1),intent(inout):: self
    !----------------------------------
    if (self%isinitiated())then
       deallocate(  self%Apri1  , self%Apri2  ,  self%OmegaP  ,self%OW  )
       call self%SetInitiated( .false. )
    endif
  endsubroutine


  subroutine SetSpectral(self)
    implicit none
    class(cptopt1),intent(inout) :: self
    !--------------------------------------------
    call SetA(self,spin1=self%spin(4),spin2=self%spin(1),A=self%Apri1)
    if (  (self%spin(2)==self%spin(4))  .and.  (self%spin(1)==self%spin(3)) )then
      self%Apri2 = self%Apri1
    else
       call SetA(self,spin1=self%spin(2),spin2=self%spin(3),A=self%Apri2)
     endif
  endsubroutine


  subroutine SetA(self,spin1,spin2,A)
    implicit none
    class(cptopt1),intent(inout) :: self
    integer,intent(in)::spin1,spin2
    complex*16,intent(inout)::A(self%ns,self%ns,self%NomegaP)
    !---------------------------------
    complex*16::Omega(2*self%NomegaP),G(2*self%NomegaP,self%ns,self%ns)
    TYPE(CEG)::GetG
    integer::Nomega

    Nomega = 2*self%NomegaP
    Omega(1:self%NomegaP       ) = self%OmegaP - self%etaP*(0._8,1._8)
    Omega(self%NomegaP+1:Nomega) = self%OmegaP + self%etaP*(0._8,1._8)

    call GetG%Initialization(solver=self%CPTPARA%GetSolverPointer(),GP=self%Gpara,&
                             PRINT_=self%getprint(),SHOW_=0)
    call getG%GetGreenMatrix(spini=spin1,spinj=spin2,NOmega=Nomega,Omega=Omega,GM=G)

    A = ( G(1:self%NomegaP,:,:) - G(self%NomegaP+1:Nomega,:,:) ) /(0._8,2._8)/3.141592653589793238462643383279_8

  endsubroutine


  subroutine GetB0XY(self,Omega_cpx,x,y,B0xy)
    implicit none
    class(cptopt1),intent(inout) :: self
    complex*16,intent(in)        :: Omega_cpx(:)
    integer,intent(in)           :: x,y
    complex*16,intent(out)        :: B0xy(:)
    !--------------------------------------------------
    integer::jop1,jop2,jo,jck,print,Nomega
    complex*16::V12(self%ns,self%ns),V34(self%ns,self%ns)
    complex*16::AQ1(self%ns,self%ns),sumw
    logical::ToPrint
    real*8::q(3),wk,ow1,tempr
    TYPE(phyfun)::f

    Nomega = size(Omega_cpx)

    tempr = self%CPTPARA%GetTemperature()

    print = self%getprint()
    if (self%getshow()>0)then
      ToPrint = .true.
    else
      ToPrint = .false.
    endif

    B0xy = ( 0._8 , 0._8 )

    do jck = 1 , self%k%getnk()

      if ( ToPrint )then
         write(print,*)"Summation in K space:", int( (jck*1._8 /  self%k%getnk() * 100)), "%"
      endif

       q  = self%k%getk(jck)
       wk = self%k%getw(jck)
       V12 = self%cpt%GetVelocityMatrix(q=q,spini=self%spin(1),spinj=self%spin(2),xy=x)
       if (  (self%spin(1)==self%spin(3)) .and. (self%spin(2)==self%spin(4))    )then
         v34=v12
       else
         V34 = self%cpt%GetVelocityMatrix(q=q,spini=self%spin(3),spinj=self%spin(4),xy=x)
       endif
       do jop1 = 1 , self%NomegaP
         !---------------------------------------
         AQ1 = matmul(  self%Apri1(:,:,jop1) ,  v12  )
         AQ1 = matmul(  v34                  ,  AQ1  )
         !---------------------------------------
         do jop2 = 1 , self%NomegaP
           !-------------------------------------
            sumw = gettr(self%ns,AQ1,self%Apri2(:,:,jop2)) * wk * self%OW(jop1)* self%OW(jop2)  &
            * (f%FermionFunc(self%OmegaP(jop1),tempr) - f%FermionFunc(self%OmegaP(jop2),tempr) )
            B0xy = B0xy - sumw / ( Omega_cpx + self%OmegaP(jop1) -self%OmegaP(jop2)  )
         enddo
       enddo
    enddo
  contains
    complex*16 function gettr(n,a,b)
      integer::n
      complex*16::a(n,n),b(n,n)
      !--------------------------
      integer::jc1,jc2
      gettr = (0._8,0._8)
      do jc1 =1 , n
        do jc2 = 1 , n
           gettr = gettr + a(jc1,jc2) * b(jc2,jc1)
        enddo
      enddo
    endfunction
  endsubroutine


  subroutine GetOpticalSigmaXYTofile(self,x,y,Nomega,o1,o2,eta,file)
    implicit none
    class(cptopt1),intent(inout) :: self
    integer,intent(in)::x,y,Nomega
    real*8,intent(in)::o1,o2,eta
    character(*),intent(in)::file
    !----------------------------------------------------------
    complex*16::omega(Nomega),Bo(Nomega)
    integer::jc

    do jc = 1 , Nomega
      omega(jc) = o1 + (o2-o1)/Nomega * (jc-1) + eta*(0._8,1._8)
    enddo

    call GetB0XY(self,omega,x,y,Bo)

    open(1234,file=trim(adjustl(file)))
    do jc = 1 , Nomega
      write(1234,*)real(omega(jc)),-imag(Bo(jc))/3.141592653589793238462643383279_8
    enddo
    close(1234)


  endsubroutine













endmodule
