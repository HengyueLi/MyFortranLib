


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : CPT_OptSpectral
! OBJECT: TYPE(cptopt1)
! USED  : CodeObject , class_integration , GCE_CPT_para , CE_Green , LaLatticeH , CreateKspace , phyfunc , class_numerical_method,functionalsubs
! DATE  : 2018-01-16
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
!           B ( w ) = sum_{s1,s2,s3,s4}  B(s1,s2,s3,s4,w)     ,  s1,s2,s3,s4 will be the input
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    NOTICE:  This model can only offer part of optical conductivity. !!!!!!
!             See the representation above:
!             B is the real optical conductivity which this model calculate B(s1,s2,s3,s4,w) only.
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!                   [sub] Initialization(CPTpara,spin4,NOmegaP,OmegaP1,OmegaP2,EtaP,nk3,print_,show_)
!                          class(GCECPTpara),target,intent(in)   :: CPTpara
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
!                            here the factor pi/N is multiplied. (That's to say, sigma(w) is obtained,see detals above.)
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
  use class_numerical_method
  use functionalsubs
  implicit none


  type,extends(object)::cptopt1
    private
    integer::spin(4)


    integer:: NomegaP
    real*8 :: etaP
    complex*16,allocatable :: G41InR(:,:,:)  !(ns,ns,Nop)     Retarded G
    complex*16,allocatable :: G41InA(:,:,:)  !(ns,ns,Nop)     Retarded G
    complex*16,allocatable :: G23InR(:,:,:)  !(ns,ns,Nop)     Retarded G
    complex*16,allocatable :: G23InA(:,:,:)  !(ns,ns,Nop)     Retarded G
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

  private::SetInverseG,SetGInverseSpinDepd,GetSpectralFromGinverse
  private::GetB0XY,GetOpticalSigmaXYTofile


contains


  subroutine Initialization(self,CPTpara,spin4,NOmegaP,OmegaP1,OmegaP2,EtaP,nk3,print_,show_)
    implicit none
    class(cptopt1),intent(inout)          :: self
    class(GCECPTpara),target,intent(in)   :: CPTpara
    integer,intent(in)                    :: spin4(4),NomegaP
    real*8,intent(in)                     :: OmegaP1,OmegaP2,EtaP
    integer,intent(in)                    :: nk3(3)
    integer,intent(in),optional           :: print_,show_
    !---------------------------------------------------------------------------
    TYPE(Integration)::itg
    real*8::lca(3,3)


    Call Finalization(self) ; call self%SetInitiated(.true.)

    if (present(print_)) call self%setprint(print_)
    if (present(show_ )) call self%setshow(show_ )

    self%spin    = spin4
    self%NomegaP = NomegaP

    self%CPTPARA => CPTpara
    self%Gpara   =  self%CPTPARA%getgp()
    self%cpt     => self%CPTPARA%GetCPTHpointer()
    self%ns      =  self%cpt%GetNs()
    self%etaP    =  EtaP

    allocate(  self%G41InR(self%ns,self%ns,self%NomegaP) )
    allocate(  self%G41InA(self%ns,self%ns,self%NomegaP) )
    allocate(  self%G23InR(self%ns,self%ns,self%NomegaP) )
    allocate(  self%G23InA(self%ns,self%ns,self%NomegaP) )

    allocate(  self%OmegaP(self%NomegaP)                )
    allocate(  self%Ow(self%NomegaP)                    )

    call itg%Initialization( N = self%NomegaP,typecase = 1 )
    call itg%get_xw(xs = OmegaP1,xe = OmegaP2, x = self%OmegaP ,w = self%ow )
    lca = self%cpt%GetVbasis()
    call self%k%Initialization(a=lca,n=nk3,meshtype=1,print_=self%getprint(),show_=0)

    call SetInverseG(self)

  endsubroutine


  subroutine Finalization(self)
    implicit none
    type(cptopt1),intent(inout):: self
    !----------------------------------
    if (self%isinitiated())then
       deallocate(  self%G41InR,self%G23InR, self%G41InA,self%G23InA , self%OmegaP ,self%OW  )
       call self%SetInitiated( .false. )
    endif
  endsubroutine


  subroutine SetInverseG(self)
    implicit none
    class(cptopt1),intent(inout) :: self
    !--------------------------------------------
    call SetGInverseSpinDepd(self,spin1=self%spin(4),spin2=self%spin(1),GR=self%G41InR,Ga=self%G41InA)
    if (  (self%spin(2)==self%spin(4))  .and.  (self%spin(1)==self%spin(3)) )then
       self%G23InA = self%G41InA
       self%G23InR = self%G41InR
    else
      call SetGInverseSpinDepd(self,spin1=self%spin(2),spin2=self%spin(3),GR=self%G23InR,Ga=self%G23InA)
    endif
  endsubroutine


  subroutine SetGInverseSpinDepd(self,spin1,spin2,GR,GA)
    implicit none
    class(cptopt1),intent(inout) :: self
    integer,intent(in)::spin1,spin2
    complex*16,intent(inout)::Gr(self%ns,self%ns,self%NomegaP),GA(self%ns,self%ns,self%NomegaP)
    !---------------------------------
    complex*16::Omega(2*self%NomegaP),G(2*self%NomegaP,self%ns,self%ns)
    TYPE(CEG)::GetG
    integer::Nomega,jco
    TYPE(nummethod)::f

    Nomega = 2*self%NomegaP
    Omega(1:self%NomegaP       ) = self%OmegaP - self%etaP*(0._8,1._8)
    Omega(self%NomegaP+1:Nomega) = self%OmegaP + self%etaP*(0._8,1._8)

    call GetG%Initialization(solver=self%CPTPARA%GetSolverPointer(),GP=self%Gpara,&
                             PRINT_=self%getprint(),SHOW_=0)
    call getG%GetGreenMatrix(spini=spin1,spinj=spin2,NOmega=Nomega,Omega=Omega,GM=G)

    do jco = 1 , self%NomegaP
      GA(:,:,jco) =  G(jco,             :,:)   ;  call f%MatrixInverse(self%ns,GA(:,:,jco))
      GR(:,:,jco) =  G(jco+self%NomegaP,:,:)   ;  call f%MatrixInverse(self%ns,GR(:,:,jco))
    enddo

  endsubroutine


  SUBROUTINE GetSpectralFromGinverse(self,q,SPINI,SPINJ,GiR,GiA,A)
    implicit none
    class(cptopt1),intent(inout) :: self
    real*8,intent(in)::q(3)
    INTEGER,INTENT(IN)::spini,spinj
    complex*16,intent(in )::Gir(self%ns,self%ns,self%NomegaP),GiA(self%ns,self%ns,self%NomegaP)
    complex*16,intent(out)::A(self%ns,self%ns,self%NomegaP)
    !---------------------------------------------------------------------------
    complex*16::GR(self%ns,self%ns),GA(self%ns,self%ns)
    TYPE(nummethod)::f
    integer::jc
    do jc = 1 , self%NomegaP
       GR = GIR(:,:,JC) - SELF%cpt%GetTqPrimatrix(q,spini,spinj)
       GA = GIA(:,:,JC) - SELF%cpt%GetTqPrimatrix(q,spini,spinj)
       call f%MatrixInverse(self%ns,GR)
       call f%MatrixInverse(self%ns,GA)
       A(:,:,jc) = (GA-GR)/(0._8,2._8)/3.141592653589793238462643383279_8
    enddo
  endsubroutine


  subroutine GetB0XY(self,Omega_cpx,x,y,B0xy)
    implicit none
    class(cptopt1),intent(inout) :: self
    complex*16,intent(in)        :: Omega_cpx(:)
    integer,intent(in)           :: x,y
    complex*16,intent(out)        :: B0xy(:)
    !--------------------------------------------------
    integer::jop1,jop2,jo,jck,print,Nomega,totNk
    complex*16::V12(self%ns,self%ns),V34(self%ns,self%ns)
    complex*16::AQ1(self%ns,self%ns),sumw
    COMPLEX*16::A41(self%ns,self%ns,self%NomegaP),A23(self%ns,self%ns,self%NomegaP)
    logical::ToPrint
    real*8::q(3),wk,ow1,tempr,Qualty
    character(128)::PrintQualty,tempC,PrintProgress,tempPP
    TYPE(phyfun)::f
    TYPE(funcsubs)::fprint



    !---------------------------------------------------------------------------------------
    ! calculate integeration quality
    !  the larger the better
    Qualty = self%etaP/(( self%OmegaP(size(self%OmegaP)) - self%OmegaP(1) ) /  self%NomegaP)
    write(PrintQualty,*)int(Qualty)
    write(tempC,*)int((Qualty - int(Qualty))*100)
    PrintQualty = trim(adjustl(PrintQualty))//"."//trim(adjustl(tempC))
    PrintQualty = "Integration quality factor ~ ["//trim(adjustl(PrintQualty))&
                    //"] (the larger the better, 1 is bad)"
    !-------------------------------------------------------------------------------------



    Nomega = size(Omega_cpx)

    tempr = self%CPTPARA%GetTemperature()

    print = self%getprint()
    if (self%getshow()>0)then
      ToPrint = .true.
    else
      ToPrint = .false.
    endif

    B0xy = ( 0._8 , 0._8 )

    if (ToPrint)then
      write(self%getprint(),*)"Start to caclulate jj correlation in Bubble level."
      write(self%getprint(),*)PrintQualty
    endif

    totNk = self%k%getnk()

    do jck = 1 , totNk

      if ( ToPrint )then
        write(PrintProgress,*)"Summation in K space:"
        write(tempPP,*)int( (jck*1._8 /  self%k%getnk() * 100))
        write(PrintProgress,*)trim(adjustl(PrintProgress))//trim(adjustl(tempPP))//"%"
        call fprint%IntegerReportProgress(print,totNk,jck,0.1_8,PrintProgress)
      endif

       q  = self%k%getk(jck)
       wk = self%k%getw(jck)

       CALL GetSpectralFromGinverse(self,q,self%spin(4),self%spin(1),self%G41InR,self%G41InA,A41)
       V12 = self%cpt%GetVelocityMatrix(q=q,spini=self%spin(1),spinj=self%spin(2),xy=x)
       V34 = self%cpt%GetVelocityMatrix(q=q,spini=self%spin(3),spinj=self%spin(4),xy=y)
       if (  (self%spin(1)==self%spin(3)) .and. (self%spin(2)==self%spin(4))    )then
         A23 = A41
       else
         CALL GetSpectralFromGinverse(self,q,self%spin(2),self%spin(3),self%G23InR,self%G23InA,A23)
       endif
       do jop1 = 1 , self%NomegaP
         !---------------------------------------
         AQ1 = matmul(  A41(:,:,jop1) ,  v12  )
         AQ1 = matmul(  v34           ,  AQ1  )
         !---------------------------------------
         do jop2 = 1 , self%NomegaP
           !-------------------------------------
            sumw = gettr(self%ns,AQ1,A23(:,:,jop2)) * wk * self%OW(jop1)* self%OW(jop2)  &
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
    real*8,parameter::pi = 3.141592653589793238462643383279_8
    real*8,parameter::Tiny=1E-10
    complex*16::omega(Nomega),Bo(Nomega)
    integer::jc


    do jc = 1 , Nomega
      omega(jc) = o1 + (o2-o1)/Nomega * (jc-1) + eta*(0._8,1._8)
    enddo



    call GetB0XY(self,omega,x,y,Bo)


    open(1234,file=trim(adjustl(file)))
    do jc = 1 , Nomega
      write(1234,*)real(omega(jc)),-imag(Bo(jc))/(real(OMega(jc)) + Tiny )/self%ns     
    enddo
    close(1234)

  endsubroutine













endmodule
