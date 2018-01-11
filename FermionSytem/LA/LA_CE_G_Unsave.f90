

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  :  MODULE
! NAME  :  LA_CE_G_unsave
! OBJECT:  TYPE(LACEGUSV)
! USED  :  LA_Subspace,LA_CE_G_PQ_unsave
! DATE  :  2018-01-07
! AUTHOR:  hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            For a given lanczos space, calculate:
!
!                                   1                                 1
!       G(i,j,ω) =    <G |c(i) ─────────── cd(j)|G > + <G |cd(j) ─────────── c(i)|G >
!                              ω - H + Eg                         ω + H - Eg
!
!          Different degeneracte states have been considered.
!
! STANDARD:
!            *CALL I
!
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  [sub] Initialization(La,i,spini,j,spinj,M_,OTH_,BZERO_,PRINT_,SHOW_)
!                        type(LASubSpace)::La
!                        integer::i,spini,j,spinj
!
!                        integer::M_     = 90
!                        logical::OTH_   = .true.
!                        real*8::bzero_  = 1.e-6
!                        integer::print_ = 6
!                        integer::show_  = 0
!
!
!
!
! avalable gets:
!                   [sub] GetG(Nomega,Omega,G)
!                         integer::Nomega
!                         complex*16::Omega(Nomega)
!                         complex*16::G(Nomega)
! avalable is :
!                  ![fun] i
! others      :
!                  [sub] S
!
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module LA_CE_G_unsave
  use LA_CE_G_PQ_unsave
  use LA_Subspace
  implicit none



  type::LACEGUSV
    private
    logical:: initiated = .false.
    class(LASubSpace),pointer ::LAsub        => null()
    class(table),pointer      ::ta           => null()
    integer::i,j,spini,spinj
    !------
    type(LAPQUSV)::p
    type(LAPQUSV)::q
    !------
    integer                   ::M            =  90
    logical                   ::oth          = .false.
    real*8                    ::bzero        = 0.000001_8
    integer                   ::eta

    integer::print = 6
    INTEGER::SHOW  = 0
    !---------------------
  contains
    procedure,pass::Initialization
    final::Finalization


    procedure,pass::GetG
  endtype


  private::Initialization,UnInitialization,Finalization

  private::GetG

contains


  subroutine Initialization(self,La,i,spini,j,spinj,M_,OTH_,BZERO_,PRINT_,SHOW_)
    implicit none
    class(LACEGUSV),intent(inout)         :: self
    class(LASubSpace),INTENT(IN),TARGET   :: La
    integer,intent(in)                    :: i , spini , j , spinj
    integer          ,intent(in),optional :: M_      ,PRINT_,SHOW_
    logical          ,intent(in),optional :: OTH_
    real*8           ,intent(in),optional :: bzero_
    !------------------------------------------------
    type(FermOper)::c,cd
    integer::para(8),ns
    logical::UseFrac

    call UnInitialization(self)  ;  self%initiated = .true.

    self%LAsub => la
    self%i     = i
    self%spini = spini
    self%j     = j
    self%spinj = spinj
    if (present(M_    )) self%M     = M_
    if (present(OTH_  )) self%OTH   = OTH_
    if (present(bzero_)) self%bzero = bzero_
    if (present(PRINT_)) self%print = print_
    if (present(SHOW_ )) self%SHOW  = SHOW_

    if ((i.eq.j) .and. (spini.eq.spinj) )then
      UseFrac = .true.
    else
      UseFrac = .false.
    endif

    !------------------------------ set cd -------------------------------
    ns = self%LAsub%GetNs()
    para = 0 ; para(1) = self%i ; para(2) = self%spini
    call c%Initialization( ns=ns,optid=1,para=para,print_=self%print)
    para = 0 ; para(1) = self%j ; para(2) = self%spinj
    call cd%Initialization(ns=ns,optid=2,para=para,print_=self%print)
    !----------------------------- set PQ --------------------------------
    call self%P%Initialization(LaSub=self%LAsub,A=c,B=cd,eta= 1,UseFrac=UseFrac,&
         M_=self%M,OTH_=self%oth,BZERO_=self%bzero,PRINT_=self%print,SHOW_=self%show)
    call self%Q%Initialization(LaSub=self%LAsub,A=cd,B=c,eta=-1,UseFrac=UseFrac,&
         M_=self%M,OTH_=self%oth,BZERO_=self%bzero,PRINT_=self%print,SHOW_=self%show)

  endsubroutine


  subroutine UnInitialization(self)
    implicit none
    class(LACEGUSV),intent(inout)         :: self
    !-------------------------------------------------
    if (self%initiated) then   ;  self%initiated = .false.

    endif
  endsubroutine

  impure elemental subroutine Finalization(self)
    implicit none
    type(LACEGUSV),intent(inout) :: self
    !-------------------------------------------------
    call  UnInitialization(self)
  endsubroutine









  subroutine GetG(self,Nomega,Omega,G)
    implicit none
    class(LACEGUSV),intent(inout)         :: self
    integer,intent(in)::Nomega
    complex*16,intent(in) ::Omega(Nomega)
    complex*16,intent(out)::G(Nomega)
    !-------------------------------------------------
    complex*16::Gp(Nomega)

    call self%LAsub%SynchronizeWithHamiltonian()

    call self%p%getg(Nomega,Omega,Gp)   ! ;write(*,*)Gp(1),"Gp"
    call self%q%getg(Nomega,Omega,G )
                                      !  ;write(*,*)G(1),"G"
    G = Gp + G                        !  ;write(*,*)G(1),"G+GP"
  endsubroutine


endmodule
