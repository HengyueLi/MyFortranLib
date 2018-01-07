


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  :  MODULE
! NAME  :  LA_GCE_G
! OBJECT:  TYPE(LAGCEG)
! USED  :  LA_CE_G,LA_GCESpace
! DATE  :  2018-01-07
! AUTHOR:  hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            Calculate the single particle Green's function of Ground Canonical Emsenbel
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
!                  [sub] Initialization(La,Nopt ,M_,OTH_,BZERO_,PRINT_,SHOW_)
!                        class(LA_GCE)::LaSub
!                        integer::Nopt(:,:)    ! Nopt(1,:) = site   Nopt(2,:) = spin
!                        integer::M_     = 90
!                        logical::OTH_   = true
!                        real*8 ::bzero_ = 1.e-7
!                        integer::PRINT_ = 6
!
!                  [sub] SynchronizeQP()
!                        check and Synchronize P&Q matrix
!---------------------------------------------------
! avalable gets:
!                   [sub] GetGmatrix(Nomega,Omega,G)
!                         integer::Nomega
!                         complex*16::Omega(Nomega),G(Nomega,ni,ni)
!                         NOTICE!!!:  The index of matrix G is not site index.
!                                     The are the same as the input Nopt.
!
! avalable is :
!                  ![fun] i
! others      :
!                  [sub] S
!
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module LA_GCE_G
  use LA_CE_G
  use LA_GCESpace
  implicit none


  type::LAGCEG
    private
    logical                :: initiated = .false.
    integer                :: Nsuib
    class(LA_GCE),pointer  :: FullLa    => null()
    TYPE(LACEG),allocatable:: subG(:)
    integer                :: ni
    integer                :: M         =  90
    logical                :: oth       =  .true.
    real*8                 :: bzero     =  1.e-7!1.e-7

    !-----
    integer::print
    INTEGER::SHOW = 0

  contains
    procedure,pass::Initialization
    final:: Finalization

    procedure,pass::SynchronizeQP
    procedure,pass::GetGmatrix
  endtype




  private::Initialization,UnInitialization,Finalization

  private::SynchronizeQP
  private::GetGmatrix

contains


  subroutine Initialization(self,La,Nopt,M_,OTH_,BZERO_,PRINT_,SHOW_)
    implicit none
    class(LAGCEG),intent(inout)         :: self
    class(LA_GCE),intent(inout),target  :: La
    integer,intent(in)                  :: Nopt(:,:)
    integer,intent(in),optional         :: M_
    logical,intent(in),optional         :: OTH_
    real*8,intent(in),optional          :: bzero_
    integer,intent(in),optional         :: PRINT_,SHOW_
    !--------------------------------
    integer::jcsub
    class(LASubSpace),pointer :: LaSub

    call UnInitialization(self) ;  self%initiated = .true.


    self%ni  =  size(Nopt(1,:))
    if(present(M_    )) self%m     = M_
    if(present(OTH_  )) self%oth   = oth_
    if(present(bzero_)) self%bzero = bzero_
    if(present(print_)) self%print = print_
    if(present(SHOW_ )) self%SHOW  = SHOW_
    self%FullLa => la

    allocate(self%subG( La%GetNsub()  ))
    do jcsub = 1 , La%GetNsub()
      LaSub => La%GetSubSpacePointer(jcsub)
      call self%subG(jcsub)%Initialization(LaSub = LaSub,Nopt = Nopt ,&
               M_ = self%M ,OTH_ = self%oth,BZERO_=self%bzero,&
               PRINT_=self%print,SHOW_=self%show )
    enddo

  endsubroutine

  subroutine UnInitialization(self)
    implicit none
    class(LAGCEG),intent(inout)::self
    !--------------------------------
    if (self%initiated)then  ; self%initiated = .false.
       deallocate(self%subG)
    endif
  endsubroutine

  impure elemental subroutine Finalization(self)
    implicit none
    type(LAGCEG),intent(inout)::self
    !--------------------------------
    call UnInitialization(self)
  endsubroutine



  subroutine SynchronizeQP(self)
    implicit none
    class(LAGCEG),intent(inout)::self
    !--------------------------------
    integer::jc,GSsubid
    do jc = 1 , self%FullLa%GetNsubGS()
       GSsubid = self%FullLa%GetSubIdGs(jc)
       call self%subG(GSsubid)%SynchronizeQP()
    enddo
  endsubroutine




   subroutine GetGmatrix(self,Nomega,Omega,G)
     implicit none
     class(LAGCEG),intent(inout)::self
     integer,intent(in)::Nomega
     complex*16,intent(in)::OMega(Nomega)
     complex*16,intent(out)::G(Nomega,self%ni,self%ni)
     !--------------------------------
     integer::jc ,GSsubid
     complex*16::dG(Nomega,self%ni,self%ni)
     G = (0._8,0._8)

     !call SynchronizeQP(self)      contained

     do jc = 1 , self%FullLa%GetNsubGS()
        GSsubid = self%FullLa%GetSubIdGs(jc)
        call self%subG(GSsubid)%GetGmatrix(Nomega,Omega,dG)
        G = G + dG * self%FullLa%GetSubDe(GSsubid)
     enddo
     G = G / self%FullLa%GetDe()

   endsubroutine








endmodule
