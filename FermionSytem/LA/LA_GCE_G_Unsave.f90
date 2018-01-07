

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  :  MODULE
! NAME  :  LA_GCE_G_Unsave
! OBJECT:  TYPE(LAGCEGUSV)
! USED  :  LA_CE_G_unsave,LA_GCESpace
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
!                  [sub] Initialization(LA,i,spini,j,spinj ,M_,OTH_,BZERO_,PRINT_,SHOW_)
!                        class(LA_GCE)::LaSub
!                        integer::i,spini,j,spinj
!                        integer::M_     = 90
!                        logical::OTH_   = true
!                        real*8 ::bzero_ = 1.e-7
!                        integer::PRINT_ = 6
!
!---------------------------------------------------
! avalable gets:
!                   [sub] GetG(Nomega,Omega,G)
!                         integer::Nomega
!                         complex*16::Omega(Nomega),G(Nomega)
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



module LA_GCE_G_Unsave
  use LA_CE_G_unsave
  use LA_GCESpace
  implicit none


  type::LAGCEGUSV
    private
    logical                   :: initiated = .false.
    class(LA_GCE),pointer     :: FullLa    => null()
    integer::Nsub
    TYPE(LACEGUSV),allocatable:: subG(:)

    integer                   :: M         =  90
    logical                   :: oth       =  .true.
    real*8                    :: bzero     =  1.e-7!1.e-7
    !-----
    integer::print
    INTEGER::SHOW = 0

  contains
    procedure,pass::Initialization
    final::Finalization

    procedure,pass::GetG
  endtype



  private::Initialization,UnInitialization,Finalization

  private::GetG


contains

  subroutine Initialization(self,La,i,spini,j,spinj ,M_,OTH_,BZERO_,PRINT_,SHOW_)
    implicit none
    class(LAGCEGUSV),intent(inout)      :: self
    class(LA_GCE),intent(in),target     :: La
    integer,intent(in)                  :: i,spini,j,spinj
    integer,intent(in),optional         :: M_
    logical,intent(in),optional         :: OTH_
    real*8,intent(in),optional          :: bzero_
    integer,intent(in),optional         :: PRINT_,SHOW_
    !--------------------------------
    integer::jc
    class(LASubSpace),pointer :: LaSub

    call UnInitialization(self) ;  self%initiated = .true.
    if(present(M_    )) self%m     = M_
    if(present(OTH_  )) self%oth   = oth_
    if(present(bzero_)) self%bzero = bzero_
    if(present(print_)) self%print = print_
    if(present(SHOW_ )) self%SHOW  = SHOW_
    self%FullLa => la


    !------------------
    self%Nsub  = self%FullLa%GetNsub()
    allocate(  self%subG(self%nsub)  )

    do jc = 1, self%Nsub
       LaSub => self%FullLa%GetSubSpacePointer(jc)
       call self%subG(jc)%Initialization(La=Lasub,i=i,spini=spini,j=j,spinj=spinj,&
           M_=self%m,OTH_=self%oth,BZERO_=self%bzero,PRINT_=self%print,SHOW_=self%show)
    enddo

  endsubroutine


  Subroutine UnInitialization(self)
    implicit none
    class(LAGCEGUSV),intent(inout)::self
    !------------------------------------
    if (self%initiated)then ; self%initiated = .false.
       deallocate( self%subG  )
    endif
  endsubroutine


  Subroutine Finalization(self)
    implicit none
    type(LAGCEGUSV),intent(inout)::self
    !------------------------------------
    call UnInitialization(self)
  endsubroutine



  subroutine GetG(self,Nomega,Omega,G)
    implicit none
    class(LAGCEGUSV),intent(inout)::self
    integer,intent(in)::Nomega
    complex*16,intent(in) ::Omega(Nomega)
    complex*16,intent(out)::G(Nomega)
    !-------------------------------------------------
    integer::jc,sid
    complex*16::dG(Nomega)
    G = (0._8,0._8)                                 !;write(*,*)"start"
    do jc = 1 , self%FullLa%GetNsubGS()
       sid = self%FullLa%GetSubIdGs(jc)
       call self%subG(sid)%getG(Nomega,Omega,dG)
       G = G + dG * self%FullLa%GetSubDe(sid)
    enddo                                       !;write(*,*)sum(G),self%FullLa%GetDe()
    G = G / self%FullLa%GetDe()                    ! ;write(*,*)"here"
  endsubroutine




endmodule
