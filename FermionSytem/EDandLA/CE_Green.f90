

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  :  MODULE
! NAME  :  CE_Green
! OBJECT:  TYPE(CEG)
! USED  :  CEsolver,CodeObject,LA_GCE_G_Unsave,ED_GCE_G
! DATE  :  2017-12-29
! AUTHOR:  hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            Canonical Emsenbel Sytem.   Calculate the single particle Green's function.
!
! STANDARD:
!            *CALL
!
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  [sub] Initialization(solver,PRINT_,SHOW_,M_,OTH_,BZERO_)
!                          class(CES),intent(in )         ,target :: solver
!                          integer   ,intent(in ),optional        :: PRINT_,SHOW_,M_
!                          logical   ,intent(in ),optional        :: oth_
!                          real*8    ,intent(in ),optional        :: Bzero_
!
!
! avalable gets:
!                   [fun] GetG(i,spini,j,spinj,Nomega,Omega,G)
!                          integer,intent(in )::i,spini,j,spinj,Nomega
!                          complex*16,intent(in )::Omega(Nomega)
!                          complex*16,intent(out)::G(Nomega)
!
!                   [sub] GetGreenMatrix(spini,spinj,NOmega,Omega,GM)
!                          integer,intent(in)::spini,spinj,Nomega
!                          complex*16,intent(in)::Omega(Nomega)
!                          complex*16,intent(out)::GM(Nomega,self%ns,self%ns)
!
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

module CE_Green
  use CodeObject
  use CEsolver
  use LA_GCE_G_Unsave
  use ED_GCE_G
  implicit none



  type::GreenPara
    ! lanczos G use
    integer :: M     = 90
    logical :: oth   = .true.
    real*8  :: bzero = 1.e-7
  endtype

  type,extends(object)::CEG

    integer::ns
    class(CES),pointer::solver => null()
    type(GreenPara)::GP
  contains
    procedure,pass::Initialization

    procedure,pass::GetG
    procedure,pass::GetGreenMatrix
  endtype


  private::Initialization
  private::getg
  private::GetGreenMatrix

contains


  subroutine Initialization(self,solver,GP,PRINT_,SHOW_)
    implicit none
    class(CEG),intent(out)                 :: self
    class(CES),intent(in )         ,target :: solver
    class(GreenPara),intent(in)            :: GP
    integer   ,intent(in ),optional        :: PRINT_,SHOW_
    !------------------------------------
    call self%SetInitiated(.true.)


    self%solver    => solver
    self%Gp        = Gp
    if( present( print_   )  )  call self%setprint(print_)
    if( present( show_    )  )  call self%setshow(show_  )

    self%ns        = self%solver%getns()

  endsubroutine


  subroutine GetG(self,i,spini,j,spinj,Nomega,Omega,G)
    implicit none
    class(CEG),intent(inout):: self
    integer,intent(in )::i,spini,j,spinj,Nomega
    complex*16,intent(in )::Omega(Nomega)
    complex*16,intent(out)::G(Nomega)
    !--------------------------------------------------
    type(LAGCEGUSV)   :: LAG
    type(EDGCEGreenf) :: EDG
    class(ED_GCE),pointer::ED
    class(LA_GCE),pointer::la
    call self%CheckInitiatedOrStop()
    select case(self%solver%GetSolverType())
    case("ED")
      ed => self%solver%GetEdPointer()
      call edg%Initialization(ED,self%getprint())
      call edg%GetG(i,spini,j,spinj,Nomega,Omega,G)
    case("LA")
      la => self%solver%GetlaPointer()
      call lag%Initialization(LA,i,spini,j,spinj ,self%gp%m,self%gp%oth,self%gp%bzero,&
            self%getprint(),self%getshow())
      call lag%GetG(Nomega,Omega,G)
    endselect
  endsubroutine

  subroutine GetGreenMatrix(self,spini,spinj,NOmega,Omega,GM)
    implicit none
    class(CEG),intent(inout):: self
    integer,intent(in)::spini,spinj,Nomega
    complex*16,intent(in)::Omega(Nomega)
    complex*16,intent(out)::GM(Nomega,self%ns,self%ns)
    !---------------------------------
    integer::jc1,jc2

    if (spini==spinj)then
      do jc1 =1 , self%ns
        do jc2 = jc1 , self%ns
          call GetG(self,jc1,spini,jc2,spinj,Nomega,Omega,GM(:,jc1,jc2))
        enddo
      enddo
      do jc1 =1 , self%ns
        do jc2 = 1 , jc1 - 1
          GM(:,jc1,jc2) = GM(:,jc2,jc1)
        enddo
      enddo
    else
      do jc1 =1 , self%ns
        do jc2 = 1 , self%ns
            call GetG(self,jc1,spini,jc2,spinj,Nomega,Omega,GM(:,jc1,jc2))
        enddo
      enddo
    endif
  endsubroutine


! GetNs


endmodule
