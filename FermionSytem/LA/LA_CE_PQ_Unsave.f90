

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  :  MODULE
! NAME  :  LA_CE_G_PQ_unsave
! OBJECT:  TYPE(LAPQUSV)
! USED  :  LA_Subspace
! DATE  :  2018-01-07
! AUTHOR:  hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            For a given lanczos space, calculate:
!
!                                    1
!           Q(ω) =    <G |A ─────────────────────── B|G >
!                            ω  - ( H - Eg ) * η
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
!                  [sub] Initialization(LaSub,A,B,eta,UseFrac,M_,OTH_,BZERO_,PRINT_,SHOW_)
!                        type(LASubSpace)::LaSub
!                        type(FermOper)  ::A,B (see DESCRIPTION)
!                        integer         ::eta (see DESCRIPTION)
!                        logical         ::UseFrac   . if true, use fractional expression.
!                                          ()
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


module LA_CE_G_PQ_unsave
  use LA_Subspace
  implicit none

  type::LAPQUSV
   private
   logical:: initiated = .false.
   class(LASubSpace),pointer :: LAsub        => null()
   class(table),pointer      :: ta           => null()
   type(FermOper)            :: A,B
   integer                   :: BSubId
   integer                   :: BSubD
   !------
   integer                   :: M            =  90
   logical                   :: oth          = .false.
   real*8                    :: bzero        = 0.000001_8
   integer                   :: eta
   logical                   :: UseFrac


   integer::print = 6
   INTEGER::SHOW  = 0
   !---------------------
 contains
   procedure,pass::Initialization
   final::Finalization

   procedure,pass::GetG
  endtype

  private::Initialization,UnInitialization,Finalization

  private::CheckBspaceSubId,GetSMA,GetG
contains


  subroutine Initialization(self,LaSub,A,B,eta,UseFrac,M_,OTH_,BZERO_,PRINT_,SHOW_)
     implicit none
     class(LAPQUSV),intent(inout)::self
     class(LASubSpace),INTENT(IN),TARGET    :: LaSub
     class(FermOper),intent(in)             :: A,B
     integer          ,intent(in)           :: eta
     logical          ,intent(in)           :: UseFrac
     integer          ,intent(in),optional  :: M_      ,PRINT_,SHOW_
     logical          ,intent(in),optional  :: OTH_
     real*8           ,intent(in),optional  :: bzero_
     !-----------------------------------------------------
     call UnInitialization(self)  ;  self%initiated = .true.
     self%LAsub => lasub
     self%A = A
     self%B = B
     self%eta = eta
     self%UseFrac = UseFrac
     self%ta => self%lasub%GetTablePoinnter()
     if (present(M_    )) self%M     = M_
     if (present(OTH_  )) self%OTH   = OTH_
     if (present(bzero_)) self%bzero = bzero_
     if (present(PRINT_)) self%print = print_
     if (present(SHOW_ )) self%SHOW  = SHOW_

     !------------------check subid of B space -------------
     self%BSubId = CheckBspaceSubId(self%Ta,self%LAsub%GetSubId(),self%B)
     if (self%BSubId.ne.-1)then
        self%BSubD = self%ta%get_sub_d(self%BSubId)
     endif
  endsubroutine

  subroutine UnInitialization(self)
    implicit none
    class(LAPQUSV),intent(inout)::self
    !-----------------------------------

  endsubroutine

  impure elemental subroutine Finalization(self)
    implicit none
    type(LAPQUSV),intent(inout)::self
    !-----------------------------------
    call UnInitialization(self)
  endsubroutine



  integer function CheckBspaceSubId(Ta,RecentSubId,B)
    implicit none
    class(table),intent(inout)::Ta
    integer,intent(in)::RecentSubId
    class(FermOper),intent(inout)::B
    !---------------------------------
    integer::jc ,sig
    INTEGER*8::si,so
    checkBspaceSubId = -1
    do jc = 1 , ta%get_sub_d(RecentSubId)
      si = ta%get_subindex_to_basis(RecentSubId,jc)
      call b%act(si,so,sig)
      if (So.ne.-1_8)then
        CheckBspaceSubId = ta%get_subid_from_basis(so)
        goto 999
      endif
    enddo
999 continue
  endfunction



  ! i = denegerate index   SubId is the dimention of Q space.
  subroutine GetSMA(self,i,Nomega,Omega,G)  ! R is operator and A is the energy vector.
    implicit none
    class(LAPQUSV),intent(inout)::self
    integer,intent(in)::i,Nomega
    complex*16,intent(in)::Omega(Nomega)
    complex*16,intent(out)::G(Nomega)
    !-----------------------------------
    class(slist),pointer::sl
    CLASS(TABLE),POINTER::TA
    CLASS(HAM),POINTER::h
    complex*16::SM(self%BSubD,self%M),AGS(self%BSubD),mags,gam,Mt(self%BSubD)!,temp
    complex*16,allocatable::am(:)
    real*8    ::A(SELF%M),B(SELF%M),Eg
    INTEGER   ::CutM,jcm,jcmt,Gd,KRM
    real*8,allocatable::Hi(:,:),Ei(:)
    logical::IsReal
    !-------------------
                                                                    !;write(*,*)"start"
    G = (0._8,0._8)
    Eg = self%LAsub%GetEg()
    Gd = self%LAsub%getD()
    IsReal = self%LAsub%IsSysReal()
    AGS = self%lasub%GetOptActState(self%B,i,self%BSubD)
    SM(:, 1) = AGS                                          !;WRITE(*,*)sum(abs(ags));write(*,*)"-------"
    TA => self%ta
    H  => self%lasub%GetHamPointer()
    allocate(sl)                                                   !;write(*,*)self%bzero

    KRM = min ( self%m  ,  self%BSubD  )

!---------------------------------------------------------------------------------------
!   set M
    ! call self%lasub%creat_krylov_space(&
    ! Ta =  TA  , H=h,sl=SL,&
    ! isreal=self%lasub%IsSysReal(),subid=self%BSubId,d=self%BSubD,&
    ! OTH1   = self%oth, OTH2=.false.,Bzero=self%bzero,&
    ! M=self%m,STATE_M=sm,A=A,B=B(2:self%M),Mcut=CutM,wtp=SELF%PRINT,SHOW=SELF%show)
!---------------------------------------------------------------------------------------
!  check max M as d
    call self%lasub%creat_krylov_space(&
    Ta =  TA  , H=h,sl=SL,&
    isreal=self%lasub%IsSysReal(),subid=self%BSubId,d=self%BSubD,&
    OTH1   = self%oth, OTH2=.false.,Bzero=self%bzero,&
    M=KRM,STATE_M=sm(:,1:KRM),A=A(1:KRM),B=B(2:KRM),Mcut=CutM,wtp=SELF%PRINT,SHOW=SELF%show)


    deallocate(sl)
    IF (CutM.eq.0) goto 999
    !------------------------------------
    !    if operator A = B , we can use A,B matrix directly   (recent result is incorrect)
    ! B(1) = 1._8
    ! if (self%UseFrac)then
    !   do jcm =   CutM , 2 ,-1
    !      G = B(jcm)**2 / (  A(jcm) - omega - Eg - G  )
    !   enddo
    !
    !   G = 1._8 /  (  A(1) - omega - Eg - G  )
    !   G = -G / self%eta
    !   goto 999
    ! endif
                                                        !  ;write(*,*)CutM,sum(abs(sm(:,8)));stop


                                                                           !;temp = 0._8
    allocate( Hi(CutM,CutM) ,Ei(CutM) )
    CALL self%lasub%Lan_Diag3Matrix(CutM ,A(1:CutM),B(1:CutM),Hi,Ei)
  !  write(*,*)Eg , Ei,6666
    allocate(am(gd))                                                       !;write(*,*)cutm
    do jcm = 1 , CutM
      !-------------------------
      !  get mt
      mt = (0._8,0._8)
      do jcmt = 1 , CutM
        mt = mt +  Hi(jcmt,jcm) * SM(:,jcmt)
      enddo
      mags = GetDotProduct(self%BSubD,IsReal, mt , AGS )
      !-------get AM ---------------------------------------
      am = self%lasub%GetOperActOnState(self%A,self%BSubId,self%BSubD,Gd,mt)

      !-----------------------------------------------------
      gam  = self%lasub%GetGsProduct(i,am)
      !---------

      G = G + mags * gam / (Omega - ( Ei(jcm) - Eg ) * self%eta  )   !;temp=temp+mags * gam
    enddo
    deallocate(am)
    !------------------------------------
    deallocate( Hi          ,Ei       )                                !;write(*,*)temp
999 continue                                !;write(*,*)"here"
  endsubroutine

  SUBROUTINE GetG(self,Nomega,Omega,G)
    implicit none
    class(LAPQUSV),intent(inout)::self
    integer,intent(in)::Nomega
    complex*16,intent(in) ::Omega(Nomega)
    complex*16,intent(out)::G(Nomega)
    !----------------------------------------------
    integer::jcde
    complex*16::dG(Nomega)
    G = (0._8,0._8)                                          !;write(*,*)"start"
    if (self%initiated)then
      if (self%BSubId.ne.-1)then
          do jcde = 1 , self%LAsub%GetDe()
             call GetSMA(self,jcde,Nomega,Omega,dG)
             G = G + dG
          enddo
      endif
      G = G / self%LAsub%GetDe()
    else
      write(self%print,*)"ERROR: Q method is used before Initialization";stop
    endif                                                        ! ;write(*,*)"here"
  endsubroutine




endmodule
