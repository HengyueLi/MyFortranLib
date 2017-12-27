

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  :  MODULE
! NAME  :  LA_CE_G_PQ
! OBJECT:  TYPE(LAPQ)
! USED  :  LA_Subspace
! DATE  :  2017-12-05
! AUTHOR:  hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            Calculate the Q and P value which will be used in Green'function calculation of Canonical Emsenbel.
!            This is general for both single particle and two particle Green's function.
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
!                  [sub] Initialization(LaSub,Ni,OPT,PQSubid,Asign,M_,OTH_,BZERO_,PRINT_,SHOW_)
!                        type(LASubSpace)::LaSub
!                        integer::Ni
!                        type(FermOper)::OPT(Ni)
!                        integer::PQSubid      after opt act on the GS, the subid of the output state.
!                        integer::Asign        For Q and P , the A differt by -1 (to check the expression)
!                                              +1 for P and -1 for Q
!                        integer::M_   = 90
!                        logical::OTH_ = .true.
!                        real*8::bzero_
!                        integer::print_
!
!                   [sub] SynchronizeQP()
!                         check and Synchronize P&Q matrix
!
!
!
!
!,,GetA
! avalable gets:
!                   [fun] GetMcut()
!                         integer  Cutoff of M
!                   [sub] GetPQ(PQ)
!                         complex*16::pq(self%M,self%Ni,self%De)
!                   [sub] GetA(A)
!                         complex*16::A(self%M)
!                   [sub] GetGpq(Tran,Nomega,Omega,G)
!                         logical::Tran
!                         Integer::Nomega
!                         complex*16::Omega(Nomega)
!                         complex*16::G(Nomega,Ni,Ni)
!                         return  (1/z)* Q^\dagger 1/(z-A)Q     if Tran = .false.
!                         return  (1/z)*(Q^\dagger 1/(z-A)Q)^T  if Tran = .true.
! avalable is :
!                  ![fun] i
! others      :
!                  [sub] S
!
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


module LA_CE_G_PQ
  use LA_Subspace
   implicit none


   type::LAPQ
      private
      logical                   ::initiated    = .false.
      class(LASubSpace),pointer ::LAsub        => null()
      class(table),pointer      ::Ta           => null()
      class(Ham)  ,pointer      ::H            => null()
      integer                   ::M            =  90
      logical                   ::oth          = .false.
      real*8                    ::bzero        = 1.e-7
      integer                   ::Ni
      TYPE(FermOper),ALLOCATABLE::OPT(:)
      integer                   ::PQSubid
      complex*16,allocatable    ::PQ(:,:,:)  !pq(m,i,De)
      real*8,allocatable        ::A(:)       !A(m)
      integer                   ::Asign
      integer                   ::CutM
      integer                   ::De
      integer::EigenId = -146519
      !-----------------------------------------------------------
      integer::print = 6
      INTEGER::SHOW  = 0

    contains
      procedure,pass::Initialization
      final::Finalization

      procedure,pass::SynchronizeQP
      procedure,pass::GetMcut
      procedure,pass::GetPQ
      procedure,pass::GetA
      procedure,pass::GetGpq
   endtype


   private::Initialization,UnInitialization,Finalization
   private::CheckDeallocateAll

   private::ReNewPQDirectly,SynchronizeQP
   private::GetMcut,GetPQ,GetA,GetGpq
   private::get_AGS


 contains


   subroutine Initialization(self,LaSub,Ni,OPT,PQSubid,Asign,M_,OTH_,BZERO_,PRINT_,SHOW_)
     implicit none
     class(LAPQ)      ,INTENT(INOUT)        :: self
     class(LASubSpace),INTENT(IN),TARGET    :: LaSub
     INTEGER          ,INTENT(IN)           :: nI
     TYPE(FermOper)   ,intent(in)           :: OPT(ni)
     integer          ,intent(in)           :: PQSubid ,Asign
     integer          ,intent(in),optional  :: M_      ,PRINT_,SHOW_
     logical          ,intent(in),optional  :: OTH_
     real*8           ,intent(in),optional  :: bzero_
     !--------------------------------------------
     call UnInitialization(self)  ;  self%initiated = .true.

     self%LaSub   => Lasub
     self%ta      => self%Lasub%GetTablePoinnter()
     self%H       => self%Lasub%GetHamPointer()
     self%Ni      =  Ni
     self%PQSubid =  PQSubid
     self%Asign   =  Asign
     allocate(self%OPT,source = OPT)
     if (present(M_    )) self%M     = M_
     if (present(OTH_  )) self%OTH   = OTH_
     if (present(bzero_)) self%bzero = bzero_
     if (present(PRINT_)) self%print = print_
     if (present(SHOW_ )) self%SHOW  = SHOW_

   endsubroutine


   subroutine UnInitialization(self)
     implicit none
     class(LAPQ),INTENT(INOUT)::self
     !--------------------------------------
     if (self%initiated)then ; self%initiated = .false.
       call CheckDeallocateAll(self)
       if (allocated(self%OPT) ) deallocate(self%OPT )
     endif
   endsubroutine

   impure elemental subroutine Finalization(self)
     implicit none
     type(LAPQ),INTENT(INOUT)::self
     !--------------------------------------
     call UnInitialization(self)
   endsubroutine


   subroutine CheckDeallocateAll(self)
     implicit none
     class(LAPQ),INTENT(INOUT)::self
     !--------------------------------------
     if (allocated(self%pq)  ) deallocate(self%pq  )
     if (allocated(self%a  ) ) deallocate(self%a   )
   endsubroutine



   !--------------------------------------------
   !renew P and Q according to LanSub without checking anything.
   Subroutine ReNewPQDirectly(self)
     use Statelist20171204   ,only: slist     => ListStru
     implicit none
     class(LAPQ),intent(inout)::self
     !---------------------------------------
     integer::De,dpq,jc,jc1,JC2,jcj
     complex*16,allocatable::SM(:,:),AGS(:)
     class(slist),pointer::sl
     real*8::A(self%M),B(self%M),Eg
     real*8,allocatable::H(:,:),E(:)
    !  TYPE(nummethod)::f
     logical::isreal

     De     = self%LaSub%GetDe()
     isreal = self%lasub%IsSysReal()
     Eg     = self%lasub%GetEg()
     self%de = de
     call CheckDeallocateAll(self)
     allocate(  self%pq(   self%M, self%Ni, De )  )  ; self%pq = (0._8,0._8)
     allocate(  self%a(    self%M              )  )  ; self%a  = 0._8
     !---------------------check zero -------------------
     if (self%PQSubid.eq.-1)then
        self%CutM = 0
        goto 999
     endif
     dpq = self%Ta%get_sub_d(self%PQSubid)
     !---------------------Creat SM---------------------
     allocate(  SM( dpq , self%M )  )
     allocate(  AGS(dpq)            )

    !  call self%Lasub%GetRadomState(dpq,SM(:,1))
     call get_AGS(SELF,dpq,SM(:,1))
     !---------------------Creat SM---------------------

     allocate(sl)
     call self%Lasub%creat_krylov_space(Ta = self%Ta,  H=self%H,  sl=sl,&
            isreal = isreal  ,  subid=self%PQSubid, d=dpq,&
            OTH1   = self%oth, OTH2=.false.,  Bzero=self%bzero,M=self%m,&
            STATE_M= SM ,A=A,B=B(2:self%M) ,  Mcut=self%cutM,wtp=self%print,SHOW=SELF%SHOW)
    deallocate(sl)
    if (self%CutM == 0) goto 999
    allocate(H(self%cutM,self%cutM),E(self%cutM))
     !----------------------diag------------------------
    ! call f%ED_tridiagonal_real(self%cutM,A(1:self%cutM),B(1:self%cutM),H,E)
     CALL self%Lasub%Lan_Diag3Matrix(self%cutM,A(1:self%cutM),B(1:self%cutM),H,E)
    !---------------------Pmi = <m|Oi|g>----------------
    ! <m| = \sum_j H_{jm} <SM_j|
    self%pq = (0._8,0._8)
    do jc = 1 , De
       do jc1 = 1 , self%ni   !;;write(*,*)self%opt(jc1)%get_para()
         AGS = self%Lasub%GetOptActState(self%OPT(jc1),jc,dpq) !;write(*,*)dpq!;pause
         do jc2 = 1 , self%CUTM
           do jcj = 1, self%CUTM
             self%pq(jc2,jc1,De) = self%pq(jc2,jc1,De) &
              + H(jcj,jc2) * GetDotProduct(dpq,isreal,SM(:,jcj),AGS)
           enddo
         end do
       enddo
    enddo
    do jc = 1 , self%cutm
      self%a(jc) = E(jc) - Eg
    enddo
    !----------------------------------------------------
    self%a = self%a * self%Asign
    !----------------------------------------------------
    deallocate(SM,ags,h,e)
999 continue
   endsubroutine



   subroutine get_AGS(SELF,dags,AGS)
     class(LAPQ),intent(inout)::self
     integer,intent(in)::dags
     complex*16,intent(out)::AGS(dags)
     !---------------------------------------
     integer::jc,jcde
     complex*16::ags1(dags)
     real*8::Rfac


    !  call self%Lasub%GetRadomState(dags,AGS)

     !----------------------------------------
     AGS = (0._8,0._8)
    !  do jc = 1 , self%Ni
    do jc = 1,self%Ni
       do jcde = 1 , self%LAsub%GetDe()
          ags1 = self%LAsub%GetOptActState(self%OPT(jc),jcde,dags)
          ! call random_number(Rfac)
          AGS  = AGS + ags1 !* ( Rfac + 0.5_8 )
       enddo
     enddo
   endsubroutine


   subroutine SynchronizeQP(self)
     implicit none
     class(LAPQ),intent(inout)::self
     !---------------------------------------
     call self%LAsub%SynchronizeWithHamiltonian()
     if (self%EigenId .ne. self%LAsub%GetEigenId())THEN
        call ReNewPQDirectly(self)
        self%EigenId = self%LAsub%GetEigenId()
     ENDIF
   endsubroutine


   integer function GetMcut(self)
     implicit none
     class(LAPQ),intent(inout)::self
     !---------------------------------------
     GetMcut = self%cutM
   endfunction

   subroutine GetPQ(self,pq)
     implicit none
     class(LAPQ),intent(inout)::self
     complex*16,intent(out)::pq(self%m,self%ni,self%de)
     !---------------------------------------
     pq = self%pq
   endsubroutine

   subroutine GetA(self,A)
     implicit none
     class(LAPQ),intent(inout)::self
     real*8,intent(out)::A(self%m)
     !---------------------------------------
     A = self%A
   endsubroutine


   subroutine GetGpq(self,Tran,Nomega,Omega,G)
     implicit none
     class(LAPQ),intent(inout)::self
     logical,intent(in)::Tran
     integer,intent(in)::Nomega
     complex*16,intent(in)::Omega(Nomega)
     complex*16,intent(out)::G(Nomega, self%Ni,self%Ni)
     !---------------------------------------
     integer::jcDe,JCO,jcm,jci,jcj
     complex*16::QdQ(self%Ni,self%Ni),Am(self%CutM)
     G = (0._8,0._8)
     !---
     do jcDe = 1 , self%De
        do jci = 1 , self%ni ; do jcj = 1 , self%ni
           do jcm = 1 , self%CutM
              G(:,jci,jcj) = G(:,jci,jcj) +&
                  conjg(self%pq(jcm,jci,jcde)) * self%pq(jcm,jcj,jcde) / (omega - self%A(jcm))

           enddo
        enddo                ; enddo
     enddo
     if (Tran) then
       do jco = 1 , nomega
         G(jco,:,:) = transpose(G(jco,:,:))
       enddo
     endif

     G = G / self%De

   endsubroutine


endmodule
