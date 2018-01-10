

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : LA_Subspace
! OBJECT: TYPE(LASubSpace)
! USED  : fermion_table,FermionHamiltonian,basic_math_functions
! DATE  : 2018-01-07
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            consider a subspace in Lanczos. This subspace is identified by a subid which is correpsonding to that in table.
!
! STANDARD:
!            *CALL Initialization(Ta,IsReal,H,subid,PRINT_,show_,PreE_,PreDe_,bzero_,M_,oth_)
!             call diagonalization()
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  [sub] Initialization(Ta,IsReal,H,subid,PRINT_,SHOW_,PreE_,PreDe_,bzero_,M_,oth_)
!                        type(table)::Ta  (will be pointed and be read only)
!                        logical    ::IsReal ! if H is a real problem.
!                        type(Ham)  ::H , will be  pointed.
!                        integer    ::subid
!                        integer::M_ = 30   ! Krylov space
!                        real*8::PreE_ = 1.e-14   precision of dE in Lanczos
!                        real*8::PreDe_ = 1.e-6  precision of checking degeneracy.  (percentage)
!                        real*8::bzero_ = 1.E-6_8 cutoff precision of b
!                        logical::oth_ = .true.  Gram－Schmidt in creating Krylov space.
!
!                        Ta shoule be initiated before this subroutine since some information would be used.
!                  ![sub] No
!
! avalable gets:
!                   [sub] SynchronizeWithHamiltonian()
!                         Check if H is diagonalized, and if not, diagonalize it.
!
!                  [fun] GetD()
!                        get the dimension of the subspace
!
!                  [fun] GetEg()
!                         get ground state energy.
!
!                  [fun] GetE(i)
!                         Get E_i
!
!                  [fun] GetDe()
!                         degeneracy
!
!                  [fun] GetSubId()
!                        integer::
!
!                  [fun] GetNs()
!                        integer::
!
!                  [fun] GetOptActState(opt,i,d)
!                        complex*16::GetOptActState(d)
!                        opt | i> = | GetOptActState(d) >
!                        the dimension, d, of output state should be known first
!                        act on the i-th GS in the subspace
!
!                  [fun] GetGsProduct(i,S)
!                        return <GS(i)|S>    complex*16
!                        size(S) should be the same as GS
!
!                  [fun] GetOperActOnState(opert,subidin,din,dout,Si)
!                        |out(dout)> = opert |Si(subidin,din)>
!
!                  [fun] GetOperateProduct(A)
!                        TYPE(FermOper)::A
!                        return complex*16::<A> = 1/D * sum_{all GS} <GS|A|GS>
!
!
! avalable is :
!                  ![fun] i
! others      :
!                  ![sub] p
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



 module LanczosState
   implicit none
    type::State
      real*8::E
      complex*16,allocatable::s(:)
    endtype
 endmodule

 module Statelist20171204
     use LanczosState , data => State
     INCLUDE "../../ListStructure/ListStructure.finc"
 endmodule



module LA_Subspace
    use fermion_table
    use FermionHamiltonian
    use LanczosState        ,only: state
    use Statelist20171204   ,only: slist     => ListStru
    implicit none
    !----------------------------

    type::LASubSpace
      private

      logical::  initiated     = .false.
      class(table),pointer::Ta => null()
      logical::  IsReal        = .false.
      class(Ham),pointer::  H  => null()
      integer:: subid !
      integer:: d     ! dimension of state in this subspace
      ! integer:: De
      TYpe(slist)::state
      !------------------
      integer::EigenId = -9135
      ! integer::
      real*8 ::Eg
      !------------------------------------------
      integer::LanM     = 30
      real*8 ::bzero    = 0.000001
      integer::Maxitera = 10000
      logical::Oth1     = .true.   !  orthonormal to previous states or not.
      real*8 ::LanPre   = 1.e-14   !
      real*8 ::DegPre   = 1.e-6    !  percentage, for checking the degenerate state.
      !------------------------------------------
      integer::show  = 0
      integer::print = 6

    contains
      procedure,pass::Initialization
      final::Finalization

      procedure,pass::SynchronizeWithHamiltonian
      procedure,pass::GetEg
      procedure,pass::GetDe
      procedure,pass::GetE
      procedure,pass::GetOptActState
      procedure,pass::GetTablePoinnter
      procedure,pass::GetHamPointer
      procedure,pass::GetSubId
      procedure,pass::GetNs


      procedure,pass::GetRadomState
      procedure,pass::IsSysReal
      procedure,nopass::creat_krylov_space
      procedure,pass::GetEigenId
      procedure,NOpass::Lan_Diag3Matrix
      procedure,pass::GetGsProduct
      ! procedure,nopass::oper_Act_On_State
      procedure,pass::GetOperActOnState
      procedure,pass::GetD
      procedure,pass::GetOperateProduct
    endtype


    private::Initialization,UnInitialization,Finalization


    private::oper_Act_On_State , H_Act_On_State
    private::ProductOperator   , ProductH
    ! private::GetDotProduct
    private::GetStateOthorgnal , GetStateOthorgnal2 , GetStateOthorgnalToList
    private::LANCZOS_ITERATE , ScanOneGroundState  ,GetRadomState
    private::Lan_Diag3Matrix
    private::Diagonalization
    private::CheckDiffSmall
    private::SynchronizeWithHamiltonian
    private::GetEg,GetDe,GetE
    private::GetOptActState
    private::GetTablePoinnter
    private::IsSysReal
    private::GetHamPointer
    private::creat_krylov_space
    PRIVATE::GetEigenId
    private::GetSubId
    private::GetGsProduct
    private::GetNs
    private::GetOperActOnState
    private::GetD
    private::GetOperateProduct

  contains

    subroutine Initialization(self,Ta,IsReal,H,subid,PRINT_,SHOW_,PreE_,PreDe_,bzero_,M_,oth_)
      implicit none
      class(LASubSpace),intent(inout)::self
      class(table),intent(inout),target::Ta
      logical,intent(in)::IsReal
      class(ham),intent(inout),target::H
      integer,intent(in)::subid
      integer,intent(in),optional::PRINT_,SHOW_,M_
      real*8,intent(in),optional::PreE_,PreDe_,bzero_
      logical,intent(in),optional::oth_
      !--------------------------------
      call UnInitialization(self)   ;  self%initiated = .true.
      if (present(PRINT_)) self%print = PRINT_
      if (present(SHOW_))  self%show  = SHOW_
      if (present(PreE_)) self%LanPre = PreE_
      if (present(PreDe_)) self%DegPre= PreDe_
      if (present(bzero_)) self%bzero = bzero_
      if (present(M_))     self%lanM  = M_
      if (present(oth_))   self%Oth1  = oth_
      !------------------------------------------------------------------------------
      !  check SubId allowed
         if ((ta%get_nsub().ge.subid) .and.(subid.ge.1) )then
         else
            write(self%print,*)"ERROR: input subid is wrong in LASubSpace."
            write(self%print,*)"This may be caused by incorrect setting of symmetry.";stop
         endif
      !-------------------------------------------------------------------------------
      self%Ta => Ta
      self%IsReal = IsReal
      self%H   => H
      self%subid = subid
      self%d   = self%ta%get_sub_d(subid)


    endsubroutine

     Impure elemental subroutine Finalization(self)
      implicit none
      type(LASubSpace),intent(inout)::self
      !--------------------------------------
      call UnInitialization(self)
    endsubroutine

    subroutine UnInitialization(self)
      implicit none
      class(LASubSpace),intent(inout)::self
      !-------------------------------------------
      if (self%initiated)then
         self%initiated = .false.
         self%Ta => null()
         self%H  => null()
      endif
    endsubroutine












        complex*16 function GetDotProduct(d,isreal,sl,sr)
           implicit none
           integer,intent(in)::d
           logical,intent(in)::isreal
           complex*16,intent(in)::sl(d),sr(d)
           !-------------------------
           GetDotProduct = DOT_PRODUCT(  sl   ,  sr     )
        endfunction

        !  S = S - <Sothor|S> * |Sothor>
        subroutine GetStateOthorgnal(d,IsReal,S,Sothor)
          implicit none
          integer::d
          logical,intent(in)::IsReal
          complex*16,intent(inout)::S(d)
          complex*16,intent(in)   ::Sothor(d)
          !------------------------------
          S = S - GetDotProduct(d,isreal,Sothor , S) * Sothor
        endsubroutine


        ! Sothor(d,n)
        !  S = S - \Sum_i^n <Sothor_i|S> * |Sothor_i>
        subroutine GetStateOthorgnal2(n,d,IsReal,S,Sothor)
          implicit none
          integer,intent(in)::n,d
          logical,intent(in)::IsReal
          complex*16,intent(inout)::S(d)
          complex*16,intent(in)   ::Sothor(d,n)
          !------------------------------
          integer::jc
          do jc = 1 , n
            call GetStateOthorgnal(d,IsReal,S,Sothor(:,jc))
          enddo
        endsubroutine

        ! s = s - <sl|s> * sl
        subroutine GetStateOthorgnalToList(d,isreal,sl,s)
          implicit none
          integer,intent(in)::d
          logical,intent(in)::isreal
          class(slist),intent(inout)::sl
          complex*16,intent(inout)::s(d)
          !------------------------------------------
          integer::jc
          class(state),pointer::p
          do jc = 1 , sl%getlen()
             p => sl%GetDataPointer(jc)
             call GetStateOthorgnal(d,IsReal,S,p%s)
          enddo
        endsubroutine


    ! for a given satate (StateIn) in subspace subid
    ! |AddStateOut(d)> = |AddStateOut> + V * oper | StateIn >
    subroutine oper_Act_On_State(Ta,SubIdIn,oper,V,IsReal,din,dout,StateIn,AddStateOut)
      implicit none
      class(table),intent(inout)::Ta
      integer,intent(in)::SubIdIn
      class(FermOper),intent(inout)::oper
      complex*16,intent(in)::V
      logical,intent(in)::IsReal
      integer,intent(in)::din,dout
      complex*16,intent(in)::StateIn(din)
      complex*16,intent(out)::AddStateOut(dout)
      !----------------------------------------------------------
      integer*8::s8,s9
      integer::jc,sid,sign,so
      complex*16::vi
      sid = SubIdIn
      vi  = V
      do jc = 1, din
        s8 = ta%get_subindex_to_basis(sid,jc)
        call oper%act(s8,s9,sign)
        if (s9.ne.-1_8)then
            so = ta%get_basis_to_sub_index(s9)
            AddStateOut(so) =  AddStateOut(so) + Vi * StateIn(jc) * sign
        endif
      enddo
    endsubroutine


    !   return complex*16      V * < StateIn | oper | StateIn >
    complex*16 function ProductOperator(Ta,subid,oper,V,IsReal,d,StateIn)
      implicit none
      class(table),intent(inout)::Ta
      integer,intent(in)::subid
      class(FermOper),intent(inout)::oper
      complex*16,intent(in)::V
      logical,intent(in)::IsReal
      integer,intent(in)::d
      complex*16,intent(in)::StateIn(d)
      !----------------------------------------------------------
      complex*16::AddStateOut(d)
      AddStateOut = (0._8,0._8)
      call oper_Act_On_State(Ta,subid,oper,V,IsReal,d,d,StateIn,AddStateOut)
      ProductOperator =V * GetDotProduct(d,isreal,StateIn , AddStateOut)
    endfunction

    !   return complex*16      V * < StateIn |H| StateIn >
    complex*16 function ProductH(Ta,subid,H,IsReal,d,StateIn,wtp)
      implicit none
      class(table),intent(inout)::Ta
      integer,intent(in)::subid
      class(Ham),intent(inout)::H
      logical,intent(in)::IsReal
      integer,intent(in)::d
      complex*16,intent(in)::StateIn(d)
      integer,intent(in)::wtp
      !----------------------------------------------------------
      complex*16::AddStateOut(d)
      AddStateOut = (0._8,0._8)
      call H_Act_On_State(Ta,subid,H,IsReal,d,StateIn,AddStateOut,wtp)
      ProductH    = GetDotProduct(d,isreal,StateIn , AddStateOut )
    endfunction



    ! for a given satate (StateIn) in subspace subid
    ! |AddStateOut> = |AddStateOut> + H | StateIn >
    subroutine H_Act_On_State(Ta,subid,H,IsReal,d,StateIn,AddStateOut,wtp)
      implicit none
      class(table),intent(inout)::Ta
      integer,intent(in)::subid
      class(Ham),intent(inout)::H
      logical,intent(in)::IsReal
      integer,intent(in)::d
      complex*16,intent(in)::StateIn(d)
      complex*16::AddStateOut(d)
      integer,intent(in)::wtp
      !----------------------------------------------------------
      class(FermOper),pointer::fa
      complex*16::v
      integer::jc
      if (.not.H%Is_Usable())then
         write(wtp,*)"When using H_Act_On_State, H is not set yet.";stop
      endif
      do jc = 1 , H%GetOptN()
         v  =  h%GetOptV(jc)
         fa => h%GetOptact(jc)
         call  oper_Act_On_State(Ta,subid,fa,V,IsReal,d,d,StateIn,AddStateOut)
      enddo
    endsubroutine








    !-------------------------------------------------------------------------------------------------------------------------
    !STATE_M: STATE_M(:,1) should be set as the input guess vector.
    !         it can be none nonmormolized.Normorlization will be considered insite.
    !A      : can be any input
    !B      : can be any input
    !oth1   : if make each states othognal to the previous states in krylov space.
    !oth2   : if make each states othognal to ground state.
    subroutine creat_krylov_space(Ta,H,sl,isreal,subid,d,OTH1,OTH2,Bzero,M,STATE_M,A,B,Mcut,wtp,SHOW)
       implicit none
       class(table),intent(inout)::Ta
       class(Ham),intent(inout)::H
       class(slist),intent(inout)::sl
       logical,intent(in)   ::isreal
       integer,intent(in)   ::subid,d
       logical,intent(in)   ::OTH1,OTH2
       real*8,intent(in)    ::Bzero
       integer,intent(in)   ::M
       complex*16,intent(inout)::STATE_M(d,M)
       real*8,intent(inout)::A(m),B(2:m)
       integer,intent(inout)::Mcut
       integer,intent(in)::wtp,SHOW
       !---------------------------------------------------------------
       integer::jc  ;real*8::b1
       !----------------------------------------------------------------
       !   check M >=3
           if (.not.(M.gt.2))then
              write(wtp,*)"ERROR: Krylov space can not be created when dimension M<3";stop
           endif
       !----------------------------------------------------------------
       !   set Mcut
           Mcut = M
       !----------------------------------------------------------------
           a = 0._8  ; b = 0._8
       !----------------------------------------------------------------
       !  oth2
          if (oth2) call GetStateOthorgnalToList(d,isreal,sl,STATE_M(:,1))
       !----------------------------------------------------------------
       call INNERPRINT(1._8/M)
       !   normolize s1
           b1 = dsqrt(REAL( GetDotProduct(d,isreal,STATE_M(:,1) , STATE_M(:,1))))
           if (b1.le.Bzero) then
            !  write(wtp,*)"the input guess state is too small,Sqrt(<G|G>)<",Bzero
            !  write(wtp,*)"the calculation may be meanness"
             !-----action
            !  stop
             Mcut = 0
             goto 999
             !-----------
           endif
           STATE_M(:,1) = STATE_M(:,1) / b1
       !----------------------------------------------------------------
       !   Initialize other SM
          STATE_M(:,2:M) = (0._8,0._8)
       !----------------------------------------------------------------
       call INNERPRINT(2._8/M)
       !   calculation at jcm = 1
          call H_Act_On_State(Ta,subid,H,IsReal,d,STATE_M(:,1),STATE_M(:,2),wtp)
          a(1) = real( GetDotProduct(d,isreal,STATE_M(:,1) , STATE_M(:,2)) )
          STATE_M(:,2) = STATE_M(:,2) - a(1) * STATE_M(:,1)
          !  oth1
          if (oth1) call GetStateOthorgnal(d,IsReal,STATE_M(:,2),STATE_M(:,1))
          !  oth2
          if (oth2) call GetStateOthorgnalToList(d,isreal,sl,STATE_M(:,2))
          !---get b(2)
          b(2) = dsqrt(REAL( GetDotProduct(d,isreal,STATE_M(:,2) , STATE_M(:,2))))
          if (b(2).le.Bzero)then
              Mcut = 1
              goto 999
          endif
          !---set SM2
          STATE_M(:,2) = STATE_M(:,2) / b(2)
        !----------------------------------------------------------------
        ! Mi>=3
          DO jc=3,M
            call INNERPRINT(jc*1._8/M)
            !------------------------------------
            !   |sn> = H |sn-1>
              call H_Act_On_State(Ta,subid,H,IsReal,d,STATE_M(:,JC-1),STATE_M(:,JC),wtp)
            !------------------------------------
            !   set a(jc-1)
              a(jc-1) = real( GetDotProduct(d,isreal,STATE_M(:,jc-1) , STATE_M(:,jc)) )
            !-----------------------------------
            !   |sn> = |sn> - a(n-1) * |sn-1> - b(n-1)^2 *  |sn-2>
              STATE_M(:,JC) = STATE_M(:,JC) - a(jc-1) * STATE_M(:,JC-1) - b(jc-1) * STATE_M(:,JC-2)
            !-----------------------------------
            ! oth1
              if (oth1) call GetStateOthorgnal2(jc-1,d,IsReal,STATE_M(:,JC),STATE_M(:,1:jc-1))
            ! oth2
              if (oth2) call GetStateOthorgnalToList(d,isreal,sl,STATE_M(:,JC))
            !-----------------------------------
            !   set b
              b(jc) = dsqrt(real(  GetDotProduct(d,isreal,STATE_M(:,JC) , STATE_M(:,JC)) ))
            !-----------------------------------
            !   check b = 0
              if (b(jc).le.bzero) then
                Mcut = jc - 1
                !----------------------------
                !  show
                if (show.ge.4)then
                    write(wtp,*)"|b| approach to 0 when M=",jc
                    write(wtp,*)"Cut off^"
                endif
                !----------------------------
                goto 999
              endif
            !------------------------------------
            !  set |sn>
              STATE_M(:,JC)=STATE_M(:,JC) /B(JC)
            !------------------------------------
            ! set a(m)
            if (jc.eq.m) a(m) = real(ProductH(Ta,subid,H,IsReal,d,STATE_M(:,m),wtp))
          ENDDO
      999 continue

    CONTAINS
      SUBROUTINE INNERPRINT(p)
        implicit none
        real*8,intent(in)::p
        !---------------------------------------
        character(64)::pR
        character(16)::R,fi
        real*8::c
        if (show .lt.4) goto 999


        c = (p*100)
        if (c.lt.10)then
          write(fi,*)".."
        elseif (c.lt.100)then
          write(fi,*)"."
        else
          write(fi,*)""
        endif
        WRITE(R,"(1f5.1)")c
        R = trim(adjustl(fi))//trim(adjustl(R))
        pr = "Creating Krylov Space:..."//trim(adjustl(r))//"%"
        WRITE(WTP,*)trim(adjustl(pr))

999     continue
      ENDSUBROUTINE
    endsubroutine


    !=================================================================================
    !input old ground state,  output new ground state and the corresponding energy.
    !input GS can be non othognal
    !---------------------------------
    ! realturn :
    !           0  succesfull
    !           1  b for GS is 0.
    !---------------------------------------------------------------------------------
    integer function LANCZOS_ITERATE(self,GS,EG)
    use class_numerical_method
    IMPLICIT NONE
    class(LASubSpace),intent(inout)::self
    complex*16,intent(inout)::gs(self%d)
    real(8),intent(inout)::eg
    !-----------------------------------------------------
    real*8::A(self%LanM),B(1:self%LanM)
    real*8,allocatable::energy(:), eigenVr(:,:)
    complex*16::SM(self%d,self%LanM)
    INTEGER::Mcut,jc1
    type(nummethod)::mt

    !-----------------------set Sm1 -----------------------------------------------
    SM(:,1)=GS
    !-----------------------Get Krylov space --------------------------------------
    call creat_krylov_space(self%Ta,self%H,self%state,self%IsReal,self%subid,&
    self%d,self%OTH1,.true.,self%bzero,self%LanM,SM,A,B(2:self%LanM),Mcut,self%print,SELF%SHOW)
    !------------------------------------------------------------------------------
    ! the input <GS|GS> =0 or GS is othognal to other ground states
    ! it is obtained by b1 = 0
    if (Mcut.eq.0_8) then
      !  write(self%print,*)"The tested GS is 0 or orthognal to other GS."
      !  write(self%print,*)"this case is not well considered."
       LANCZOS_ITERATE = 1
       goto 999
    endif
    !-----------diag Krylov space--------
    allocate(  energy (Mcut     )  )
    allocate(  eigenVr(Mcut,Mcut)  )
    !--------------
    ! CALL mt%ED_tridiagonal_real(Mcut ,A(1:Mcut),B(1:Mcut),eigenVr,energy)
    CALL Lan_Diag3Matrix(Mcut ,A(1:Mcut),B(1:Mcut),eigenVr,energy)
    GS = (0._8,0._8)
    DO JC1=1 , Mcut
        GS = GS + SM(:,JC1) * eigenVr(JC1,1_8)
    ENDDO

    EG=energy(1)
    deallocate(  energy   )
    deallocate(  eigenVr  )
    LANCZOS_ITERATE = 0

    999 continue
  ENDfunction


  subroutine Lan_Diag3Matrix(N,A,B,H,E)
    use class_numerical_method
    implicit none
    integer,intent(in)::N
    real*8,intent(in)::A(N),B(N)
    real*8,intent(inout)::H(N,N),E(N)
    !-----------------------
    integer::jc,ifo
    type(nummethod)::mt

   if (n.eq.0)then
     !------------------
     ! do nothing
     !------------------
    !  write(*,*)"???"
    !  stop
     goto 999
   elseif (n.eq.1)then
      ! dimension = 1
      H(1,1) = 1._8
      E(1)   = A(1)
    else
      call mt%ED_tridiagonal_real(N ,A,B,H,E,ifo)
      if (ifo.ne.0)then
        H = 0._8
        H(1,1) = A(1)
        do jc = 2 , n
          H(jc,jc) =A(jc)
          H(jc,jc-1) = b(jc)
          H(jc-1,jc) = b(jc)
          call mt%ED_Hermitian_matrix("U",N,H,E)
        enddo
      endif
    endif
999 continue
  endsubroutine


    ! without checing anything. We scan space and find out a GS(othognal to the exsited GS)
    integer function ScanOneGroundState(self,GS,EG)
      implicit none
      class(LASubSpace),intent(inout)::self
      complex*16,intent(out)::GS(self%d)
      real*8,intent(out)::EG
      !--------------------------------------
      real*8::Etemp
      integer::con,ierr
      logical::Converged
      !--------------------------------
      call GetRadomState(self,self%d,GS)
      !--------------------------------
      con        = 0
      Converged  = .false.
      Etemp      = huge(Etemp)/2
      do while (.not.converged)
         ierr = LANCZOS_ITERATE(self,GS,EG)
         ScanOneGroundState = ierr
         con = con + 1
         !--------------------show --------------------------
         if (self%show .ge.2)then
            if (ierr.eq.0)then
             write(self%print,*)"Intera: ",int(con),". Eg = ",Eg
           elseif(ierr.eq.1)then
             write(self%print,*)"invalid input GS( othognal to existed GS)"
             write(self%print,*)"maybe M>d"
            endif
         endif
         !---------------------------------------------------
         if (ierr.eq.1)then
            goto 999
         endif
         if (abs(EG-Etemp).le.self%LanPre)then
           converged = .true.
         else
           Etemp = EG
         endif
      enddo
  999 continue
    endfunction


    !  <GS|GS> .ne. 1  !!!!
    subroutine GetRadomState(self,d,GS)
      use basic_math_functions
      class(LASubSpace),intent(inout)::self
      integer,intent(in)::d
      complex*16,intent(out)::GS(d)
      !----------------------------------------
      TYPE(bmathf)::f
      real*8::GSR(d)
      if (self%isreal)then
         call f%get_random_array(GSR)
         GS = GSR
      else
         call f%get_random_array(GS)
      endif                    !;write(*,*)gs;stop
    endsubroutine


    subroutine Diagonalization(self)
      implicit none
      class(LASubSpace),intent(inout)::self
      !---------------------------------------------------------
      class(State),pointer::S
      logical::Finished
      real*8::Eg
      integer::ierr
      !---------------------------------------------------------
      !   Initialization
          call self%state%Initialization(print_=self%print)
          Finished = .false.
      !---------------------------------------------------------
      !--------------------show --------------------------
      if (self%show .ge.2)then
         write(self%print,*)"Searching for the first ground state..."
      endif
      !---------------------------------------------------
      !  find the first one directly
         allocate( s           )
         allocate( s%s(self%d) )
         ierr = ScanOneGroundState(self,s%s,s%e)
         self%Eg = s%e
         if (ierr.ne.0)then
            write(self%print,*)"finding for the first GS is faild. Unknow reason";stop
         endif
         call self%state%append(s)
      !----------------------------------------------------------
      do while (.not.Finished)
        allocate( s           )
        allocate( s%s(self%d) )
        !--------------------show --------------------------
        if (self%show .ge.2)then
          write(self%print,*)" Searching for the",self%state%getlen()+1,"-th ground state..."
        endif
        !---------------------------------------------------
        ierr = ScanOneGroundState(self,s%s,s%e)
        if (ierr.eq.1)then
          ! maybe M > d
          Finished = .true.
          goto 200
        endif

        if (   CheckDiffSmall(self,s%e,self%EG,self%DegPre)  )then
          call self%state%append(s)
        else
          deallocate(s)
          Finished = .true.
        endif
200     continue
        !--------------------show --------------------------
        if (self%show .ge.2)then
          if (Finished)then
            write(self%print,*)"All GS have been found,Degeneracy=",self%state%getlen()
          else
            write(self%print,*)"New GS is found,Degeneracy=",self%state%getlen()
          endif
        endif
        !---------------------------------------------------
      enddo
      self%EigenId = self%H%GetEigenId()
    endsubroutine

    logical function CheckDiffSmall(self,E1,E2,rePre)
      use basic_math_functions
      implicit none
      class(LASubSpace),intent(inout)::self
      real*8,intent(in)::E1,E2,rePre
      !----------------------------------------
      real*8::r1,r2
      TYPE(bmathf)::f
      r1 = e1 ;if (abs(r1)<=self%LanPre) r1 = 0._8
      r2 = e2 ;if (abs(r2)<=self%LanPre) r2 = 0._8
      CheckDiffSmall = f%IsTwoValuePercentageTheSame(r1,r2,rePre)
    endfunction


    impure elemental subroutine SynchronizeWithHamiltonian(self)
      implicit none
      class(LASubSpace),intent(inout)::self
      !---------------------------------------------------------
      if (self%EigenId .ne. self%H%GetEigenId())then
         call Diagonalization(self)
      endif
    endsubroutine






  real*8 function GetEg(self)
    implicit none
    class(LASubSpace),intent(inout)::self
    !---------------------------------------------------------
    GetEg = self%Eg
  endfunction

  integer function GetDe(self)
    implicit none
    class(LASubSpace),intent(inout)::self
    !---------------------------------------------------------
    GetDe = self%state%getlen()
  endfunction

  real*8 function GetE(self,i)
    implicit none
    class(LASubSpace),intent(inout)::self
    integer,intent(in)::i
    !---------------------------------------------------------
    class(State),pointer::p
    p => self%state%GetDataPointer(i)
    GetE = p%E
  endfunction


  function GetOptActState(self,opt,i,d)  result(r)
    implicit none
    class(LASubSpace),intent(inout)::self
    class(FermOper),intent(inout)::opt
    integer,intent(in)::d,i
    complex*16::r(d)
    !---------------------------------------------------------
    complex*16::v
    integer::di
    class(State),pointer::p
    v = (1._8,0._8)
    p => self%State%GetDataPointer(i)
    di = size( p%s  )
                                           !!!!!write(*,*)d
    r = (0._8,0._8)
    call oper_Act_On_State(Ta=self%Ta,subidin=self%subid,oper=opt,V=v,&
                        IsReal=self%IsReal,din=di,dout=d,StateIn=p%s,AddStateOut=r)

    !write(*,*)DOT_PRODUCT(p%s,p%s)
    !write(*,*)r
  endfunction


  function GetTablePoinnter(self)  result(r)
    implicit none
    class(LASubSpace),intent(inout)::self
    class(table),pointer::r
    !---------------------------------------------------------
    r => self%Ta
  endfunction

  function GetHamPointer(self) result(r)
    implicit none
    class(LASubSpace),intent(inout)::self
    class(Ham),pointer::r
    !---------------------------------------------------------
    r => self%H
  endfunction


  logical function IsSysReal(self)
    implicit none
    class(LASubSpace),intent(inout)::self
    !---------------------------------------------------------
    IsSysReal = self%isreal
  endfunction

  INTEGER FUNCTION GetEigenId(SELF)
    implicit none
    class(LASubSpace),intent(inout)::self
    !---------------------------------------------------------
    GetEigenId = SELF%EIGENID
  ENDFUNCTION


  integer function GetSubId(self)
    implicit none
    class(LASubSpace),intent(inout)::self
    !---------------------------------------------------------
    GetSubId = self%subid
  endfunction


  complex*16 function GetGsProduct(self,i,S)
    implicit none
    class(LASubSpace),intent(inout)::self
    integer,intent(in)::i
    complex*16,intent(in)::S(self%d)
    !---------------------------------------------------------
    class(State),pointer::p
    p => self%State%GetDataPointer(i)
    GetGsProduct =  GetDotProduct( self%d , self%IsReal,p%s , S)
                                    !write(*,*)self%d,666;stop
  endfunction


  integer function GetNs(self)
    implicit none
    class(LASubSpace),intent(inout)::self
    !--------------------------------------------
    GetNs = self%Ta%get_ns()
  endfunction


  function GetOperActOnState(self,opert,subidin,din,dout,Si) result(r)
    implicit none
    class(LASubSpace),intent(inout)::self
    class(FermOper),intent(inout)::opert
    integer,intent(in)::subidin,din,dout
    complex*16,intent(in)::Si(din)
    complex*16::r(dout)
    !--------------------------------------------
    r = (0._8,0._8)
    call oper_Act_On_State(Ta=self%Ta,SubIdIn=subidin,oper=opert,V=(1._8,0._8)&
        ,IsReal=self%isreal,din=din,dout=dout,StateIn=si,AddStateOut=r)
  endfunction


  integer function GetD(self)
    implicit none
    class(LASubSpace),intent(inout)::self
    !-----------------------------------------
    GetD = self%d
  endfunction


  complex*16 function GetOperateProduct(self,opt)
    implicit none
    class(LASubSpace),intent(inout)::self
    class(FermOper),intent(inout)::opT
    !-----------------------------------------
    integer::jc
    complex*16::S(SELf%D)
    GetOperateProduct = (0._8,0._8)
    do jc = 1 , self%state%getlen()
      S = SELF%GetOptActState(opt,JC,SELF%D)
      GetOperateProduct = GetOperateProduct + SELF%GetGsProduct(JC,S)
    enddo
    GetOperateProduct = GetOperateProduct / self%state%getlen()
  endfunction





endmodule
