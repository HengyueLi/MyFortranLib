






!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : FermionOperators
! OBJECT: TYPE(FermOper)
! USED  :
! DATE  : 2017-12-06
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!              offer the operation on fermion spin 1/2 system. THe maxmuer size of the system may be L=18.
!            This size is determined by the max size subspace  L = ( C_N^(N/2) ) ^ 2
!              The fermion states are saved in integer*8 of the form (spinup,spindown).
!             for a system of size L = Ns , the index is from 0 to Ns-1
!             spin = 0 for up and 1 for down
!
!
! STANDARD:
!           [sub] Initialization (see "avalable sets" for details)
!
!            after the type is initiated, type(FermOper)::o for instance, use it as :
!            call o%act(BasisIn,BasisOut,sign_) (  see that in "avalable others"  )
!  ╔═════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╗
!  ║  optid = -1:                                                                                                                    ║
!  ║              always returen unavalabel state   ( output = -1)                                                                   ║
!  ║  optid = 0 :                                                                                                                    ║
!  ║              do nothing                                                                                                         ║
!  ║  optid = 1 : c_{i,spin}                                                                                                         ║
!  ║              i=para(1) ; spin=para(2)                                                                                           ║
!  ║                                                                                                                                 ║
!  ║  optid = 2 : c^+_{i,spin}                                                                                                       ║
!  ║              i=para(1) ; spin=para(2)                                                                                           ║
!  ║                                                                                                                                 ║
!  ║  optid = 3 : c^+_{i,sini} * c_{j,spinj}                                                                                         ║
!  ║              i=para(1) ; spini=para(2) ; j=para(3) ; spinj=para(4)                                                              ║
!  ║                                                                                                                                 ║
!  ║  optid = 4 : c_{i,sini} * c_{j,spinj}                                                                                           ║
!  ║              i=para(1) ; spini=para(2) ; j=para(3) ; spinj=para(4)                                                              ║
!  ║                                                                                                                                 ║
!  ║  optid = 5 : c^+_{i,sini} * c^+_{j,spinj}                                                                                       ║
!  ║              i=para(1) ; spini=para(2) ; j=para(3) ; spinj=para(4)                                                              ║
!  ║                                                                                                                                 ║
!  ║  optid = 6 : c^+_{i,sini} * c_{j,spinj} * c^+_{k,sini} * c_{l,spinj}                                                            ║
!  ║              i=para(1) ; spini=para(2) ; j=para(3) ; spinj=para(4) ; k=para(5) ; spink=para(6) ; l=para(7) ; spinl=para(8) ;    ║
!  ║                                                                                                                                 ║
!  ║  optid = 7 : c^+_{i,sini} * c^+_{j,spinj} * c_{k,sini} * c_{l,spinj}                                                            ║
!  ║              i=para(1) ; spini=para(2) ; j=para(3) ; spinj=para(4) ; k=para(5) ; spink=para(6) ; l=para(7) ; spinl=para(8) ;    ║
!  ║                                                                                                                                 ║
!  ║  optid = 8 : n_{i,spini} * n_{j,spinj}                                                                                          ║
!  ║              i=para(1) ; spini=para(2) ; j=para(3) ; spinj=para(4)                                                              ║
!  ║                                                                                                                                 ║
!  ║  optid = 9 : n_{i,spini}                                                                                                        ║
!  ║              i=para(1) ; spini=para(2)                                                                                          ║
!  ╚═════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝
!  usually for para(i), odd i=2n+1 represent site index and even i=2n represent spin
!
! USING LIST:
!
!
!
! avalable sets:
!                  [sub] Initialization(ns,optid,para,print_)
!                        integer::ns  !size of the system
!                        integer::optid  !chose operators
!                        integer::para(8)!information for operators
!                        integer::print_  ! where to print, default is 6
!                  [sub] UnInitialization()
!
! avalable gets:
!
!                  [fun] get_para()
!                        integer::get_para(8)
!
!                  [fun] get_optid()
!                        integer::get_optid
! avalable IS  :
!                  [fun] IsInitiated()
!
! avalable others:
!                  [sub] self%act(Bin,Bout,sign)
!                        integer(8)::Bin(2),Bout(2),sign
!                       Bin is the input basis, Bout is the output basis, sign= 1 / -1 is the correpsonding sign.
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

















module FermionOperators
     implicit none

    type FermOper
         private
         logical                          :: Initialized=.false.
         integer                          :: ns                       !total correlated sites in the system.
         integer                          :: para(8) !operator informations.
         integer                          :: optid
         integer                          :: print=6
         procedure(action),pointer,nopass :: opt => null()
         !--------------------------------------------------------
         ! internal use,to simplify para(8)
         ! para(8) will not do not directly inputted input subs, while Bsite(8) do.
         integer                          :: Bsite(4)
         !  the following two (opttemp1,opttemp2) are for bulting two particle operators.
         procedure(action),pointer,nopass :: opttemp1 => null() !
         procedure(action),pointer,nopass :: opttemp2 => null()
         integer                          :: Bsite2(4) ! cccc use
    contains
         !---------------------------------init,final
         procedure,pass::Initialization
         procedure,pass::UnInitialization
         final         ::Finalization

         procedure,pass::act
        !  procedure,pass::IsTheSameOperator
        !  procedure,pass::get_symmetry
         procedure,pass::get_para
         procedure,pass::get_optid
         procedure,pass::IsInitiated

         procedure::CheckTwoTheSame
         generic  :: operator(.eq.) => CheckTwoTheSame
    end type

    !--------------------------------------------------------------------------
    abstract interface
       subroutine  action(self,Bsite,basisin,basisout,sign_) ! para is saved in self
          import FermOper
          implicit none
          class(FermOper),intent(inout)::self
          integer   ,intent(in)::Bsite(4)
          integer*8,intent(in) ::basisin
          integer*8,intent(out)::basisout
          integer,intent(out)  ::sign_
       end subroutine action
    end interface
   !----------------------------------------------------------------------------





    private:: Initialization,UnInitialization,Finalization,act

    private:: select_operator

    !---------------------optid = -1
    private:: Forbi
    !---------------------optid = 0
    private:: DoNothing
    !---------------------optid = 1
    private:: SingleC
    !---------------------optid = 2
    private:: SingleCdagger
    !---------------------optid = 3
    private:: CPC_I_ne_J
    !---------------------optid = 4
    private:: CC_I_gt_J   , CC_I_lt_J
    !---------------------optid = 5
    private:: CpCp_I_gt_J , CpCp_I_lt_J
    !---------------------optid = 8
    private:: nn
    !---------------------optid = 9
    private:: particleN
    !---------------------optid = 6,7
    private:: cccc




    private::Check_Sign,get_optid

    private::checksite,checkspin

    private::IsInitiated,CheckTwoTheSame
  contains


    subroutine Initialization(self,ns,optid,para,print_)
      implicit none
      class(FermOper),intent(inout)::self
      integer,intent(in)::ns
      integer::optid
      integer,intent(in)::para(8)
      integer,optional::print_
      !----------------------------------------------------------
      call UnInitialization(self)
      !---------------------------
      self%Initialized = .true.
      self%ns          =  ns
      self%para        =  para
      self%optid       =  optid
      !------
      if (present(print_)) self%print = print_
      !------------------ select operator -----------
      call select_operator(self%ns,self%optid,self%para,self%Bsite,self%Bsite2&
                          ,self%opt,self%opttemp1,self%opttemp2,self%print)
      !----------------------------------------------
    endsubroutine

    subroutine UnInitialization(self)
      implicit none
      class(FermOper),intent(inout)::self
      !----------------------------------------------------------
      if (self%Initialized)then
        self%opt => null()
        self%Initialized = .false.
      endif
    endsubroutine

    impure elemental subroutine Finalization(self)
      implicit none
      type(FermOper),intent(inout)::self
      !----------------------------------------------------------
      call UnInitialization(self)
    endsubroutine


    recursive subroutine select_operator(ns,optid,para,Bsite,Bsite2,opt,opt1,opt2,print)
      implicit none
      integer,intent(in)::ns,optid,para(8)
      integer,intent(inout)::Bsite(4),Bsite2(4)
      procedure(action),pointer::opt,opt1,opt2
      integer,intent(in)::print
      !----------------------------------------------------------
      integer::optidcc1,optidcc2,para1(8),para2(8)
      !----------------  setting  bsite -------------------------
      Bsite(1) = para(2) * ns + para(1)
      Bsite(2) = para(4) * ns + para(3)
      Bsite(3) = min(Bsite(1) , Bsite(2)) + 1 !Bsite3: starting checking point
      Bsite(4) = abs(Bsite(1) - Bsite(2)) - 2 ! Bsite4: checking length
      !----------------------------------------------------------
      !----------------  selection  ------------------------------
      select case(optid)
      case(-1)!===========   optid = -1  ================================================
        opt => Forbi
      case(0)!===========    optid = 0   ================================================
        opt => DoNothing
      case(1)!===========    optid = 1   ================================================
        call checksite(ns,para(1))
        call checkspin(   para(2))
        opt => SingleC
      case(2)!===========    optid = 2   ================================================
        call checksite(ns,para(1))
        call checkspin(   para(2))
        opt => SingleCdagger
      case(3)!===========    optid = 3   ================================================
        call checksite(ns,para(1))  ; call checksite(ns,para(3))
        call checkspin(   para(2))  ; call checkspin(   para(4))
        if     (Bsite(1).ne.Bsite(2)) then
          opt => CPC_I_ne_J
        else
          opt => particleN
        endif
      case(4)!===========    optid = 4   ================================================
        if (Bsite(1).gt.Bsite(2))then
          opt => CC_I_gt_J
        elseif (Bsite(1).lt.Bsite(2))then
          opt => CC_I_lt_J
        else
          opt => Forbi
        endif
      case(5)!===========    optid = 5   ================================================
        if (Bsite(1).gt.Bsite(2))then
          opt => CpCp_I_gt_J
        elseif (Bsite(1).lt.Bsite(2))then
          opt => CpCp_I_lt_J
        else
          opt => Forbi
        endif
      case(6,7)!===========  optid = 6,7  ================================================
        if (optid.eq.6)then
          optidcc1 = 3  ;  optidcc2 = 3
        else
          optidcc1 = 5  ;  optidcc2 = 4
        endif
        para1(1:4) = para(1:4)  ;para2(1:4) = para(5:8)
        call select_operator(ns,optidcc1,para1,Bsite,Bsite2,opt1,opt,opt2,print)
        call select_operator(ns,optidcc2,para2,Bsite2,Bsite,opt2,opt,opt1,print)
        opt => cccc
      case(8)!===========    optid = 8   ================================================
        opt => nn
      case(9)!===========    optid = 9   ================================================
        opt => particleN
      case default
        write(print,*)"ERROR: Unknow optid input @ FermionOperators";stop
      endselect
    endsubroutine


    subroutine act(self,BasisIn,BasisOut,sign_)
      implicit none
      class(FermOper),intent(inout)::self
      integer*8,intent(in)::BasisIn
      integer*8,intent(out)::BasisOut
      integer,intent(out)::sign_
      !------------------------------------
      call self%opt(self,self%Bsite,BasisIn,BasisOut,sign_)
    endsubroutine



    !===========================================================================
    !  selected operators

    subroutine cccc(self,Bsite,basisin,basisout,sign_)
      implicit none
      class(FermOper),intent(inout)::self
      integer   ,intent(in)::Bsite(4)
      integer*8,intent(in) ::basisin
      integer*8,intent(out)::basisout
      integer,intent(out)  ::sign_
      !------------------------------------------
      integer::sign__
      integer*8::stemp
      call self%opttemp2(self,self%Bsite2,basisin,basisout,sign__)
      if (basisout.ne.-1_8)then
         stemp = basisout
         call self%opttemp1(self,self%Bsite,stemp,basisout,sign_)
         sign_ = sign_ * sign__
      endif
    endsubroutine cccc



    subroutine Forbi(self,Bsite,basisin,basisout,sign_)
      implicit none
      class(FermOper),intent(inout)::self
      integer   ,intent(in)::Bsite(4)
      integer*8,intent(in) ::basisin
      integer*8,intent(out)::basisout
      integer,intent(out)  ::sign_
      !------------------------------------------
      basisout = -1_8
      sign_    = 0
    endsubroutine Forbi



    subroutine DoNothing(self,Bsite,basisin,basisout,sign_)
      implicit none
      class(FermOper),intent(inout)::self
      integer   ,intent(in)::Bsite(4)
      integer*8,intent(in) ::basisin
      integer*8,intent(out)::basisout
      integer,intent(out)  ::sign_
      !------------------------------------------
      sign_    = 1_8
      basisout = basisin
    endsubroutine DoNothing

    ! the first position is started from 0
    ! check from pos to pos+len
    integer function Check_Sign(basis,pos,len)
        implicit none
        integer*8,intent(in)::basis
        integer,intent(in)::pos,len
        !-----------------------
        integer::jc
        Check_Sign = 1
        do jc = 0,len
           if (btest(basis,pos+jc)) Check_Sign=-Check_Sign
        enddo
    endfunction Check_Sign

    ! c_{i,spin}
    subroutine SingleC(self,Bsite,basisin,basisout,sign_)
      implicit none
      class(FermOper),intent(inout)::self
      integer   ,intent(in)::Bsite(4)
      integer*8,intent(in) ::basisin
      integer*8,intent(out)::basisout
      integer,intent(out)  ::sign_
      !------------------------------------------
      if (btest(basisin,Bsite(1)))then
         basisout = IBCLR(basisin,Bsite(1))
         sign_    = Check_Sign(basisin,0,Bsite(1)-1)
      else
         basisout = -1_8
      endif
    endsubroutine SingleC

    ! c^+_{i,spin}
    subroutine SingleCdagger(self,Bsite,basisin,basisout,sign_)
      implicit none
      class(FermOper),intent(inout)::self
      integer   ,intent(in)::Bsite(4)
      integer*8,intent(in) ::basisin
      integer*8,intent(out)::basisout
      integer,intent(out)  ::sign_
      !------------------------------------------
      if (.not.btest(basisin,Bsite(1)))then
         basisout = IBSET(basisin,Bsite(1))
         sign_    = Check_Sign(basisin,0,Bsite(1)-1)
      else
         basisout = -1_8
      endif
    endsubroutine SingleCdagger


    ! c^+_{i,spin}c_{j,spinj} where i .NE. j
    subroutine CPC_I_ne_J(self,Bsite,basisin,basisout,sign_)
      implicit none
      class(FermOper),intent(inout)::self
      integer   ,intent(in)::Bsite(4)
      integer*8,intent(in) ::basisin
      integer*8,intent(out)::basisout
      integer,intent(out)  ::sign_
      !------------------------------------------
      if ( btest(basisin,Bsite(2)) .and. (.not.btest(basisin,Bsite(1))) )then
         basisout = IBSET(basisin ,Bsite(1))
         basisout = IBCLR(basisout,Bsite(2))
         sign_    = Check_Sign(basisin,Bsite(3),Bsite(4))
      else
         basisout = -1_8
      endif
    endsubroutine CPC_I_ne_J

    ! c_{i,spin}c_{j,spinj} where i .lt. j
    subroutine CC_I_lt_J(self,Bsite,basisin,basisout,sign_)
      implicit none
      class(FermOper),intent(inout)::self
      integer   ,intent(in)::Bsite(4)
      integer*8,intent(in) ::basisin
      integer*8,intent(out)::basisout
      integer,intent(out)  ::sign_
      !------------------------------------------
      if ( btest(basisin,Bsite(2)) .and. btest(basisin,Bsite(1)) )then
         basisout = IBCLR(basisin ,Bsite(1))
         basisout = IBCLR(basisout,Bsite(2))
         sign_    = -Check_Sign(basisin,Bsite(3),Bsite(4))
      else
         basisout = -1_8
      endif
    endsubroutine CC_I_lt_J

    ! c_{i,spin}c_{j,spinj} where i .gt. j
    subroutine CC_I_gt_J(self,Bsite,basisin,basisout,sign_)
      implicit none
      class(FermOper),intent(inout)::self
      integer   ,intent(in)::Bsite(4)
      integer*8,intent(in) ::basisin
      integer*8,intent(out)::basisout
      integer,intent(out)  ::sign_
      !------------------------------------------
      if ( btest(basisin,Bsite(2)) .and. btest(basisin,Bsite(1)) )then
         basisout = IBCLR(basisin ,Bsite(1))
         basisout = IBCLR(basisout,Bsite(2))
         sign_    = Check_Sign(basisin,Bsite(3),Bsite(4))
      else
         basisout = -1_8
      endif
    endsubroutine CC_I_gt_J



    ! c^+_{i,spin}c^+_{j,spinj} where i .lt. j
    subroutine CpCp_I_lt_J(self,Bsite,basisin,basisout,sign_)
      implicit none
      class(FermOper),intent(inout)::self
      integer   ,intent(in)::Bsite(4)
      integer*8,intent(in) ::basisin
      integer*8,intent(out)::basisout
      integer,intent(out)  ::sign_
      !------------------------------------------
      if ( (.not.btest(basisin,Bsite(2))) .and. (.not.btest(basisin,Bsite(1))) )then
         basisout = IBSET(basisin ,Bsite(1))
         basisout = IBSET(basisout,Bsite(2))
         sign_    = Check_Sign(basisin,Bsite(3),Bsite(4))
      else
         basisout = -1_8
      endif
    endsubroutine CpCp_I_lt_J

    ! c^+_{i,spin}c^+_{j,spinj} where i .gt. j
    subroutine CpCp_I_gt_J(self,Bsite,basisin,basisout,sign_)
      implicit none
      class(FermOper),intent(inout)::self
      integer   ,intent(in)::Bsite(4)
      integer*8,intent(in) ::basisin
      integer*8,intent(out)::basisout
      integer,intent(out)  ::sign_
      !------------------------------------------
      if ( (.not.btest(basisin,Bsite(2))) .and. (.not.btest(basisin,Bsite(1))) )then
         basisout = IBSET(basisin ,Bsite(1))
         basisout = IBSET(basisout,Bsite(2))
         sign_    = -Check_Sign(basisin,Bsite(3),Bsite(4))
      else
         basisout = -1_8
      endif
    endsubroutine CpCp_I_gt_J



    ! N{i,spin}
    subroutine nn(self,Bsite,basisin,basisout,sign_)
      implicit none
      class(FermOper),intent(inout)::self
      integer   ,intent(in)::Bsite(4)
      integer*8,intent(in) ::basisin
      integer*8,intent(out)::basisout
      integer,intent(out)  ::sign_
      !------------------------------------------
      if ( BTEST(basisin,Bsite(1)) .and. BTEST(basisin,Bsite(2)) ) then
        basisout = basisin
        sign_    = 1_8
      else
        basisout = -1_8
      endif
    endsubroutine nn


    ! N{i,spin}
    subroutine particleN(self,Bsite,basisin,basisout,sign_)
      implicit none
      class(FermOper),intent(inout)::self
      integer   ,intent(in)::Bsite(4)
      integer*8,intent(in) ::basisin
      integer*8,intent(out)::basisout
      integer,intent(out)  ::sign_
      !------------------------------------------
      if ( BTEST(basisin,Bsite(1)) ) then
        basisout = basisin
        sign_    = 1_8
      else
        basisout = -1_8
      endif
    endsubroutine particleN











    integer function get_optid(self)
      implicit none
      class(FermOper),intent(inout)::self
      !----------------------------------------
      get_optid = self%optid
    endfunction

    function  get_para(self) result(r)
      implicit none
      class(FermOper),intent(inout)::self
      integer::r(8)
      !----------------------------------------
      r = self%para
    endfunction


    logical function IsInitiated(self)
      implicit none
      class(FermOper),intent(inout)::self
      !----------------------------------------
      IsInitiated = self%Initialized
    endfunction


    impure elemental subroutine checksite(ns,c)!,
      implicit none
      integer,intent(in)::ns
      integer,intent(in)::c
      !----------------------------------------
      if ((c.ge.0) .and.  (c.le.ns)  )then
      else
         write(*,*)"In definition of operator, site out of rang.";stop
      endif
    endsubroutine

    impure elemental subroutine checkspin(c)!,checkspin
      implicit none
      integer,intent(in)::c
      !----------------------------------------
      if ((c==0) .or.  (c==1)  )then
      else
         write(*,*)"In definition of operator, spin out of rang.";stop
      endif
    endsubroutine


    logical function CheckTwoTheSame(a,b)
      implicit none
      class(FermOper),intent(in)::a,b
      !-------------------------------------
      CheckTwoTheSame = .false.
      if (a%optid.eq.b%optid)then
         if (  sum(abs(a%para - b%para)).eq.0   )then
            CheckTwoTheSame = .true.
         endif
      endif
    endfunction



end module
