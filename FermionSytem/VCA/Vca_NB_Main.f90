
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : VCA_NB_MAIN
! OBJECT: TYPE(VCANB)
! USED  : CodeObject,fermion_table,LaPrimaryH,LatticeConfig,LaLatticeH,VCA_DeltaH,VCA_WaldFun,CEsolver,VCA_Variation
! DATE  : 2017-12-30
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            Vational Cluster Approach without bath sites.
!
!                       ╭────╮             ╭────╮
!                       │PriH│             │LaCo│
!                       ╰──┬─╯             ╰┬─┬─╯
!                          │                │ │
!                          │┌───────────────┘ │
!                          ││                 │
!   ╭────╮  ╭────╮      ╭──┴┴╮             ╭──┴─╮
!   │ GP │  │ EP │      │CPTH│             │DelH│
!   ╰─┬──╯  ╰─┬┬─╯      ╰──┬┬╯             ╰┬┬──╯
!     │       │└──────────┐│└───┐           ││
!     │       │          ╭┴┴──╮ │           ││
!     │       │          │EDTA│ │           ││
!     │       │          ╰──┬─╯ │           ││
!     │       │             └──┐│           ││
!     │       └───────────────┐││┌──────────┘│
!     └──────────────────────┐││││           │
!                           ╭┴┴┴┴┤           │
!                           │Wafu│           │
!                           ╰──┬─╯           │
!                              │             │
!                              └──────┬──────┘
!                                     │
!                                   ╭─┴──╮
!                                   │Vari│
!                                   ╰────╯
! STANDARD:
!            *CALL Initialization( )
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                   [sub] Initialization( print_,show_ )
!
! avalable gets:
!                   [fun] GetOverLapCPTH()
!                         reutrn a type(LH)::GetOverLapCPTH in which both H and dH have been considered.
!                         used for CPT calculation.
!
! avalable is :
!                  ![fun] i
! others      :
!                  ![sub] TestVariational(filepath,dH_Disc,mode,rmode,R,N)
!                          character(*),intent(in)::filepath
!                          character(32),intent(in):: dH_Disc
!                          integer,intent(in)::mode
!                          real*8,intent(in) ::rmode
!                          real*8,intent(in) ::R
!                          integer,intent(in)::N
!
!                          test a stationary point (on one direction.).
!                          the testing range is in [o-R,o+R] where o is the recent value.
!
!                          mode = 0 : Homoginouse distribution
!                          mode = 1 :
!
!
!                   [sub] TestAllVariational(FolderPath,mode,rmode,R,N)
!                         class(VCANB),intent(inout)::self
!                         character(*),intent(in)::FolderPath
!                         integer,intent(in)::mode
!                         real*8,intent(in) ::rmode
!                         real*8,intent(in) ::R
!                         integer,intent(in)::N
!
!                   [fun] CMsearching(ResetDH)
!                         logical::ResetDH
!                         integer::CMsearching    return the error code of searching process.
!                         if (ResetDH) : before CM, reset all parameters in dH to 0.
!
!                   [sub] ReportParameters()
!
!
!!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! override part
! ! procedure,pass::SetEDPara
! ! procedure,pass::SetPrimaryCellHamiltonian
! ! procedure,pass::SetLatticeConfiguration
! ! procedure,pass::SetDeltaMatrix
! ! procedure,pass::SetWaldPara
! ! procedure,pass::SetVariationalPara
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module VCA_NB_MAIN
  use CodeObject
  use fermion_table
  use LaPrimaryH
  use LatticeConfig
  use LaLatticeH
  use VCA_DeltaH
  use VCA_WaldFun
  use CEsolver
  use VCA_Variation
  IMPLICIT NONE


  type,extends(object)::VCANB
    ! private
    !-----------------------------
    ! solver parameters : this should be fixed for all calculation
      type(SolverPara) :: EDPA
    !-----------------------------
    ! primary Hamiltonian
      type(pH)         :: PriH
    !-----------------------------
    ! Lattice configuration
      TYPE(LaCon)      :: LaCo
    !-----------------------------
    ! CPTH
      type(LH)         :: CPTH
    !----------------------------
    ! meanfiled
      type(VCAdH)      :: DelH
    !---------------------------
    ! using for ED
      type(table)      :: EDTa
    !---------------------------
    ! waldfunctional
      type(waldf)      :: Wafu
    !---------------------------
    ! variational method
      type(VCAva)      :: vari

  contains
    procedure,pass::Initialization

    procedure,pass::TestVariational
    procedure,pass::TestAllVariational
    procedure,pass::CMsearching
    procedure,pass::ReportParameters
    procedure,pass::GetOverLapCPTH


   !!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
   ! override part
    procedure,pass::SetEDPara
    procedure,pass::SetPrimaryCellHamiltonian
    procedure,pass::SetLatticeConfiguration
    procedure,pass::SetDeltaMatrix
    procedure,pass::SetWaldPara
    procedure,pass::SetVariationalPara
  endtype

  private::Initialization


  private::TestVariational,TestAllVariational
  private::CMsearching,ReportParameters
  private::GetOverLapCPTH



  !!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! override part
  private::SetEDPara
  private::SetPrimaryCellHamiltonian
  private::SetLatticeConfiguration
  private::SetDeltaMatrix
  private::SetWaldPara
  private::SetVariationalPara
contains

  subroutine Initialization(self,print_,show_)
    implicit none
    class(VCANB),intent(inout)::self
    integer,intent(in),optional::print_,show_
    !---------------------------------
    type(GreenPara) :: GP
    integer         :: Wdjobi(20)
    real*8          :: Wdjobr(20)
    integer         :: Vajobi(20)
    real*8          :: Vajobr(20)

    call self%SetInitiated(.true.)
    if (present(print_)) call self%setprint(print_)
    if (present(show_ )) call self%setshow(show_  )


    !--------------------------------------------------------------------------------------------
    ! set solver parameters
    Call self%SetEDPara(self%edpa)
    !--------------------------------------------------------------------------------------------
    !  initiate PriH
    Call self%PriH%Initialization()
    call self%SetPrimaryCellHamiltonian(self%PriH)
    !--------------------------------------------------------------------------------------------
    ! inititate LaCo
    call self%SetLatticeConfiguration(self%LaCo)
    !--------------------------------------------------------------------------------------------
    ! initiate CPTH
    call self%CPTH%Initialization(self%PriH,self%LaCo,self%getprint(),self%getshow()+1 )
    !--------------------------------------------------------------------------------------------
    ! initiate DelH
    Call self%DelH%Initialization( self%LaCo,self%getprint(),self%getshow()+1)
    Call self%SetDeltaMatrix(self%DelH)
    !--------------------------------------------------------------------------------------------
    ! initiate EDTa
    Call self%EDTa%Initialization( self%LaCo%GetNs() ,   self%EDPA%symmetry   , self%getprint()   )
    !--------------------------------------------------------------------------------------------
    ! initiate Wafu
    Call self%SetWaldPara(GP,wdjobi,wdjobr)
    Call self%Wafu%Initialization(Ta=self%EDTa,SP=self%edpa,GP=gp,CPTH=self%cpth,&
                dH=self%DelH,jobi=Wdjobi,jobr=Wdjobr,print_=self%getprint(),show_=self%getshow()+1)
    !--------------------------------------------------------------------------------------------
    ! initiate vari
    Call self%SetVariationalPara(Vajobi,Vajobr)
    call self%vari%Initialization(dh=self%DelH,wf=self%Wafu,jobi=Vajobi,jobr=Vajobr,&
                                  print_=self%getprint(),show_=self%getshow()+1)


  endsubroutine


  ! test a stationary point (on one direction.).
  ! the testing range is in [o-R,o+R] where o is the recent value.
  !
  ! mode = 0 : Homoginouse distribution
  ! mode = 1 :
  subroutine TestVariational(self,filepath,dH_Disc,mode,rmode,R,N)
    implicit none
    class(VCANB),intent(inout)::self
    character(*),intent(in)::filepath
    character(32),intent(in):: dH_Disc
    integer,intent(in)::mode
    real*8,intent(in) ::rmode
    real*8,intent(in) ::R
    integer,intent(in)::N
    !-----------------------------------------------------
    real*8::mid,r1,r2,v
    integer::jc,print
    print = self%getprint()

    mid = self%DelH%GetValueByDiscription(dH_Disc)
    r1 = mid - abs(r)
    r2 = mid + abs(r)

    open(99,file=trim(adjustl(filepath)))
    do jc = 1 , n
       v = getr(mode,rmode,r1,mid,r2,jc,n,print)
       call self%DelH%SetValueByDiscription(dH_Disc,v)
       write(99,*)v,self%Wafu%GetLatticeOmegaPerSite()
    enddo
    close(99)
    !-----------------------
    ! recover ogiginal value
    call self%DelH%SetValueByDiscription(dH_Disc,mid)
  contains
    real*8 function getr(mode,rmode,r1,mid,r2,jc,n,print)
      implicit none
      integer,intent(in)::mode,jc,n,print
      real*8,intent(in) ::rmode,r1,r2,mid
      !-----------------------------------
      real*8::l
      real*8::s
      getr = (r2-r1)/n * (jc-0.5_8) + r1
      select case(mode)
      case(0)
      case(1)
        l = abs(getr - mid)
        s = sign(1._8,getr - mid)
        getr = ( l / (mid - r1) ) ** rmode
        getr = getr * l * s + mid
      case DEFAULT
        write(print,*)"ERROR: Unknow mode =",mode," in TestVariational@VCA";stop
      endselect
    endfunction
  endsubroutine

  ! see the illstration of TestVariational
  ! file will be named by discription and put on FolderPath
  subroutine TestAllVariational(self,FolderPath,mode,rmode,R,N)
    implicit none
    class(VCANB),intent(inout)::self
    character(*),intent(in)::FolderPath
    integer,intent(in)::mode
    real*8,intent(in) ::rmode
    real*8,intent(in) ::R
    integer,intent(in)::N
    !-----------------------------------------------------
    character(64)::filename
    character(32)::disc
    integer::jc
    do jc = 1 , self%DelH%GetNumOfVarialtinalTerms()
       disc = self%DelH%GetDisciption(jc)
       write(self%getprint(),"(A38,A32)")">>>Test Stationary Point at Direction:",disc
       filename = trim(adjustl(FolderPath))//"/"//trim(adjustl(disc))
       call TestVariational(self,filename,disc,mode,rmode,R,N)
       write(self%getprint(),"(A9)")"--->done."
    enddo
  endsubroutine


  integer function CMsearching(self,ResetDH)
    implicit none
    class(VCANB),intent(inout)::self
    logical,intent(in)::ResetDH
    !-----------------------------------
    CMsearching = self%vari%CrossOverSearching(ResetDH)
  endfunction


  subroutine ReportParameters(self)
    implicit none
    class(VCANB),intent(inout)::self
    !-----------------------------------
    Call self%CPTH%report(self%getprint())
    Call self%DelH%report(self%getprint())
  endsubroutine



  function GetOverLapCPTH(self) result(r)
    implicit none
    class(VCANB),intent(inout)::self
    type(LH)::r
    !-----------------------------------
    r = self%CPTH
    call r%AbsorbMeanFiled(self%DelH%GetIdataArray())
  endfunction


















  !!############################################################################
  !! override part
  !!############################################################################

  subroutine SetEDPara(self,edpa)
    implicit none
    class(VCANB),intent(inout)::self
    class(SolverPara),intent(inout)::edpa
    !--------------------------------------------
    write(self%getprint(),*)"ERROR: SetEDPara is not defined yet."
    write(self%getprint(),*)"Read document to set it. There is no procedure needed to be called."
    stop
  endsubroutine


  subroutine SetPrimaryCellHamiltonian(self,priH)
    implicit none
    class(VCANB),intent(inout)::self
    class(PH),intent(inout)::priH
    !--------------------------------------------
    write(self%getprint(),*)"ERROR: SetPrimaryCellHamiltonian is not defined yet."
    write(self%getprint(),*)"for a input priH(already intiated.) Read document to set it."
    stop
  endsubroutine

  subroutine SetLatticeConfiguration(self,laco)
    implicit none
    class(VCANB),intent(inout)::self
    class(LaCon),intent(inout)::LaCo
    !--------------------------------------------
    write(self%getprint(),*)"ERROR: SetLatticeConfiguration is not defined yet."
    write(self%getprint(),*)"for a input laco(not intiated yet). Read document to initiate and set it."
    stop
  endsubroutine

  subroutine SetDeltaMatrix(self,delh)
    implicit none
    class(VCANB),intent(inout)::self
    class(VCAdH),intent(inout)::delh
    !--------------------------------------------
    write(self%getprint(),*)"ERROR: SetDeltaMatrix is not defined yet."
    write(self%getprint(),*)"for a input laco(already intiated). Read document to set it."
    stop
  endsubroutine

  subroutine SetWaldPara(self,gp,jobi,jobr)
    implicit none
    class(VCANB),intent(inout)::self
    class(GreenPara),intent(inout)::gp
    integer,intent(inout)::jobi(:)
    real*8,intent(inout) ::jobr(:)
    !-----------------------------------------
    write(self%getprint(),*)"ERROR: SetWaldPara is not defined yet."
    write(self%getprint(),*)" Read document to set it."
    stop
  endsubroutine

  subroutine SetVariationalPara(self,Vajobi,Vajobr)
    implicit none
    class(VCANB),intent(inout)::self
    integer,intent(inout)::vajobi(:)
    real*8,intent(inout) ::vajobr(:)
    !-----------------------------------------
    write(self%getprint(),*)"ERROR: SetVariationalPara is not defined yet."
    write(self%getprint(),*)" Read document to set it."
    stop
  endsubroutine


endmodule



















!
! module VCA_NB_MAIN
!   use CodeObject
!   use CEsolver    , only: VCASP => SolverPara
!   use CE_Green    , only: VCAGP => GreenPara
!   use LaPrimaryH
!   use LatticeConfig
!   use VCA_DeltaH
!   use LaLatticeH
!   use fermion_table
!   use VCA_Types
!   use VCA_WaldFun
!   implicit none
!
!
!   type,extends(object)::VCANB
!     ! private
!
!     !--------------------------------------------------------------------
!     ! saving a Hamitonian in a primary cell. No geometry is contained.
!       TYPE(PH),pointer    :: PrimaryH
!     !--------------------------------------------------------------------
!     ! the lattice configuration of the system.
!       TYPE(LaCon)         :: LattConf
!     !--------------------------------------------------------------------
!     ! variational meanfiled
!       TYPE(VCAdH)         :: DeltaH
!     !--------------------------------------------------------------------
!     ! wald functional
!       TYPE(waldf)         :: WaFun
!
!
!
!
!     !====================================================================
!     !  inner usage
!     !====================================================================
!     ! Frequently be used by solver. Thus we save the value.
!       Type(table) :: SolverTa
!     ! Cpt H
!       TYPE(LH) :: CPTH
!     !--------
!       type(VCASP)::sp
!
!
!   contains
!     !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     !  override
!
!     !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!
!
!
!     procedure,pass::Initialization
!     final::Finalization
!
!   endtype
!
!
!   private::Initialization,Finalization
!
!
!
!
! contains
!
!
!
!
!
!   subroutine Initialization(self,WaldPara,LattPara,PriH,DHarray ,print_,show_)
!     implicit none
!     class(VCANB),intent(inout)        :: self
!     class(VCA_Wald_Para),intent(in) :: WaldPara
!     class(VCA_Latt_Para),intent(in) :: LattPara
!     class(PH),intent(inout),target  :: PriH
!     class(VCA_DelH_Para),intent(in) :: DHarray(:)
!     integer,intent(in),optional     :: print_,show_
!   !--------------------------------------------------------------------
!     type(VCASP)::sp
!     type(VCAGP)::gp
!     integer::jc
!     !------------------------------------
!     call Finalization(self)
!     call self%SetInitiated(.true.)
!     if (present(print_)) call self%setprint(print_)
!     if (present(show_ )) call self%setshow(show_  )
!     !------------------------------------
!
!     sp%SvType       =  WaldPara%SolverType
!     sp%IsReal       =  WaldPara%IsReal
!     sp%Temperature  =  WaldPara%Temperature
!     sp%pre          =  WaldPara%pre
!     sp%bzero        =  WaldPara%bzero
!     sp%m            =  WaldPara%M
!     sp%oth          =  WaldPara%oth
!     sp%DegPre       =  WaldPara%DegPre
!     self%sp = sp
!
!
!     gp%M            =  WaldPara%GM
!     gp%oth          =  WaldPara%Goth
!     gp%bzero        =  WaldPara%Gbzero
!
!     allocate(self%PrimaryH ,source= PriH)
!     if (.not.self%PrimaryH%IsInitiated())then
!        write(self%getprint(),*)"ERROR: input PrimaryH is not initated in VCA";stop
!     endif
!
!     !-------------------------------------------------------------------------------
!     ! set Lattice configuration
!     call self%LattConf%Initialization(  Vp    = LattPara%PCBasis,&
!                                         PC    = LattPara%SitePos(:,1:LattPara%NSiteInPC),&
!                                         PCi   = LattPara%OrbIndx(1:LattPara%NSiteInPC)  ,&
!                                         Vl    = LattPara%LCBasis   )
!     !-------------------------------------------------------------------------------
!     ! set CPT H
!     call self%CPTH%Initialization( PrH = self%PrimaryH,LaC = self%LattConf,&
!                                    print_= self%getprint(),show_=self%getshow()+1 )
!     !--------------------------------------------------------------------------------
!     ! set table
!     call self%SolverTa%Initialization( self%LattConf%GetNs()  ,   WaldPara%symmetry  )
!     !--------------------------------------------------------------------------------
!     ! set dH
!     call self%DeltaH%Initialization( self%LattConf  )
!     call self%DeltaH%StartAppending()
!     !**************************************
!     do jc = 1 , size(DHarray)
!       call self%DeltaH%AppendVariationalTerm(DHarray(jc)%Discription,DHarray(jc)%Type  ,&
!                                    DHarray(jc)%jobi, DHarray(jc)%jobr  )
!     enddo
!     !**************************************
!     call self%DeltaH%EndAppending()
!     !--------------------------------------------------------------------------------
!     call self%WaFun%Initialization(Ta=self%SolverTa,SP=sp,GP=gp,CPTH=self%CPTh,dH=self%DeltaH,&
!                                  jobi=WaldPara%jobi,jobr=WaldPara%jobr,&
!                                  print_= self%getprint(),show_=self%getshow()+1 )
!
!   endsubroutine
!
!
!   subroutine Finalization(self)
!     implicit none
!     type(VCANB),intent(inout)  :: self
!     if (self%IsInitiated())then
!        call self%SetInitiated(.false.)
!        deallocate(self%PrimaryH)
!     endif
!   endsubroutine
!
!
!
!
!
!
!
!
!
!
!
! endmodule
!














!
! subroutine SetPrimaryCellHamiltonian(self,PrimaryH)
!   implicit none
!   class(VCANB),intent(inout)::self
!   class(PH),intent(inout)   ::PrimaryH
!   !------------------------------------
!   write(self%getprint(),*)"ERROR: SetPrimaryCellHamiltonian is not defined yet."
!   write(self%getprint(),*)"-------------------------------------------------------------------"
!   write(self%getprint(),*)"use Subroutine PrimaryH%Append(Disc,Pos,Itype,Ipara,Ivalu) where"
!   write(self%getprint(),*)"character(DiscLen)  :: Disc          ( DiscLen = 32?  )"
!   write(self%getprint(),*)"integer             :: pos(3)         denote cluster position"
!   write(self%getprint(),*)"character(ITypeLen) :: Itype         ( ITypeLen = 16? )"
!   write(self%getprint(),*)"integer             :: Ipara(8)       1-basis for site index.  0,1 for spin up and down."
!   write(self%getprint(),*)"complex*16          :: Ivalu"
!   write(self%getprint(),*)"Itype can be checked in module FermionHamiltonian."
!   write(self%getprint(),*)"And only for the terms that the first two elements in para(8) represent sites can be used."
!   stop
! endsubroutine
!
