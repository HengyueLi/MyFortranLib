


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : LaLatticeH
! OBJECT: TYPE(LH)
! USED  : CodeObject,LaPrimaryH,LatticeConfig,FermionHamiltonian,functionalsubs
! DATE  : 2017-12-27
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            Fully discribe a Hamitonian of a Lattice system by choosing a lattice cell.
!
! STANDARD:
!            *CALL Initialization( PrH,LaC,print_,show_ )
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                   [sub] Initialization( PrH,LaC,print_,show_ )
!                          class(PH),target,intent(in)   :: PrH
!                          class(LaCon),target,intent(in):: Lac
!                          integer,intent(in),optional   :: print_,show_
!
!                          PrH and  Lac should be initiated previously.
!
!                   [sub] SetValueByDiscription(Dis,v)
!                         character(DiscLen)::Dis     ( DiscLen = 32?  )
!                         complex*16::V
!
!                         scan all interacting and reset value of whose discription is Dis
!
!                   [sub] SetMeanFieldState(tf)         !spin are 0-basis, sites are 1-basis.
!                         set meanfil or not.
!                         logical::tf
!
!                   [sub] SetMeanField(M,spini,spinj)   !spin are 0-basis, sites are 1-basis.
!                         complex*16::M(self%ns,self%ns)
!                         integer::spini,spinj
!
!
!
! avalable gets:
!                   [fun] GetNs()
!                         get total sites in LC
!
!                   [fun] GetVp()
!                         return PC basis
!
!                   [fun] GetVbasis()
!                         real*8::GetVbasis(3,3)  the basis of LC.
!
!                   [fun] GetNnearLC()
!                         total number of connected LC.
!
!                   [fun] GetLocalHMatix(spini,spinj)
!                         return complex*16::M(ns,ns)
!                         where ns is the LC size
!
!                   [fun] GetSpinSuppresedLocalHMatrix()
!                         2ns 2ns matrix
!
!                   [fun] GetNearLCpos(i)
!                         integer::i
!                         return real*8::GetNearLCpos(3)  the position of near LC in real space.
!
!                   [fun] GetNearTMatix(i,spini,spinj)
!                         return complex*16::M(ns,ns)
!                         the i-th connected T   (exp part is not included)
!
!                   [fun] GetNearTexpMatrix(q,i,spini,spinj)
!                         GetNearTMatix multiply by exp( i r.q )
!
!                   [fun] GetTqMatrix(q,spini,spinj)
!                         summation of GetNearTexpMatrix,  this is Tq
!
!                   [fun] GetSpinSuppresedTq(q)
!                         2ns X 2ns matrix where both spin have been contained.
!
!                   [fun] GetLatticeConfigPointer()
!                         return class(LaCon)   lattice configration
!
!                   [fun] GetMeanFiled(spini,spinj)  !spin are 0-basis, sites are 1-basis.
!                         complex*16::GetMeanFiled(self%ns,self%ns)
!
!                   [fun] GetMeanFieldState()
!                         logical
!
!                   [fun] GetVelocityOneCluster(q,spini,spinj,xy,i)  result(r)! i = [0,self%GetNnearLC()]
!                          class(LH),intent(inout) :: self
!                          real*8,intent(in)       :: q(3)
!                          integer,intent(in)      :: spini,spinj,xy,i
!                          complex*16              :: r(self%Ns,self%Ns)
!
!                          r = -   ( T_R * exp(i q R) * (  R - X   )     )
!                          R is the cluster position, X is ri-rj in the cluster.
!
!                   [fun] GetVelocityMatrix(q,spini,spinj,xy)
!                         see  'GetVelocityOneCluster'. All the near connected clusters are sumed.
!
!
!
!
! avalable is :
!                  ![fun] i
! others      :
!                   [sub] AppendLocalDataToHam(H)
!                        for a input Type(Ham)::H, append all the local interacting (type idata) into it.
!                        H maybe contains some other terms already, but here it does not check that.
!
!                   [sub] Report(wtp)
!                         integer::wtp
!                         show the recently paramters
!
!                   [sub] AbsorbMeanFiled(idataarry)
!                         for a idata array class(idata)::idataarry(:)
!                         put it into the local part.
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module LaLatticeH
  use CodeObject
  use LaPrimaryH
  use LatticeConfig
  use LaPrimaryHUsedList    , only : Ilist => ListStru
  use CPTInterType
  use FermionHamiltonian    , only : Ham
  implicit none



  integer,parameter::NSurroundedCLuster = 30  ! max number of surrounded Lattice Cell.


  type,extends(Object)::LH
    private

    integer :: Ns ! number of total sites in a LC

    class(PH),pointer      :: PrH  => null()    ! primary H
    class(LaCon),pointer   :: LaC  => null()    ! lattice Configration
    !--------------------------------
    !  the size is approximated as  Ni * Np  where Ni is number of interacing in PC and Np is number of PC
    integer                :: NHin
    type(idata),allocatable:: Hin(:)         ! list for saving all Hamiltonian in Lattice.Inner use
    !--------------------------------
    ! approximated as Ni * Np
    integer                :: NLC              ! number of surrounded LCs.
    ! save all the position of surrounded Lattice Cell.
    integer :: LaPos( 3,NSurroundedCLuster)   !In the basis of PC's
    real*8  :: LaposR(3,NSurroundedCLuster)   !In real space
    integer                :: NHou(NSurroundedCLuster)
    type(idata),allocatable:: Hou(:,:)         ! list for saving all Hamiltonian out of Lattice.Inner use

    ! class(Ham),pointer   :: Hm   => NULL()    ! Ham used for cluster solver.

    !------------------------------------------------
    ! the possible meanfiled .
    ! For pure CPT, there is not such kind of term.
    ! for other, VCA and DMFT, we have this term.
    logical::HaveMeanFiled
    complex*16,allocatable::DeltaM(:,:,:,:) !   (ns,ns,0:1,0:1) the later two are spin index.
    !-------------------------------------------------
    !  used for velocity matrix
    real*8,allocatable :: DeltaX(:,:,:)     ! (ns,ns,3)
    !
  CONTAINS
    procedure,pass::Initialization
    final::Finalization

    procedure,pass::GetNs
    procedure,pass::GetNearTMatix
    procedure,pass::GetNearTexpMatrix
    procedure,pass::GetNearLCpos
    procedure,pass::GetNnearLC
    procedure,pass::SetValueByDiscription
    procedure,pass::GetLocalHMatix
    procedure,pass::GetTqMatrix
    procedure,pass::AppendLocalDataToHam
    procedure,pass::GetVbasis
    procedure,pass::GetSpinSuppresedTq
    procedure,pass::GetSpinSuppresedLocalHMatrix
    procedure,pass::Report
    procedure,pass::AbsorbMeanFiled
    procedure,pass::GetVp
    procedure,pass::GetLatticeConfigPointer
    procedure,pass::SetMeanFieldState
    procedure,pass::SetMeanField
    procedure,pass::GetMeanField
    procedure,pass::GetMeanFieldState
    procedure,pass::GetVelocityOneCluster
    procedure,pass::GetVelocityMatrix
  endtype


  private::Initialization,Finalization

  private::SettingLatticeHbyPrimaryH
  private::AppendOneInterFromPrimaryCellToLatticeCell
  private::GetSurroundLC_IdFromPos,ExtendNegativePartH
  private::SetValueByDiscription
  ! private::SetMatixOneTerm

  private::GetNs
  private::GetNnearLC
  private::GetLocalHMatix
  private::GetNearTMatix
  private::GetNearTexpMatrix
  private::GetTqMatrix
  private::AppendLocalDataToHam
  private::GetVbasis
  private::GetSpinSuppresedTq
  private::GetSpinSuppresedLocalHMatrix
  private::GetNearLCpos
  private::Report
  private::AbsorbMeanFiled
  private::GetVp
  private::GetLatticeConfigPointer
  private::SetMeanFieldState,SetMeanField,getMeanField,GetMeanFieldState
  private::GetVelocityOneCluster
  private::GetVelocityMatrix

contains




  subroutine Initialization(self,PrH,LaC,print_,show_)
    implicit none
    class(LH),intent(inout)       :: self
    class(PH),target,intent(in)   :: PrH
    class(LaCon),target,intent(in):: Lac
    integer,intent(in),optional   :: print_,show_
    !----------------------------------------------
    integer::jc1,jc2
    call Finalization(self)
    call self%SetInitiated(.true.)
    if (present(print_)) call self%SetPrint(print_ )
    if (present(show_))  call self%SetShow( show_  )

    self%PrH => PrH
    self%Lac => Lac
    !-----------------------------------------
    !  check initiated : PrH, Lac
    if (    self%prh%IsInitiated()    .and.   self%Lac%IsInitiated()    ) then
    else
      write(self%GetPrint(),*)"ERROR: In Initialization of LH, PH and LaCon should be initiated first."
      stop
    endif
    !-----------------------------------------
    self%Ns = self%LaC%GetNs()
    !-----------------------------------------
    allocate( self%Hin( self%prH%GetLen() * self%LaC%GetNp()                      )   )
    self%NHin = 0
    allocate( self%Hou( self%prH%GetLen() * self%LaC%GetNp() ,NSurroundedCLuster  )   )
    self%NHou = 0
    self%NLC  = 0
    !-----------------------------------------
    call SettingLatticeHbyPrimaryH(self)
    call ExtendNegativePartH(self)

    !----------------------------
    self%HaveMeanFiled = .false.
    allocate(self%DeltaM(self%Ns,self%Ns,0:1,0:1))
    self%DeltaM = (0._8,0._8)
    !---------------------------------------------------
    ! velocity matrix
    allocate(self%DeltaX(self%ns,self%ns,3))
    do jc1= 1 , self%ns
      do jc2 = 1 , self%ns
         self%DeltaX(jc1,jc2,:) = self%LaC%GetSiteRealP(jc1) - self%LaC%GetSiteRealP(jc2)
      enddo
    enddo
    !
  endsubroutine

  subroutine Finalization(self)
    implicit none
    type(LH),intent(inout)::self
    !-----------------------------
    if (self%IsInitiated())then
       deallocate(   self%Hin  , self%Hou   , self%DeltaM , self%DeltaX )
       call self%SetInitiated( .false. )
    endif
  endsubroutine



  subroutine SettingLatticeHbyPrimaryH(self)
    implicit none
    class(LH),intent(inout)         :: self
    !--------------------------------------
    integer::jcp,jci
    type(idata)::Inter



    if (self%PrH%getstate()==3)then
    else
      write(self%getprint(),*)"ERROR: Primary Cell H is not set, only after which Lattice Cell H can be used"
      stop
    endif

    !--------------------------------
    ! for all primary cell
    do jcp = 1 , self%LaC%GetNp()
      !------------------------------
      ! For all interacting in one PC
      do jcI = 1 , self%PrH%getlen()
        !----------------------------
        ! For one interacting
        Inter = self%PrH%GetIdata(jci)                                      !;write(*,*)jcI;call inter%report(6)
        Call AppendOneInterFromPrimaryCellToLatticeCell(self,Inter,jcp)
      enddo                                                                  !;stop
    enddo
    !--------------------------------
  endsubroutine


    subroutine AppendOneInterFromPrimaryCellToLatticeCell(self,Inter,PCindex)
      implicit none
      class(LH),intent(inout)         :: self
      type(idata),intent(in)        :: Inter
      integer,intent(in)            :: PCindex
      !--------------------------------------
      type(idata)::Iin
      integer::i,j,PCid,LCid,LCpos(3)
      integer::Rj(3),xj(3),bj(3)    !  R = A . x + b

      Iin         = Inter

      !-----------------------------------------------------------------------
      ! Get fist index, the first index is always in the cluster.
      Iin%para(1) = self%LaC%GetSiteIdFromPC(PCindex,  inter%para(1)  )
      if (  (Iin%para(1) .lt.1)  .or.  (Iin%para(1) .gt.self%ns)      ) then
        write(self%getprint(),*)"ERROR: 20171227";stop
      endif
      !-----------------------------------------------------------------------
      !Get PC position of J cluster
      Rj = self%LaC%GetPcPos(PCindex) + Inter%P
      !-----------------------------------------------------------------------
      !decompose j
      call self%LaC%GetDecomposeR(Rj,xj,bj)
      if (sum(abs(xj))==0)then
        !-----------------------------local term ---------------
        !Even j is in LC, its index may be different from input
        !  get PC id
          Iin%P = Iin%P + self%LaC%GetPcPos(PCindex)
          PCid  = self%LaC%GetPCidFromPos( Iin%P )
          Iin%P = 0 ! reset
        ! get site index in LC
          Iin%para(2) = self%LaC%GetSiteIdFromPC(PCid,  inter%para(2)  )
        !---append
        self%NHin = self%NHin + 1
        if ( self%NHin.gt.size(self%Hin) )then
          write(self%getprint(),*)"self%NHi is too small in LH";stop
        endif
        self%Hin( self%NHin ) = Iin
      else
        !-----------------------------None local term-----------
           LCpos = matmul( self%LaC%GetVl()  , xj  )                  ! ;write(*,*)666,LCpos,666
           LCid = GetSurroundLC_IdFromPos(self,LCpos)
           if (LCid.eq.-1)then
             self%NLC = self%NLC + 1
             LCid = self%NLC
               !------------
               if (LCid .gt.NSurroundedCLuster )then
                  write(self%getprint(),*)"ERROR: NSurroundedCLuster is too small";stop
               endif
            self%LaPos(:,LCid ) =  LCpos
            self%LaPosR(:,LCid) = self%LaC%GetRealPosFromPCBsis(LCpos)
           endif
           !------------
           PCid = self%LaC%GetPCidFromPos(bj)
           Iin%para(2) = self%LaC%GetSiteIdFromPC(PCid,  inter%para(2)  )!;write(*,*)Iin%para(1),Iin%para(2)
           !-------------append inter
           self%NHou(LCid) = self%NHou(LCid) + 1                        ! ;write(*,*)LCid,self%NHou(LCid)
           if (self%NHou(LCid).gt.  size(self%Hou(:,LCid))  )then
             write(self%getprint(),*)"Size of Hou(1-th,*) is too small";stop
           endif
           self%Hou( self%NHou(LCid)  ,LCid) = Iin
      endif
    endsubroutine




  subroutine ExtendNegativePartH(self)
    implicit none
    class(LH),intent(inout)         :: self
    !--------------------------------------
    integer::NLChalf,jc,iNew,jc2

    NLChalf = self%NLC
    self%NLC  = self%NLC  * 2
    if (self%NLC.gt.NSurroundedCLuster)then
      write(self%getprint(),*)"ERROR: NSurroundedCLuster is too small, recently it=",NSurroundedCLuster
      write(self%getprint(),*)"The number we may need is ",self%NLC
      stop
    endif

    do jc = 1 , NLChalf
       iNew = jc + NLChalf
       self%LaPos(:,iNew)  =  - self%LaPos(:,jc)
       self%LaPosR(:,iNew) =  - self%LaPosR(:,jc)
       self%NHou(iNew)     =    self%NHou(jc)
       do jc2 = 1 , self%NHou(jc)
          self%Hou(jc2,iNew) = self%Hou(jc2,jc)
          self%Hou(jc2,iNew)%para(1) = self%Hou(jc2,jc)%para(2)
          self%Hou(jc2,iNew)%para(2) = self%Hou(jc2,jc)%para(1)
          self%Hou(jc2,iNew)%v = conjg(self%Hou(jc2,jc)%v)
       enddo
    enddo


  endsubroutine


  ! return -1 if it is not existed.
  integer function GetSurroundLC_IdFromPos(self,Pos)
    implicit none
    class(LH),intent(inout) :: self
    integer,intent(in)      :: Pos(3)
    !-------------------------------------
    integer::jc
    GetSurroundLC_IdFromPos = -1
    do jc = 1 , self%NLC
       if (  sum(abs(pos - self%lapos(:,jc))) .eq.0 ) then
         GetSurroundLC_IdFromPos = jc
         goto 999
       endif
    enddo
999 continue
  endfunction






  subroutine SetValueByDiscription(self,Dis,v)
    use functionalsubs
    implicit none
    class(LH),intent(inout) :: self
    character(DiscLen),intent(in)::Dis
    complex*16,intent(in)::V
    !------------------------------------------
    TYPE(funcsubs)::f
    character(DiscLen)::Cin,Ccheck
    integer::jc,jc1

    Cin = f%get_string_upper(Dis)
    Cin = adjustl(trim(Cin))
    !-------------------------------------------------------
    !  local part
    do jc =1 , self%NHin
      Ccheck = f%get_string_upper(  self%Hin(jc)%Disc  )
      Ccheck = adjustl(trim(Ccheck))
      if (Ccheck.eq.Cin) self%Hin(jc)%v = v
    enddo
    !-------------------------------------------------------
    !  connected part   (first part, positive part)
    do jc = 1 , self%NLC/2
      do jc1 = 1 , self%NHou(jc)
        Ccheck = f%get_string_upper(  self%Hou(jc1,jc)%Disc  )
        Ccheck = adjustl(trim(Ccheck))
        if (Ccheck.eq.Cin) self%Hou(jc1,jc)%v = v
      enddo
    enddo
    !-------------------------------------------------------
    !  connected part   (second part, negative part)
    do jc = self%NLC/2 , self%NLC
      do jc1 = 1 , self%NHou(jc)
        Ccheck = f%get_string_upper(  self%Hou(jc1,jc)%Disc  )
        Ccheck = adjustl(trim(Ccheck))
        if (Ccheck.eq.Cin) self%Hou(jc1,jc)%v = conjg(v)
      enddo
    enddo

  endsubroutine

  !
  ! !  add one signle particle term into matrix
  ! !  NOTICE: Matrix = Matrix + data
  ! !  recognized nteracting:
  ! !     [
  ! !     "SpinOnSite"  ,  "OnSite"  ,  "SpinHopping"  ,  "Hopping"
  ! !     ]
  ! !  logcal::IsLocal  ->  put on the conjugate term.
  ! subroutine SetMatixOneTerm(ns,Matrix,spini,spinj,data,IsLocal)
  !   use functionalsubs
  !   implicit none
  !   ! class(LH),intent(inout) :: self
  !   integer,intent(in)      :: Ns
  !   complex*16,intent(inout):: Matrix(Ns,Ns)
  !   integer,intent(in)      :: spini,spinj
  !   type(idata),intent(in)  :: data
  !   logical,intent(in)      :: IsLocal
  !   !-----------------------------------
  !   TYPE(funcsubs)::f
  !
  !   select case(adjustl(trim(f%get_string_upper(data%Itype))))
  !   case("SPINONSITE")               !     "SpinOnSite"
  !     if ( (spini==spinj) .and.  (spini==data%Para(3)) ) then
  !        Matrix(data%Para(1),data%Para(1)) = Matrix(data%Para(1),data%Para(1)) + data%v
  !     endif
  !   case("ONSITE")                   !     "OnSite"
  !     if ( (spini==spinj)                              ) then
  !        Matrix(data%Para(1),data%Para(1)) = Matrix(data%Para(1),data%Para(1)) + data%v
  !     endif
  !   case("SPINHOPPING")             !   "SpinHopping"
  !     if ( (spini==spinj).and.  (spini==data%Para(3))  ) then
  !       Matrix(data%Para(1),data%Para(2)) = Matrix(data%Para(1),data%Para(2)) + data%v
  !       if ( IsLocal ) then!------local term
  !         Matrix(data%Para(2),data%Para(1)) = Matrix(data%Para(2),data%Para(1)) + conjg(data%v)
  !       endif
  !     endif
  !   case("HOPPING")                !      "Hopping"
  !     if ( (spini==spinj)                             ) then
  !       Matrix(data%Para(1),data%Para(2)) = Matrix(data%Para(1),data%Para(2)) + data%v
  !       if ( IsLocal )then !------local term
  !         Matrix(data%Para(2),data%Para(1)) = Matrix(data%Para(2),data%Para(1)) + conjg(data%v)
  !       endif
  !     endif
  !   end select
  !
  ! endsubroutine




  function GetLocalHMatix(self,spini,spinj) result(r)
    implicit none
    class(LH),intent(inout) :: self
    integer,intent(in)::spini,spinj
    complex*16::r(self%ns,self%ns)
    !-----------------------------------
    integer::jc
    call self%CheckInitiatedOrStop()
    r = (0._8,0._8)
    do jc = 1 , self%NHin
       call self%Hin(jc)%SetIntoMatix(self%ns,r,spini,spinj,.True.)
    enddo
  endfunction

  function GetSpinSuppresedLocalHMatrix(self,spini,spinj) result(r)
    implicit none
    class(LH),intent(inout) :: self
    complex*16::r(self%ns*2,self%ns*2)
    !-----------------------------------
    integer::spini,spinj,l1,l2,r1,r2
    integer::jc
    do spini = 0 , 1
      do spinj = 0 , 1
         l1 = spini * self%ns + 1
         l2 = l1 + self%ns - 1
         r1 = spinj * self%ns + 1
         r2 = r1 + self%ns - 1
         r(l1:l2,r1:r2) =  GetLocalHMatix(self,spini=spini,spinj=spinj)
      enddo
    enddo
  endfunction




  ! Get near connected T (exp part is not included.)
  function GetNearTMatix(self,i,spini,spinj) result(r)
    implicit none
    class(LH),intent(inout) :: self
    integer,intent(in)::i,spini,spinj
    complex*16::r(self%ns,self%ns)
    !------------------------------------
    integer::jc
    call self%CheckInitiatedOrStop()
    r = (0._8,0._8)
    if ( (i.gt.0) .and. (i.le.self%NLC)   )then
       Do jc = 1 , self%NHou(i)
          ! call SetMatixOneTerm(self%ns,r,spini,spinj,self%Hou(jc,i),.false.)
          call self%Hou(jc,i)%SetIntoMatix(self%ns,r,spini,spinj,.false.)
       enddo
    else
      write(self%getprint(),*)"ERROR: input i =",i,"is an illegal value";stop
    endif
  endfunction

   integer function GetNs(self)
     implicit none
     class(LH),intent(inout)  :: self
     !-------------------------------------
     GetNs = self%LaC%GetNs()
   endfunction

    integer function GetNnearLC(self)
      implicit none
      class(LH),intent(inout)  :: self
      !-------------------------------------
      call self%CheckInitiatedOrStop()
      GetNnearLC = self%NLC
    endfunction

    function GetNearLCpos(self,i) result(r)
      implicit none
      class(LH),intent(inout)  :: self
      integer,intent(in)       :: i
      real*8::r(3)
      !---------------------------------
      call self%CheckInitiatedOrStop()
      if ( i>0 .and. i<= self%NLC )then
        r = self%LaposR(:,i)
      else
        write(self%getprint(),*)"ERROR: input i should in range [1,self%NLC] where self%NLC=",self%NLC
        write(self%getprint(),*)"Recently i=",i   ;stop
      endif
    endfunction


    function GetNearTexpMatrix(self,q,i,spini,spinj) result(r)
      implicit none
      class(LH),intent(inout)  :: self
      real*8,intent(in)::q(3)
      integer,intent(in)::i,spini,spinj
      complex*16::r(self%ns,self%ns)
      !-------------------------------------
      call self%CheckInitiatedOrStop()

      r = GetNearTMatix(self,i,spini,spinj) * Zexp( sum(self%LaPosR(:,i) * q) * (0._8,1._8)   )
    endfunction

    function GetTqMatrix(self,q,spini,spinj) result(r)
      implicit none
      class(LH),intent(inout)  :: self
      real*8,intent(in)::q(3)
      integer,intent(in)::spini,spinj
      complex*16::r(self%ns,self%ns)
      !-------------------------------------
      integer::jc
      r = (0._8,0._8)
      do jc = 1 , self%NLC
         r = r + GetNearTexpMatrix(self,q,jc,spini,spinj)
      enddo
    endfunction

    function GetSpinSuppresedTq(self,q) result(r)
      implicit none
      class(LH),intent(inout)  :: self
      real*8,intent(in)::q(3)
      complex*16::r(self%ns*2,self%ns*2)
      !-------------------------------------
      integer::spini,spinj,l1,l2,r1,r2
      integer::jc
      do spini = 0 , 1
        do spinj = 0 , 1
           l1 = spini * self%ns + 1
           l2 = l1 + self%ns - 1
           r1 = spinj * self%ns + 1
           r2 = r1 + self%ns - 1
           r(l1:l2,r1:r2) =  GetTqMatrix(self,q,spini=spini,spinj=spinj)
        enddo
      enddo
    endfunction




    subroutine AppendLocalDataToHam(self,H)
      implicit none
      class(Lh),intent(inout)::self
      class(Ham),intent(inout)::H
      !----------------------------------------
      integer::n,jc
      do jc = 1 , self%NHin
        Call self%Hin(jc)%AppendToHam(H)
      enddo
    endsubroutine

    function GetVbasis(self) result(r)
      implicit none
      class(Lh),intent(inout)::self
      real*8::r(3,3)
      !--------------------------------
      r = self%LaC%GetVlReal()
    endfunction


    ! subroutine Report(self,wtp)
    !   implicit none
    !   class(Lh),intent(inout)::self
    !   integer,intent(in)::wtp
    !   !------------------------------
    !   Call self%PrH%report(wtp)
    ! endsubroutine



   ! mode = 0 : good looking
   ! mode = 1 : easy to read from file
    subroutine report(self,wtp,mode)
      implicit none
      class(Lh),intent(inout)::self
      integer,intent(in)::wtp
      integer,intent(in),optional::mode
      !---------------------------------
      integer::jc,jc1,Pmode
      if (present(mode)) then
        Pmode = mode
      else
        Pmode = 0
      endif
      write(wtp,"(A67)",advance='no')"    +=============================================================+"
      write(wtp,*)" "
      write(wtp,'(A68)')"    |                 Hamitonian parameters:                      | "
      do jc = 1 , self%NHin
         do jc1 = 1 , jc - 1
            if ( SELF%Hin(JC1)%Disc== SELF%Hin(JC)%Disc) goto 100
         enddo
         call reportone(wtp,SELF%Hin(JC)%Disc,SELF%Hin(JC)%v,Pmode)
     100 continue
      enddo
    contains

      subroutine reportone(wtp,disc,v,Pmode)
        implicit none
        integer::wtp
        character(32)::disc
        complex*16::v
        integer::Pmode
        !------------------------------
        select case(Pmode)
        case(0)
          write(wtp,100,advance='no')"    |",trim(adjustl(disc)) ," = " ,real(v), " , " ,imag(v),"   | "
          write(wtp,"(A)")" "
        case(1)
          write(wtp,99 ,advance='no')"    |",trim(adjustl(disc)) ," =                           |"
          write(wtp,"(A)")" "
          write(wtp,"( A20 , ES22.14 )")"         RealPart = ",real(v)
          write(wtp,"( A20 , ES22.14 )")"         ImagPart = ",Imag(v)
        case default
          write(*,*)"Pmode in report@lalatticeH.f90 error"
        endselect
        write(wtp,'(A68)',advance='no')"    +.............................................................+ "
        write(wtp,"(A)")" "
     99 format( A5  ,A32 , A30 )
    100 format( A5 , A32 ,A3 ,  ES10.3  ,A3,  ES10.3 ,A5)
    101 format( A5 , A32 )
      endsubroutine

    endsubroutine











    subroutine AbsorbMeanFiled(self,idataarry)
      implicit none
      class(Lh),intent(inout)::self
      class(idata),intent(in)::idataarry(:)
      !---------------------------------------
      type(idata),allocatable::temp(:)
      integer::lin
      lin = size(idataarry)
      allocate(temp(self%NHin))
      temp = self%Hin(1:self%NHin)
      deallocate(self%Hin)
      allocate(self%Hin(self%NHin+lin))
      self%Hin(1:self%NHin) = temp
      self%Hin(self%NHin+1:self%NHin+lin) = idataarry
      self%NHin = self%NHin + lin
      DEALLOCATE(temp)
    endsubroutine


      function GetVp(self) result(r)
        implicit none
        class(LH),intent(in)   :: self
        real*8::r(3,3)
        !-------------------------------
        r = self%LaC%GetVp()
      endfunction


        function GetLatticeConfigPointer(self) result(r)
          implicit none
          class(LH),intent(in) :: self
          class(LaCon),pointer::r
          !---------------------------------------
          call self%CheckInitiatedOrStop()
          r => self%Lac
        endfunction


        logical function GetMeanFieldState(self)
          implicit none
          class(LH),intent(inout) :: self
          !---------------------------------------
          GetMeanFieldState = self%HaveMeanFiled
        endfunction

        subroutine SetMeanFieldState(self,tf)
          implicit none
          class(LH),intent(inout) :: self
          logical,intent(in)::tf
          !---------------------------------------
          self%HaveMeanFiled = tf
        endsubroutine

        subroutine SetMeanField(self,M,spini,spinj)
          implicit none
          class(LH),intent(inout) :: self
          complex*16,intent(in)::M(self%ns,self%ns)
          integer,intent(in)::spini,spinj
          !---------------------------------------
          self%DeltaM(:,:,spini,spinj) = M
        endsubroutine

        function getMeanField(self,spini,spinj) result(r)
          implicit none
          class(LH),intent(inout) :: self
          integer,intent(in)::spini,spinj
          complex*16::r(self%ns,self%ns)
          !---------------------------------------
          r = self%DeltaM(:,:,spini,spinj)
        endfunction


        function GetVelocityOneCluster(self,q,spini,spinj,xy,i)  result(r)! i = [0,self%GetNnearLC()]
          implicit none
          class(LH),intent(inout) :: self
          real*8,intent(in)       :: q(3)
          integer,intent(in)      :: spini,spinj,xy,i
          complex*16              :: r(self%Ns,self%Ns)
          !---------------------------------------------------------------------
          real*8::CR(3),LoX(self%ns,self%ns)
          integer::jc1,jc2
          if (i==0)then
            r  = self%GetLocalHMatix(spini,spinj) - self%getMeanField(spini,spinj)
            CR = 0._8
          else
            r  = self%GetNearTexpMatrix(q,i,spini,spinj)
            CR = self%LaposR(:,i)
          endif
          r = r * ( self%DeltaX(:,:,xy) - CR(xy) )
        endfunction


        function GetVelocityMatrix(self,q,spini,spinj,xy) result(r)
          implicit none
          class(LH),intent(inout) :: self
          real*8,intent(in)       :: q(3)
          integer,intent(in)      :: spini,spinj,xy
          complex*16              :: r(self%Ns,self%Ns)
          !---------------------------------------------------------------------
          integer::jc
          r = (0._8,0._8)
          do jc = 0 , self%GetNnearLC()
             r = r + GetVelocityOneCluster(self,q,spini,spinj,xy,jc)
          enddo
        endfunction


    ! function GetHamList(self) result(r)
    !   implicit none
    !   class(VCAdH),intent(inout)::self
    !   TYPE(Ham)::r
    !   !----------------------------------------
    !   type(idata)::idataarray(Nmax)
    !   integer::n,jc,para(8)
    !
    !   call r%Initialization( self%ns  , self%getprint() )
    !
    !   call GetTotalIdataArray(self,idataarray,n)
    !
    !   call r%StartAppendingInteraction()
    !   do jc = 1 , n
    !     para = idataarray(jc)%Para
    !     para(1:2) = para(1:2) - 1
    !     call r%AppendingInteraction(  idataarray(jc)%Itype  , Para , idataarray(jc)%v)
    !   enddo
    !   call r%EndAppendingInteraction()
    !
    ! endfunction












  !
  ! subroutine test(self)
  !   implicit none
  !   class(LH),intent(inout) :: self
  !   !-----------------------------------
  !   integer::jc,jc1,jc2,id
  !   real*8::q(3)
  !   complex*16::H(self%ns,self%ns)
  !
  !   q(1) = 0.6; q(2) = 0.9 ; q(3) = 0._8   ;q=0._8
  !
  !   H = self%GetTqMatrix(q,spini=0,spinj=0)
  !   write(*,*) H - transpose(conjg(H))
  !
  !
  !
  !   ! H = self%GetNearTMatix(1,0,0)
  !   ! do jc1 =1 , 4 ; do jc2 = 1 , 4
  !   !   if (abs(H(jc1,jc2)).gt.0.0001 ) write(*,*)jc1,jc2,H(jc1,jc2)
  !   ! enddo ;enddo
  !
  !
  !
  ! endsubroutine
  !



endmodule
