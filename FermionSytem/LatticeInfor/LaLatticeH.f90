


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : LaLatticeH
! OBJECT: TYPE(LH)
! USED  : CodeObject,LaPrimaryH,LatticeConfig,FermionHamiltonian,functionalsubs
! DATE  : 2017-12-25
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
! avalable gets:
!                   [fun] G
!
! avalable is :
!                  ![fun] i
! others      :
!                  ![sub] p
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module LaLatticeH
  use CodeObject
  use LaPrimaryH
  use LatticeConfig
  use LaPrimaryHUsedList    , only : Ilist => ListStru
  use LaprimaryHusedDatatype
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
    ! save all the position of surrounded Lattice Cell. In the basis of PC's
    integer :: LaPos(3,NSurroundedCLuster)
    integer                :: NHou(NSurroundedCLuster)
    type(idata),allocatable:: Hou(:,:)         ! list for saving all Hamiltonian out of Lattice.Inner use

    ! class(Ham),pointer   :: Hm   => NULL()    ! Ham used for cluster solver.

  CONTAINS
    procedure,pass::Initialization
    final::Finalization


    procedure,pass::SetValueByDiscription
  endtype


  private::Initialization,Finalization

  private::SettingLatticeHbyPrimaryH
  private::AppendOneInterFromPrimaryCellToLatticeCell
  private::GetSurroundLC_IdFromPos,ExtendNegativePartH
  private::SetValueByDiscription
  private::SetMatixOneTerm


contains


  subroutine Initialization(self,PrH,LaC,print_,show_)
    implicit none
    class(LH),intent(out)         :: self
    class(PH),target,intent(in)   :: PrH
    class(LaCon),target,intent(in):: Lac
    integer,intent(in),optional   :: print_,show_
    !----------------------------------------------
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
  endsubroutine

  subroutine Finalization(self)
    implicit none
    type(LH),intent(inout)::self
    !-----------------------------
    if (self%IsInitiated())then
       deallocate(   self%Hin  , self%Hou    )
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
        Inter = self%PrH%GetIdata(jci)
        Call AppendOneInterFromPrimaryCellToLatticeCell(self,Inter,jcp)
      enddo
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
    integer::i,j,PCpos(3),PCdx(3),PCdp(3),LCid,LCpos(3),PCid
    logical::IsLClocal


    Iin         = Inter
    !Get fist index, the first index is always in the cluster.
    Iin%para(1) = self%LaC%GetSiteIdFromPC(PCindex,  inter%para(1)  )
    !Get second index
    PCpos = Inter%P + self%LaC%GetPcPos(PCindex)
    PCid = self%LaC%GetPCidFromPos( PCpos )
    if (PCid.eq.-1)then
      IsLClocal = .False.
    else
      Iin%para(2) = self%LaC%GetSiteIdFromPC(PCid,  inter%para(2)  )
      if ( (Iin%para(2) .le. self%LaC%GetNs())  .and.  (Iin%para(2).ge.0)  )then
        IsLClocal = .true.
      else
        IsLClocal = .false.
      endif
    endif

                         !
                         !
                        !  write(*,*)PCID
                        !  write(*,*)Iin%Disc
                        !  write(*,*)iin%Itype
                        !  write(*,*)"PC=",Inter%Para(1:2)
                        !  write(*,*)"LC=",Iin%para(1:2)
                        !  write(*,*)"----------------------"



    if ( IsLClocal ) then
      !---local term
      self%NHin = self%NHin + 1
      if ( self%NHin.gt.size(self%Hin) )then
        write(self%getprint(),*)"self%NHi is too small in LH";stop
      endif
      self%Hin( self%NHin ) = Iin
    else
      !---Non Local---------------------------------
      !  get the postition of connected PC
         PCpos = self%lac%GetPcPos(PCindex) + Iin%p
      !  decompoase Position
         call self%LaC%GetDecomposeR(PCpos,PCdx,PCdp)
      !------------
         LCpos = matmul( self%LaC%GetVl()  , PCdx  )
         LCid = GetSurroundLC_IdFromPos(self,LCpos)
         if (LCid.eq.-1)then
           self%NLC = self%NLC + 1
           LCid = self%NLC
             !------------
             if (LCid .gt.NSurroundedCLuster )then
                write(self%getprint(),*)"ERROR: NSurroundedCLuster is too small";stop
             endif
          self%LaPos(:,LCid) =  LCpos
         endif

         !------------
         PCid = self%LaC%GetPCidFromPos(PCdp)
         Iin%para(2) = self%LaC%GetSiteIdFromPC(PCid,  inter%para(2)  )
         !-------------append inter
         self%NHou(LCid) = self%NHou(LCid) + 1
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
       self%LaPos(:,iNew) = - self%LaPos(:,jc)
       self%NHou(iNew)    =   self%NHou(jc)
       do jc2 = 1 , self%NHou(jc)
          self%Hou(jc2,iNew) = self%Hou(jc2,jc)
          self%Hou(jc2,iNew)%para(1) = self%Hou(jc2,iNew)%para(2)
          self%Hou(jc2,iNew)%para(2) = self%Hou(jc2,iNew)%para(1)
          self%Hou(jc2,iNew)%v = conjg(self%Hou(jc2,iNew)%v)
       enddo
    enddo


  endsubroutine


  ! return -1 if it is not existed.
  integer function GetSurroundLC_IdFromPos(self,Pos)
    implicit none
    class(LH),intent(inout) :: self
    integer,intent(in)    :: Pos(3)
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


  !  add a two particle term into matrix
  !  NOTICE: Matrix = Matrix + data
  subroutine SetMatixOneTerm(self,Matrix,spini,spinj,data)
    use functionalsubs
    implicit none
    class(LH),intent(inout) :: self
    complex*16,intent(inout):: Matrix(self%Ns,self%Ns)
    integer,intent(in)      :: spini,spinj
    type(idata),intent(in)  :: data
    !-----------------------------------
    TYPE(funcsubs)::f

    select case(adjustl(trim(f%get_string_upper(data%Itype))))
    case("SPINONSITE")               !     "SpinOnSite"
      if ( (spini==spinj) .and.  (spini==data%Para(3)) ) then
         Matrix(data%Para(1),data%Para(1)) = Matrix(data%Para(1),data%Para(1)) + data%v
      endif
    case()

    case

    end select

  endsubroutine




  subroutine test(self)
    implicit none
    class(LH),intent(inout) :: self
    !-----------------------------------
    integer::jc
    do jc = 1 , self%NHin
      write(*,*) self%hin(jc)%Disc
      write(*,*) self%hin(jc)%para(1:3)
      write(*,*) self%hin(jc)%Itype
      write(*,*) self%hin(jc)%v
      write(*,*)"------------------"
    enddo
  endsubroutine




endmodule
