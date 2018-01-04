

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : LatticeConfig
! OBJECT: TYPE(LaCon)
! USED  : CubicLatticeCell
! DATE  : 2017-12-23
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            Fully discribe a geometry information of a lattice system.
!
! STANDARD:
!            *CALL Initialization(Vp,PC,PCi,Vl,PRINT_,show_)
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  [sub] Initialization(Vp,PC,PCi,Vl,PRINT_,show_)
!                        real*8 ::Vp(3,3)    the basis of the primary cell.
!                        real*8 ::PC(3,Np)   Np is the number of sites in Primary cell. PC(3,i) is the position
!                        integer::PCi(Np)    orbital index.          !1-basis
!                        integer::Vl(3,3)    use Vp as basis, the three vector of lattice cell
!
! avalable gets:
!                   [fun] GetNs()
!                         integer       get total number of site in the Lattice Cell
!
!                   [fun] GetSiteP(i)
!                         return integer::r(3)  is the position of the site
!
!                   [fun] GetSiteRealP(i)
!                         return real*8::r(3)  is the position of site in real space.
!
!                   [fun] GetOrbitIndex(i):
!                         integer::
!                         return the orbital index of site i
!
!                   [fun] GetPcPos(i)
!                         return integer::r(3) is the position of i-th PC.
!
!                   [fun] GetPCidFromPos(p)
!                         integer::p(3),GetPCidFromPos
!                         return the PC id from its position.
!
!                   [sub] GetDecomposeR(R,p,x)
!                         integer::R(3),p(3),x(3)
!                         decompose   R = Vl . p + X    where Vl is the basis of LC.
!
!
!                   [fun] GetNp()
!                         integer    get the total number of primary cell
!
!                   [fun] GetSiteIdFromPC(Pcid,SiteInPC)
!                         for a given PC index, and also a site index in one PC, return the site id in lattice.
!
!                   [fun] GetVl()
!                         return integer::vl(3,3)
!
!                   [fun] GetVlReal()
!                         return real*8::vl(3,3)  (in real space.)
!
!                   [fun] GetRealPosFromPCBsis(p)
!                         return real*8::r(3)
!                         for a input integer::p(3) -> represent a position by using PC basis
!                         return the position of p in real space.
!
!
! avalable is :
!                  [fun] IsInitiated
! others      :
!                  ![sub] p
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module LatticeConfig
  use CubicLatticeCell
  implicit none


  type::LaCon
    private
    logical             :: initiated = .false.
    !-------------
    real*8              :: Vp(3,3)   ! basis of PC
    integer             :: Nc        ! number of sites in one PC
    real*8,allocatable  :: pc(:,:)   ! save all sites position in one PC
    integer,allocatable :: PCi(:)    ! save all orbital index
    !-------------
    integer             :: Np        ! number of PC in LC
    integer,allocatable :: Pp(:,:)   ! position of these PCs
    !-------------
    integer             :: Vl(3,3)   ! basis of LC  (in the basis of PC)
    integer             :: Ns        ! total number of sites in the lattice cell.
    real*8,allocatable  :: lc(:,:)   ! all the position of sites in LC
    integer,allocatable :: lci(:)    ! all orbital index
    !-------------
    TYPE(CubLattCe)     :: Cubic
    !=========
    integer::print = 6
    integer::show  = 0

  contains
    procedure,pass::Initialization
    final::Finalization

    procedure,pass::GetNs
    procedure,pass::GetSiteP
    procedure,pass::IsInitiated

    procedure,pass::GetNp
    procedure,pass::GetSiteIdFromPC
    procedure,pass::GetPcPos
    procedure,pass::GetPCidFromPos
    procedure,pass::GetDecomposeR
    procedure,pass::GetVl
    procedure,pass::GetRealPosFromPCBsis
    procedure,pass::GetSiteRealP
    procedure,pass::GetOrbitIndex
    procedure,pass::GetVlReal
  endtype


  private::Initialization,UnInitialization,Finalization

  private::CheckInitiatedStop
  private::GetNs,GetSiteP

  private::IsInitiated

  private::GetNp
  private::GetSiteIdFromPC
  private::GetPcPos,GetPCidFromPos
  private::GetDecomposeR
  private::GetVl
  private::GetRealPosFromPCBsis
  private::GetSiteRealP
  private::GetOrbitIndex
  private::GetVlReal

contains



  subroutine Initialization(self,Vp,PC,PCi,Vl)
    implicit none
    class(LaCon),intent(inout)::self
    real*8,intent(in) ::Vp(3,3)
    real*8,intent(in) ::PC(:,:)
    integer,intent(in)::PCi(:)
    integer,intent(in)::Vl(3,3)
    !--------------------------------------------------
    integer::jc,jc1,jc2
    real*8 ::R(3)

    call UnInitialization(self)   ;  self%initiated = .true.


    self%vp = vp
    !--------------------------------------------------------------------------------
    !   allocate  PC
      self%Nc = size(PC(1,:))
      if (size(PC(:,1)).ne.3)then
         write(self%print,*)"ERROR: format of input PC in LaCon is wrong." ; stop
      endif
      allocate(self%pc , source = PC )
      allocate(self%pci, source = pci)
    !-------------------------------------------------------------------------------
    !  allocate Cubic
       self%vl = vl
       call self%Cubic%Initialization( Vl )
    !-------------------------------------------------------------------------------
    !  allocate position of PCs
       self%Np = self%Cubic%GetNc()
       allocate(   self%Pp(3,self%Np)   )
       do jc = 1 , self%Np
          self%pp(:,jc) = self%Cubic%GetSite(jc)
       enddo
    !-------------------------------------------------------------------------------
       self%Ns = self%Np * self%Nc
       allocate(  self%lc(3,self%Ns)    )
       allocate(  self%lci(self%Ns)     )
       jc2 = 0
       do jc = 1 , self%Np
         R = matmul(  self%Vp  ,  self%Pp(:,jc)*1._8  )
          do jc1 = 1 , self%Nc
            !------------------------------------
            !  jc2 = jc2 + 1
            !------------------
            jc2 = GetSiteIdFromPC(self,jc,jc1)    !;write(*,*)jc2,jc,jc1
            !-------------------------------------
             self%lc(:,jc2) = self%pc(:,jc1) + R
             self%lci( jc2) = self%pci(jc1 )
          enddo
       enddo
    !-------------------------------------------------------------------------------


  endsubroutine



  subroutine UnInitialization(self)
    implicit none
    class(LaCon),intent(inout)::self
    !-------------------------------------
    if (self%initiated)then ; self%initiated = .false.
       deallocate(   self%pc  , self%pp  , self%lc ,self%lci,  self%pci    )
    endif
  endsubroutine


  impure elemental subroutine Finalization(self)
    implicit none
    type(LaCon),intent(inout)::self
    !-------------------------------------
    call UnInitialization(self)
  endsubroutine


  subroutine CheckInitiatedStop(self)
    implicit none
    class(LaCon),intent(inout)::self
    !-------------------------------------
    if (.not.self%initiated)then
       write(self%print,*)"LaCon is not initiated yet" ; stop
    endif
  endsubroutine



  integer function GetNs(self)
    implicit none
    class(LaCon),intent(inout)::self
    !-------------------------------------
    call CheckInitiatedStop(self)
    GetNs = self%Ns
  endfunction

  integer function GetNp(self)
    implicit none
    class(LaCon),intent(inout)::self
    !-------------------------------------
    call CheckInitiatedStop(self)
    GetNp = self%Np
  endfunction




  function GetSiteP(self,i) result(r)
    implicit none
    class(LaCon),intent(inout)::self
    integer,intent(in)::i
    real*8::r(3)
    !-------------------------------------
    call CheckInitiatedStop(self)
    r = self%lc(:,i)
  endfunction


  integer function GetSiteIdFromPC(self,PCid,Sid)
    implicit none
    class(LaCon),intent(inout)::self
    integer,intent(in)::PCid,Sid
    !--------------------------------------------
    GetSiteIdFromPC = self%Nc * ( PCid - 1 ) + Sid
  endfunction


  logical function IsInitiated(self)
    implicit none
    class(LaCon),intent(in)::self
    !-------------------------------------
    IsInitiated = self%initiated
  endfunction

  function GetPcPos(self,i) result(r)
    implicit none
    class(LaCon),intent(inout)::self
    integer,intent(in)::i
    integer::r(3)
    !-------------------------------------
    r = self%pp(:,i)
  endfunction

  integer function GetPCidFromPos(self,p)
    implicit none
    class(LaCon),intent(inout)::self
    integer,intent(in)::p(3)
    !------------------------------------
    GetPCidFromPos = self%Cubic%GetIdFromP(p)
  endfunction


  subroutine GetDecomposeR(self,R,P,X)
    implicit none
    class(LaCon),intent(inout)::self
    integer,intent(in)::R(3)
    integer,intent(out)::P(3),x(3)
    !------------------------------------
    call self%Cubic%DecomposeR(R,p,x)
  endsubroutine

  Function GetVl(self) result(r)
    implicit none
    class(LaCon),intent(inout)::self
    integer::r(3,3)
    !-----------------------------------
    r = self%vl
  endfunction

  Function GetVlReal(self) result(r)
    implicit none
    class(LaCon),intent(inout)::self
    real*8::r(3,3)
    !-----------------------------------
    r = matmul(self%vp , self%vl*1._8)
  endfunction



  function GetRealPosFromPCBsis(self,p) result(r)
    implicit none
    class(LaCon),intent(inout)::self
    integer,intent(in)::p(3)
    real*8::r(3)
    !-----------------------------------
    r = matmul(  self%Vp   , p * 1._8   )
  endfunction


  function GetSiteRealP(self,i) result(r)
    implicit none
    class(LaCon),intent(inout)::self
    integer,intent(in)::i
    real*8::r(3)
    !-----------------------------------
    r = self%lc(:,i)
  endfunction


  integer function GetOrbitIndex(self,i)
    implicit none
    class(LaCon),intent(inout)::self
    integer,intent(in)::i
    !-----------------------------------
    GetOrbitIndex = self%lci(i)
  endfunction





endmodule
