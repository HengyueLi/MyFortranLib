!
! !!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! ! TYPE  : MODULE
! ! NAME  : VCA_Types
! ! OBJECT: TYPE(VCATP)
! ! USED  :
! ! DATE  : 2018-01-05
! ! AUTHOR: hengyueli@gmail.com
! !--------------
! ! Open-Source : No
! !------------------
! ! DESCRIPTION:
! !              Use to Initialize VCA.
! !
! ! STANDARD:
! !
! !
! ! USING LIST:
! !            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! !
! !
! ! avalable sets:
! !
! !
! ! avalable gets:
! !
! !
! ! avalable is :
! !
! ! others      :
! !
! !
! !
! !
! !!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!
! module VCA_Types
!   implicit none
!
!
!
!   type,public::VCA_Wald_Para
!     !=====================================================
!     ! Wald functional parameters
!     integer::jobi(25)
!     real*8 ::jobr(25)
!     !=====================================================
!     ! solver parameters
!     character(2):: SolverType  = "ED"        ! = "ED" / "LA"
!     logical     :: IsReal      = .false.
!     integer     :: symmetry    =  0
!     real*8      :: Temperature =  0.000001_8
!     !--------------------------
!     ! ED use
!     !
!     !--------------------------
!     !  Lanczos use
!     real*8  :: pre      = 1.e-14
!     real*8  :: bzero    = 0.000001
!     integer :: M        = 30
!     logical :: oth      = .true.
!     real*8  :: DegPre   = 1.e-6
!     !=====================================================
!     ! Green's function parameters
!     !---------------------------
!     ! ED
!     !---------------------------
!     ! LANCZOS
!     integer :: GM     = 90
!     logical :: Goth   = .true.
!     real*8  :: Gbzero = 1.e-7
!   endtype
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
!   integer,parameter,private::NMaxPCsite = 20
!
!   Type,public::VCA_Latt_Para
!     integer:: NSiteInPC                    ! number of sites in Primary Cell (PC)
!     real*8 :: PCBasis(3,3 )                ! the basis of the PC.
!     real*8 :: SitePos(3,NMaxPCsite)        ! SitePos(:,i) is the position of site in PC.
!     integer:: OrbIndx(NMaxPCsite  ) = 1    ! for each site, we can give a orbital index (if needed).
!     integer:: LCBasis(3,3 )                ! In PCBasis representation, the basis of chosen Lattice Cell(LC).
!
!   endtype
!
!
!   type,public::VCA_DelH_Para
!     character(32):: Type
!     character(32):: Discription
!     integer      :: jobi(20)
!     real*8       :: jobr(20)
!   endtype
!
!
!
!
!
!
! endmodule
