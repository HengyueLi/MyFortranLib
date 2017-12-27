

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : CubicLatticeCell
! OBJECT: TYPE(CubLattCe)
! USED  : class_numerical_method
! DATE  : 2017-12-22
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            For a cubic lattice of dimension 3, assume the length is 1, get information
!        of a lattice cell we chosen.
!
! STANDARD:
!            *CALL Initialization(A,PRINT_,show_)
!
!
!          5 ────────────────────┐7
!           ╱┊                  ╱│
!          ╱ ┊                 ╱ │
!         ╱  ┊                ╱  │
!        ╱   ┊               ╱   │
!       ╱    ┊ 2            ╱    │
!      ╱     .╌╌╌╌╌╌╌╌╌╌╌╌╌╱╌╌╌╌╌╱  6
!     3 ─────────────────┐4     ╱
!     │   .              │     ╱
!     │  .               │    ╱
!     │ .                │   ╱
!     │.                 │  ╱
!     │                  │ ╱
!    0 ──────────────────── 1

!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  [sub] Initialization(A,PRINT_,show_)
!                        integer::A(3,3)
!                        where A(:,i) represent the i-th basis of lattice cell.
!
! avalable gets:
!                   [fun] GetNc()
!                         integer,  return the total number of site in the CellLattice
!
!                   [fun] GetSite(i)
!                         return integer::GetSite(3) is the position of the i-th site in the CellLattice
!
!                   [fun] GetIdFromP(p)
!                         input the position of P(3) , return the index in the CellLattice
!
!                   [sub] DecomposeR(R,x,b)
!                         integer::R(3),x(3),b(3)
!                         for any given vector R  we can compose it as :
!                              R = A . x + b
!
! avalable is :
!                  ![fun] i
! others      :
!                  ![sub] p
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module CubicLatticeCell

  implicit none


  type::CubLattCe
    private
    logical             :: Initiated = .false.
    integer             :: A(3,3)
    integer             :: Nc     ! number os sites in the shape
    integer,allocatable :: P(:,:) ! store all the sites in  the shape.
    !------------------------------------------------------------------
    real*8:: shirink = 0.00000001_8  ! (ratio)
    !============================
    integer :: print = 6
    integer :: show  = 0


  contains
    procedure,pass::Initialization
    final::Finalization

    procedure,pass::GetNc
    procedure,pass::GetSite
    procedure,pass::GetIdFromP
    procedure,pass::DecomposeR
  endtype
  private::Initialization,UnInitialization,Finalization

  private::GetNc
  private::Det
  private::check_point_in_volume
  private::DecomposeRbyBasis
  private::Find_All_point_closed_by_volume
  private::CheckInitiatedStop
  private::GetSite
  private::GetIdFromP
  private::DecomposeR

contains


  subroutine Initialization(self,A)
    implicit none
    class(CubLattCe),intent(inout)::self
    integer,intent(in)::a(3,3)
    !-----------------------------------
    call UnInitialization(self) ; self%Initiated = .true.


    self%A = a
    !-------------------------------------
    !         set volume
      self%nc = abs(Det(a(:,1),a(:,2),a(:,3)))     
    !-------------------------------------
      allocate( self%p(3,self%nc)  ) ; self%p = 0
    !-------------------------------------
    !  searching points
      call Find_All_point_closed_by_volume(self)

  endsubroutine

  subroutine UnInitialization(self)
    implicit none
    class(CubLattCe),intent(inout)::self
    !----------------------------------------
    if (self%Initiated)then ; self%Initiated = .false.
       deallocate(self%p)
    endif
  endsubroutine

  impure elemental subroutine Finalization(self)
    implicit none
    type(CubLattCe),intent(inout)::self
    !----------------------------------------
    call UnInitialization(self)
  endsubroutine


  subroutine CheckInitiatedStop(self)
    implicit none
    class(CubLattCe),intent(inout)::self
    !----------------------------------------
    if (.not.self%Initiated)then
      write(self%print,*)"ERROR: CubicLatticeCell is not initiated yet."
    endif
  endsubroutine


  !  make a 3D volume by vector a1,a2,a3,  check if point r is in this volume or not.
  ! they are all in a 3D Orthonormal basis.
  ! zero is a small value used to treat the case when the point is on the length of the volume.
  logical function check_point_in_volume(self,r)
    implicit none
    class(CubLattCe),intent(inout)::self
    real*8,intent(in)::r(3)
    !-------------------------------------------
    real*8::a(3,3),RP(3)
    INTEGER::JC
    a = self%A * ( 1._8 - self%shirink )
    RP = DecomposeRbyBasis( R , a )     !
    check_point_in_volume = .TRUE.
    DO JC = 1,3
       check_point_in_volume = check_point_in_volume .AND. ( RP(JC).LE.1._8 ).AND. ( RP(JC).GE.0._8 )
    ENDDO
  endfunction

 ! decompose vector R by basis b
  function DecomposeRbyBasis(R,b) result(a)
    use class_numerical_method
    implicit none
    real*8,intent(in)::R(3),b(3,3)
    real*8::a(3)
    !--------------------------------------
    TYPE(nummethod)::f
    real*8::M(3,3)
    M = b
    call f%MatrixInverse( 3 , M )
    a = MATMUL( M , R )
  endfunction

  ! subroutine SearchingAllInclosedPoints(self)
  !   implicit none
  !   class(CubLattCe),intent(inout)::self
  !   !-------------------------------------------
  !   real*8::
  ! endsubroutine


  ! for a given basis, the volume is represented by three integer vector
  ! find out all the posible point(which are also integer) in it.
  subroutine Find_All_point_closed_by_volume(self)
    implicit none
    class(CubLattCe),intent(inout)::self
    !-----------------------
    integer::np
    integer::a1(3),a2(3),a3(3)
    integer::ScaleOfC(3,0:7),pmin(3),pmax(3),jc,jc1,jc2,jc3
    real*8::r(3),checkboudary


    checkboudary = ( 1._8 - self%shirink )
    !
    a1 = self%a(:,1)
    a2 = self%a(:,2)
    a3 = self%a(:,3)
    ScaleOfC(:,0) = 0
    ScaleOfC(:,1) = a1
    ScaleOfC(:,2) = a2
    ScaleOfC(:,3) = a3
    ScaleOfC(:,4) = a1 + a2
    ScaleOfC(:,5) = a1 + a3
    ScaleOfC(:,6) = a2 + a3
    ScaleOfC(:,7) = a1 + a2 + a3
    !-----------------------------------
    do jc = 1,3
      pmin(jc) = min(ScaleOfC(jc,0),ScaleOfC(jc,1),ScaleOfC(jc,2),ScaleOfC(jc,3),&
                 ScaleOfC(jc,4),ScaleOfC(jc,5),ScaleOfC(jc,6),ScaleOfC(jc,7)     )
      pmax(jc) = max(ScaleOfC(jc,0),ScaleOfC(jc,1),ScaleOfC(jc,2),ScaleOfC(jc,3),&
                 ScaleOfC(jc,4),ScaleOfC(jc,5),ScaleOfC(jc,6),ScaleOfC(jc,7)     )
    enddo
                                  !;write(*,*)a2!,a2,a3
    np = 0
    do jc1= pmin(1),pmax(1)
      do jc2 = pmin(2),pmax(2)
        do jc3 = pmin(3),pmax(3)
           r(1) = jc1*1._8 ;r(2) = jc2*1._8  ;r(3) = jc3*1._8
           if (check_point_in_volume(self,r))then
             np = np + 1
             if (np.gt.self%Nc)then
               write(self%print,*)"ERROR: A number of included sites is larger than volume";stop
             endif
             self%p(1,np) = jc1 ;self%p(2,np) = jc2 ;self%p(3,np) = jc3 !; write(*,*)jc3
           endif
        enddo
      enddo
    enddo

  endsubroutine


  integer function GetNc(self)
    implicit none
    class(CubLattCe),intent(inout)::self
    !----------------------------------
    call CheckInitiatedStop(self)
    GetNc = self%Nc
  endfunction

  function GetSite(self,i) result(r)
    implicit none
    class(CubLattCe),intent(inout)::self
    integer,intent(in)::i
    integer::r(3)
    !-----------------------
    call CheckInitiatedStop(self)
    r = self%p(:,i)
  endfunction

  integer function GetIdFromP(self,p)
    implicit none
    class(CubLattCe),intent(inout)::self
    integer,intent(in)::p(3)
    !-----------------------
    integer::jc
    call CheckInitiatedStop(self)
    GetIdFromP = -1
    do jc = 1 , self%Nc
       if (  sum(abs(self%p(:,jc) - p)) .eq.  0  )then
         GetIdFromP = jc
         goto 999
       endif
    enddo
999 continue
  endfunction


  subroutine DecomposeR(self,p,xi,r)
    implicit none
    class(CubLattCe),intent(inout)::self
    integer,intent(in) ::  p(3)
    integer,intent(out):: xi(3)   , r(3)
    !-----------------------------------
    integer::jc
    real*8::rlcs(3,3) , rp(3)
    rlcs = self%A * 1._8
    rp   = p      * 1._8
    rp = DecomposeRbyBasis(rp,rlcs)
    do jc = 1,3
      if (rp(jc).ge.0._8)then
          rp(jc) = rp(jc) + self%shirink
      else
         if ( abs(rp(jc)-int(rp(jc))).gt.self%shirink ) rp(jc) = rp(jc) -1._8
      endif
    enddo
    xi = int(rp)
    r  = p - matmul( self%A , xi )
  endsubroutine







  integer function Det(v1,v2,v3)
    implicit none
    integer,intent(in)::v1(3),v2(3),v3(3)
    !------------------------------------------
    Det = &
         v1(1) * v2(2) * v3(3) - v1(1) * v2(3) * v3(2) + &
         v1(2) * v2(3) * v3(1) - v1(2) * v2(1) * v3(3) + &
         v1(3) * v2(1) * v3(2) - v1(3) * v2(2) * v3(1)
  endfunction

  ! integer function Vof



endmodule
