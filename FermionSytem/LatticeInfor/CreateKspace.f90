



!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : CreateKspace
! OBJECT: TYPE(Kspace)
! USED  : CodeObject
! DATE  : 2017-08-9
! AUTHOR: Hengyue Li
!--------------
! DESCRIPTION:
!            creat the reciprocal space of a lattice system.
!            step1 *call self%Initialization(a,n,meshtype)
!                  meshtype=0 :creat equal distant k points
!                  meshtype=1 :creat GS k points
!
!
! STANDARD:
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  [sub] Initialization(a,n,meshtype)
!                        real*8::a(3,3)
!                        integer::n(3) ! represent the sampling points in each axis.
!                        integer::meshtype
!                        if one direction is suppresed, set n = 1
!
! avalable gets:
!                  [fun] getnk()
!                        return the total Nk
!
!                  [fun] getk(i)
!                        real*8::getk(3),  return jc-th k point
!
!                  [fun] getw(i)
!                        real*8::getw   return the weight
!
!                  [fun] GetA(i)
!                        return a(:,i) the basis of cell
!
!                  [fun] GetB(i)
!                        return the reciprocal vector
!
!                  [fun] getvolr()
!                        real*8::real    return the volue of unit cell
!
!                  [fun] getvolb()
!                        real*8::real    return the volue of Brillouin zone
!
! avalable  IS  :
!                  [fun]
!
! others        :
!                  [sub]
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$






module CreateKspace
  use CodeObject
  use class_integration
  implicit none

  type,extends(object)::Kspace
    private
    real*8 :: a(3,3)
    real*8 :: b(3,3)
    integer:: n(3)
    integer:: meshtype
    !--------------------
    integer::nk
    real*8,allocatable::k(:,:)
    real*8,allocatable::w(:)

  contains
    procedure,pass::Initialization
    final::Finalization

    procedure,pass::getk
    procedure,pass::getw
    procedure,pass::getnk
    procedure,pass::geta
    procedure,pass::GetB
    procedure,pass::getvolr
    procedure,pass::getvolb
  endtype


  private::Initialization,Finalization
  private::vector3D_multiply,cala_volue
  private::SetKw
  private::CheckCanBeGetOrStop,CheckDimentionGetAllowed
  private::getk,getw,getnk,geta,GetB,getvolr,getvolb

contains

  subroutine Initialization(self,a,n,meshtype,print_,show_)
    implicit none
    class(Kspace),intent(out)   :: self
    real*8,intent(in)           :: a(3,3)
    integer,intent(in)          :: n(3)
    integer,intent(in)          :: meshtype
    integer,intent(in),optional :: print_,show_
    !------------------------------------------------
    call self%SetInitiated(.true.)
    self%a        = a
    self%n        = n
    self%meshtype = meshtype
    if (present(print_)) call self%setprint(print_)
    if (present(show_ )) call self%SetShow( show_ )

    if (  min(  n(1),n(2),n(3) )  .lt.1 )then
       write(self%getprint(),*)"ERROR: min(n)>0 must be kept in Kspace. Now N=",n
       stop
    endif

    self%nk = n(1) * n(2) * n(3)
    allocate(  self%k(3,self%nk)   )
    allocate(  self%w(  self%nk)   )

    call SetKw(self)

  endsubroutine

  subroutine Finalization(self)
    implicit none
    type(Kspace),intent(inout) :: self
    if (self%IsInitiated())then
       call self%SetInitiated(.false.)
       deallocate(  self%k  ,self%w   )
    endif
  endsubroutine

  subroutine SetKw(self)
    implicit none
    class(Kspace),intent(inout) :: self
    !------------------------------------------------
    select case(self%meshtype)
    case(0,1)
      call creatk_mesh_01(self)
    case default
      write(self%getprint(),*)"ERROR: Unknow meshtype in Kspace,mesh=",self%meshtype;stop
    endselect
  endsubroutine




  subroutine vector3D_multiply(a,b,c)
      implicit none
      real(8)::a(3),b(3),c(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)                    ! ;write(*,*)a  ;write(*,*)b ;write(*,*)c ;write(*,*)"---------"
  end subroutine


  real*8 function cala_volue(a1,a2,a3)
          implicit none
          real*8,intent(in)::a1(3),a2(3),a3(3)
          !--------------------------
          real*8::temp(3)
          call vector3D_multiply(a1,a2,temp)
          cala_volue=dabs(sum( a3*temp   ))
  end function

  subroutine creatk_mesh_01(self)
      implicit none
      type(Kspace),intent(inout)::self
      !--------------------------------
      real*8,parameter::pi = 3.141592653589793238462643383279_8
      real*8::vol,tempV(3),tempV2(3),mid(3),G(3,3),J,lenth
      real*8,allocatable::x(:,:),w(:,:)
      type(Integration)::Intlist
      integer::jc,nmax,c(3),jcd
      !--------------------------------
      Vol=cala_volue(self%a(:,1),self%a(:,2),self%a(:,3))
      ! calculate reciprocal vector G1,G2 and G3
      call vector3D_multiply(self%a(:,2),self%a(:,3),G(:,1)) ;self%b(:,1)=G(:,1)*Pi*2._8/Vol
      call vector3D_multiply(self%a(:,3),self%a(:,1),G(:,2)) ;self%b(:,2)=G(:,2)*Pi*2._8/Vol
      call vector3D_multiply(self%a(:,1),self%a(:,2),G(:,3)) ;self%b(:,3)=G(:,3)*Pi*2._8/Vol
      !-----------------------
      !--calculate  J
      tempV = self%b(:,1)/dsqrt(sum(self%b(:,1)**2._8))
      tempV2= self%b(:,2)/dsqrt(sum(self%b(:,2)**2._8))
      call vector3D_multiply(tempV,tempV2,mid)
      J = sum (  mid *   self%b(:,3)/dsqrt(sum(self%b(:,3)**2._8))  )
      !-----------------------
      nmax=max(self%n(1),self%n(2),self%n(3))
      allocate(x(nmax,3),w(nmax,3))
      do jc=1,3
         call Intlist%Initialization(self%n(jc),  self%meshtype )
         call Intlist%get_xw(-1._8,1._8,x(1:self%n(jc),jc),w(1:self%n(jc),jc))
      end do
      !-------------------
      c(1)=0_8;c(2)=1_8;c(3)=1_8
      DO JC=1_8,self%nk
       !----------------------
          c(1)=c(1)+1_8
          if (c(1).gt.self%n(1))then
            c(2)=c(2)+1_8
            if (c(2).gt.self%n(2))then
               c(3)=c(3)+1_8
               c(2)=1_8
            end if
            c(1)=1_8
          end if
       !---------(c1,c2,c3)----------
       self%w(  jc)=1._8
       self%k(:,jc)=0._8
       do jcd=1,3
           lenth=dsqrt(sum(self%b(:,JCd)**2))
           self%k(:,jc)=self%k(:,jc)+x(c(jcd),jcd)*self%b(:,JCd)/2._8
           self%w(  jc)=self%w( jc)*w(c(jcd),jcd)*lenth   /2._8
       end do
       self%w(  jc) = self%w(  jc) * J
      END DO
      self%w=self%w/(2._8*pi)**3*Vol
      deallocate(x,w)
  end subroutine


  subroutine CheckCanBeGetOrStop(self,i)
    implicit none
    class(Kspace),intent(inout)::self
    integer,intent(in)::i
    !--------------------------------
    call self%CheckInitiatedOrStop()
    if (  (i.ge.1) .and. (i.le.self%nk)  )then
    else
      write(self%getprint(),*)"ERROR: illegal index i=",i
      write(self%getprint(),*)"i should >0 and < nk=",self%nk
      stop
    endif
  endsubroutine


  function getk(self,i) result(r)
    implicit none
    class(Kspace),intent(inout)::self
    integer,intent(in)::i
    real*8::r(3)
    !--------------------------------
    call CheckCanBeGetOrStop(self,i)
    r = self%k(:,i)
  endfunction

  real*8 function getw(self,i)
    implicit none
    class(Kspace),intent(inout)::self
    integer,intent(in)::i
    !--------------------------------
    call CheckCanBeGetOrStop(self,i)
    getw = self%w(i)
  endfunction

  integer function getnk(self)
    implicit none
    class(Kspace),intent(inout)::self
    !--------------------------------
    call self%CheckInitiatedOrStop()
    getnk = self%nk
  endfunction

  subroutine CheckDimentionGetAllowed(self,i)
    implicit none
    class(Kspace),intent(inout)::self
    integer,intent(in)::i
    !--------------------------------
    call self%CheckInitiatedOrStop()
    if (  (i.ge.1) .and. (i.le.3)  )then
    else
      write(self%getprint(),*)"ERROR: dimension index =",i," is out of range [1,3]"
      stop
    endif
  endsubroutine

  function geta(self,i) result(r)
    implicit none
    class(Kspace),intent(inout)::self
    integer,intent(in)::i
    real*8::r(3)
    !--------------------------------
    call CheckDimentionGetAllowed(self,i)
    r = self%a(:,i)
  endfunction

  function GetB(self,i) result(r)
    implicit none
    class(Kspace),intent(inout)::self
    integer,intent(in)::i
    real*8::r(3)
    !--------------------------------
    call CheckDimentionGetAllowed(self,i)
    r = self%b(:,i)
  endfunction



  function getvolr(self) result(r)
    implicit none
    class(Kspace),intent(inout)::self
    real*8::r
    !--------------------------------
    r =  cala_volue(self%a(:,1),self%a(:,2),self%a(:,3))
  endfunction

  function getvolb(self) result(r)
    implicit none
    class(Kspace),intent(inout)::self
    real*8::r
    !--------------------------------
    r =  cala_volue(self%b(:,1),self%b(:,2),self%b(:,3))
  endfunction



endmodule
