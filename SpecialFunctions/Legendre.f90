


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : LegendreFun
! OBJECT: TYPE(func)
! USED  : CodeObject
! DATE  : 2018-01-02
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!               Special Funtionals-  Legendre Polynomials :
!
!                       p(n,x)
!                       p(0,x) = 1
!                       p(1,x) = x
!                       k * p(k,x) = (2k-1) * x * p(k-1,x) - (k-1) * p(k-2,x) (k>=1)
!
! STANDARD:
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  []
!
!
! avalable gets:
!

!
!                 --------------
!                  [fun] GetLegendreP(n,x)
!                        integer::n
!                        real*8 ::x
!
!                  [fun] GetLegendreDxP(n,x)
!                        derivative of GetLegendreP
!
!                  [fun] GetZerosOfLegendre(n,i)
!                        for P(n,x), it contains n zeros. Return the i-th zero (ascending order)
! avalable Is:
!
! others:
!
!
!
!
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module LegendreFun
  use CodeObject
  implicit none

  type,extends(object)::func
    !
  contains
    procedure,pass::GetLegendreP
    procedure,pass::GetLegendreDxP
    procedure,pass::GetZerosOfLegendre
  endtype


  private::GetLegendreP,GetLegendreDxP



contains


  real*8 function GetLegendreP(self,n,x)
    implicit none
    class(func),intent(in)::self
    integer,intent(in)::n
    real*8,intent(in)::x
    !-------------------------------
    integer::jc
    REAL*8::pk_1,pk_2
    if (n.lt.0)then
      write(self%getprint(),*)"MathError: n<0 is not defined in GetLegendreP";stop
    elseif(n==0)then
      GetLegendreP = 1._8
    elseif(n==1)then
      GetLegendreP = x
    else
      pk_1 = x
      pk_2 = 1._8
      do jc = 2 , n
         GetLegendreP = (  (2*jc-1) * x * pk_1 - (jc-1) * pk_2 ) / jc
         pk_2 = pk_1
         pk_1 = GetLegendreP
      enddo
    endif
  endfunction

  real*8 function GetLegendreDxP(self,n,x)
    implicit none
    class(func),intent(in)::self
    integer,intent(in)::n
    real*8,intent(in)::x
    !-------------------------------
    if (n.lt.0)then
      write(self%getprint(),*)"MathError: n<0 is not defined in GetLegendreP.";stop
    elseif(n==0)then
      GetLegendreDxP = 0._8
    else
      GetLegendreDxP = n * ( x *GetLegendreP(self,n,x) - GetLegendreP(self,n-1,x) )/(x**2.0_8-1._8)
    endif
  endfunction


  function GetZerosOfLegendre(self,n,i) result(r)
    implicit none
    class(func),intent(in)::self
    integer,intent(in)::n,i
    real*8::r
    !-------------------------------
    real*8,parameter::pi=3.14159265358979323846264338327950288439_8
    real*8::Pre,temp1

    Pre   = 0.000000000000001_8
    temp1 = 800.0_8  ! random large number

    if (n.lt.1)then
      write(self%getprint(),*)"MathError: For n<1 case, there is no zero for Legendre function.";stop
    else
      if ( (i.gt.0) .and. (i.le.n) )then
        !r = DCOS( (4*I-1)*pi/(4*N+2) )
        r = -DCOS( (4*I-1)*pi/(4*N+2) )   ! re-order
        do while ( abs( r - temp1) .gt.Pre )
           temp1 = r
           r     = r - GetLegendreP(self,n,temp1)/GetLegendreDxP(self,n,temp1)
        enddo
      else
        write(self%getprint(),*)"MathError: i should in range of [1,n] for GetZerosOfLegendre"
        write(self%getprint(),*)"However,now (i,n)=(",i,",",n,")";stop
      endif
    endif
  endfunction






endmodule
