


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : SpelFunc
! OBJECT: TYPE(funcs)
! USED  : CodeObject
! DATE  : 2018-01-02
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!               Special Funtionals
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
!                 1.Legendre Polynomials :
!                       p(n,x)
!                       p(0,x) = 1
!                       p(1,x) = x
!                       k * p(k,x) = (2k-1) * x * p(k-1,x) - (k-1) * p(k-2,x) (k>=1)
!
!                 --------------
!                  [fun] GetLegendreP(n,x)
!                        integer::n
!                        real*8 ::x
!                  [fun] GetLegendreDxP(n,x)
!                        derivative of GetLegendreP
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



module SpelFunc
  use CodeObject
  implicit none

  type,extends(object)::funcs
    !
  contains
    procedure,pass::GetLegendreP
    procedure,pass::GetLegendreDxP
  endtype


  private::GetLegendreP,GetLegendreDxP



contains


  real*8 function GetLegendreP(self,n,x)
    implicit none
    class(funcs),intent(in)::self
    integer,intent(in)::n
    real*8,intent(in)::x
    !-------------------------------
    integer::jc
    REAL*8::pk_1,pk_2
    if (n.lt.0)then
      write(self%getprint(),*)"ERROR: n<0 is not defined in GetLegendreP"
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
    class(funcs),intent(in)::self
    integer,intent(in)::n
    real*8,intent(in)::x
    !-------------------------------
    if (n.lt.0)then
      write(self%getprint(),*)"ERROR: n<0 is not defined in GetLegendreP"
    elseif(n==0)then
      GetLegendreDxP = 0._8
    else
      GetLegendreDxP = n * ( x *GetLegendreP(self,n,x) - GetLegendreP(self,n-1,x) )/(x**2.0_8-1._8)
    endif 
  endfunction






endmodule
