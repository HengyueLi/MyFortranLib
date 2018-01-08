






!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : phyfunc
! OBJECT: TYPE(phyfun)
! USED  :
! DATE  : 2017-12-01
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            offer some frequently used functions in physics.
!
! STANDARD:
!           [sub] I
!
!
! USING LIST:
!
!
!
! avalable sets:
!                  [sub] I
! avalable gets:
!
!                  [fun] FermionFunc(Omega,T)
!                        real*8           ::T  ! temperature
!                        real*8/complex*16:: Omega
!                        get the fermin function : f = 1/( exp(omega/T) + 1 )
!
!                  [fun] GetPauliMatrix(direction)
!                        integer::direction     ! = 0,1,2,3  for  I,x,y,z respectivily
!                        return a 2 by 2 matrix complex*16::GetPauliMatrix(2,2)
! avalable IS  :
!                  [fun] I
!
! avalable others:
!                  [sub] s
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


module phyfunc
       implicit none



       type::phyfun


       contains
            procedure,nopass::FermionFuncR
            procedure,nopass::FermionFuncC
            generic::FermionFunc=>FermionFuncR,FermionFuncC

            procedure,nopass::GetPauliMatrix

       endtype



       private::FermionFuncR,FermionFuncC,GetPauliMatrix


     contains

     real*8 function FermionFuncR(Omega,T)
             implicit none
             real*8,intent(in)::Omega,T
             !--------------------------------------
             real*8::temp
             temp = omega/T
             if (omega.ge.0)then
               FermionFuncR = dexp(-temp) / ( dexp(-temp) + 1._8  )
             else
               FermionFuncR = 1._8/(dexp(temp)+1._8)
             endif
           endfunction

   complex*16 function FermionFuncC(Omega,T)
           implicit none
           complex*16,intent(in)::Omega
           real*8,intent(in)::T
           !--------------------------------------
           complex*16::temp
           temp = omega/T
           if (real(temp).ge.0)then
             FermionFuncC = zexp(-temp) / ( zexp(-temp) + 1._8  )
           else
             FermionFuncC = 1._8/(zexp(temp)+1._8)
           endif
         endfunction


  function GetPauliMatrix(Direction) result(r)
    implicit none
    integer,intent(in)::Direction
    complex*16::r(2,2)
    !---------------------------------------
    select case(Direction)
    case(0)
      r(1,1) = ( 1._8 , 0._8 )  ;   r(1,2) = ( 0._8 , 0._8 )
      r(2,1) = ( 0._8 , 0._8 )  ;   r(2,2) = ( 1._8 , 0._8 )
    case(1)
      r(1,1) = ( 0._8 , 0._8 )  ;   r(1,2) = ( 1._8 , 0._8 )
      r(2,1) = ( 1._8 , 0._8 )  ;   r(2,2) = ( 0._8 , 0._8 )
    case(2)
      r(1,1) = ( 0._8 , 0._8 )  ;   r(1,2) = ( 0._8 ,-1._8 )
      r(2,1) = ( 0._8 , 1._8 )  ;   r(2,2) = ( 0._8 , 0._8 )
    case(3)
      r(1,1) = ( 1._8 , 0._8 )  ;   r(1,2) = ( 0._8 , 0._8 )
      r(2,1) = ( 0._8 , 0._8 )  ;   r(2,2) = (-1._8 , 0._8 )
    case default
      write(*,*)"Inccorect input for getting Pauli Matrix, input can only be [0,3]"
      write(*,*)"Now it is:",Direction
      stop
    endselect
  endfunction





endmodule
