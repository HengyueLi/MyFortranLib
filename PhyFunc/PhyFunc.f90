






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

       endtype



       private::FermionFuncR,FermionFuncC


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





endmodule
