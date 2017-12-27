!
!
!
!
! !!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! ! TYPE  : MODULE
! ! NAME  : ThreeDSpaceCal
! ! OBJECT: TYPE(TdCal)
! ! USED  :
! ! DATE  : 2017-12-22
! ! AUTHOR: hengyueli@gmail.com
! !--------------
! ! Open-Source : No
! !------------------
! ! DESCRIPTION:
! !            Offer some methods for the calculatins in 3D space.
! !
! ! STANDARD:
! !
! !
! ! USING LIST:
! !            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! !
! !
! ! avalable sets:
! !                   [sub] I
! !
! ! avalable gets:
! !                   [fun] GetProduct(a,b)
! !                         return GetProduct = <a|b>
! !
! !                   [fun] GetCross(a,b)
! !                         return GetProduct(3) = a X b
! !
! !                   [fun] GetTripleproduct(v1,v2,v3)
! !                         return real*8 = v1.( v2 X v3 ) = ( v1 X v2 ) . v3
! !
! !                   [fun] GetVofTriclinic(V1,V2,V3)
! !                         for three given vectors v1,v2 and v3, we can unitly fix a  Triclinic shape, get the volume
! !
! !                   [fun] GetVofTetrahedron(v1,v2,v3)
! !                         for three given vectors v1,v2 and v3, we can unitly fix a  Tetrahedron, get the volume
! !
! ! avalable is :
! !                  ![fun] i
! ! others      :
! !                  ![sub] p
! !
! !
! !
! !!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!
!
! module ThreeDSpaceCal
!   implicit none
!
!
!   type::TdCal
!     private
!
!
!
!   contains
!     procedure,nopass::GetProduct
!     procedure,nopass::GetCross
!     procedure,nopass::GetVofTriclinic
!     procedure,nopass::GetTripleproduct
!     procedure,nopass::GetVofTetrahedron
!   endtype
!
!   private::GetProduct,GetCross,GetTripleproduct
!   private::GetVofTetrahedron,GetVofTriclinic
!
! contains
!
!
!   Function GetProduct(a,b)     result(c)
!     implicit none
!     real*8,intent(in)::a(3),b(3)
!     real*8::c
!     !--------------------------------
!     c = sum( a * b )
!   endfunction
!
!   Function GetCross(a,b)     result(c)
!     implicit none
!     real*8,intent(in)::a(3),b(3)
!     real*8::c(3)
!     !--------------------------------
!     c(1) = a(2) * b(3) - a(3) * b(2)
!     c(2) = a(3) * b(1) - a(1) * b(3)
!     c(3) = a(1) * b(2) - a(2) * b(1)
!   endfunction
!
!   real*8 function GetTripleproduct(v1,v2,v3)
!     implicit none
!     real*8,intent(in)::v1(3),v2(3),v3(3)
!     !------------------------------------------
!     GetTripleproduct = &
!                        v1(1) * v2(2) * v3(3) - v1(1) * v2(3) * v3(2) + &
!                        v1(2) * v2(3) * v3(1) - v1(2) * v2(1) * v3(3) + &
!                        v1(3) * v2(1) * v3(2) - v1(3) * v2(2) * v3(1)
!   endfunction
!
!
!   real*8 function GetVofTriclinic(v1,v2,v3)
!     implicit none
!     real*8,intent(in)::v1(3),v2(3),v3(3)
!     !------------------------------------------
!     GetVofTriclinic = dabs(  GetTripleproduct(v1,v2,v3)  )
!   endfunction
!
!
!   real*8 function GetVofTetrahedron(V1,V2,V3)
!     implicit none
!     real*8,intent(in)::v1(3),v2(3),v3(3)
!     !------------------------------------------
!     GetVofTetrahedron = GetVofTriclinic(v1,v2,v3) / 2._8
!   endfunction
!
!
!
! endmodule
