


!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : VCA_WaldFun
! OBJECT: TYPE(waldf)
! USED  : CodeObject
! DATE  : 2018-01-02
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!
!             Evaluate the Wald functional of the system. That is :
!                               \Omega = \Omega' + Trln[-G] - Trln[-G'] = \Omega' - I
!             where  I = T * \sum_{omega_n}\sum_k lndet[1-Vkd * G']
!
!             basically, we seperate the Integration into two part: I = I_i + I_s
!             where
!                  I_i is an integration part and I_s is an summation part.
!             The different of I_i and I_s can be concludes as the only difference of integration weight.
!
! STANDARD:
!            *CALL Initialization(  )
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                   [sub] I
!
! avalable gets:
!
!                   [fun] G
!
!
! avalable is :
!                  ![fun] i
! others      :
!                   [sub] A
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


module VCA_WaldFun
  use CodeObject
  implicit none


  integer,parameter::NparaMax = 50

  type,extends(object)::waldf
    private
    integer::method
    integer::jobi(NparaMax)
    real*8 ::jobr(NparaMax)



  endtype


  private::Initialization


contains

  subroutine Initialization(self,method,jobi,jobr,print_,show_)
    implicit none
    class(waldf),intent(inout) :: self
    integer,intent(in)         :: method
    integer,intent(in)         :: jobi(:)
    real*8 ,intent(in)         :: jobr(:)
    integer,intent(in),optional:: print_
    integer,intent(in),optional:: show_
    !-------------------------------------------------------------



  endsubroutine



endmodule
