!
!
! !!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! ! TYPE  : MODULE
! ! NAME  : FermionActions
! ! OBJECT: TYPE(FerAct)
! ! USED  : fermion_table , FermionOperators
! ! DATE  : 2017-12-10
! ! AUTHOR: hengyueli@gmail.com
! !--------------
! ! Open-Source : No
! !------------------
! ! DESCRIPTION:
! !            For any fermion state, collect all actions
! !
! ! STANDARD:
! !            *CALL Initialization(Ta,IsReal,PRINT_,show_)
! !
! ! USING LIST:
! !            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! !
! !
! ! avalable sets:
! !                  [sub] Initialization(Ta,IsReal,PRINT_,show_)
! !
! ! avalable gets:
! !                   [sub] S
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
! module FermionActions
!   use fermion_table
!   use FermionOperators
!   implicit none
!
!   type::FerAct
!     private
!     logical::Initiated = .false.
!
!
!     integer::print = 6
!     integer::show  = 0
!   endtype
!
!
!   ! private::Initialization,UnInitialization,Finalization
!
!
! contains
!
!   ! subroutine
!
!
! endmodule
