
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : ED_phyCal
! OBJECT: TYPE(EDPhyCal)
! USED  : FermionOperators,FermionHamiltonian,ED_WholeSpace
! DATE  : 2017-11-28
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            consider a subspace in ED. This subspace is identified by a subid which is correpsonding to that in table.
!
! STANDARD:
!            *CALL Initialization(T,ClustH,subid)
!             call diagonalization()
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  [sub] Initialization(ED,*PRINT)
!                        class(ED_GCE),pointer::ED
!
! avalable gets:
!
!                  [fun] GetTraceValue(Intertype,para)
!                        complex*16::GetTraceValue
!                        Intertype = character(16) / integer
!                        integer::para(8)
!                        if   Intertype = character(16), the usage is as the Interaction type in module FermionHamiltonian
!                        if   Intertype = integer  , the usage is input a single operator. See module FermionOperators
!
! avalable is :
!                  [fun] i
! others      :
!                  [sub] printE()
!                        print out eigenenergy
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


module ED_phyCal
  use ED_WholeSpace
   implicit none


   type::EDPhyCal
     private
     logical::Initiated        = .false.
     class(ED_GCE),pointer::ED => null()


     integer::print = 6
   contains
     procedure,pass::Initialization
     final::Finalization

     procedure,pass::GetTraceValue
   endtype

   private::Initialization,UnInitialization,Finalization

   private::GetTraceValue

 contains


   subroutine Initialization(self,ED,print_)
     implicit none
     class(EDPhyCal),intent(inout)::self
     class(ED_GCE),intent(inout),target::ed
     integer,intent(in),optional::print_
     !-----------------------------------
     call UnInitialization(self)
     self%initiated = .true.

     self%ed => ed
     if (present(print_)) self%print = print_
   endsubroutine


   subroutine UnInitialization(self)
     implicit none
     class(EDPhyCal),intent(inout)::self
     !-----------------------------------
     if (self%initiated) then
       self%initiated = .false.
       self%ed => null()
     endif
   endsubroutine


   subroutine Finalization(self)
     implicit none
     TYPE(EDPhyCal),intent(inout)::self
     !-----------------------------------
     call UnInitialization(self)
   endsubroutine



   complex*16 function GetTraceValue(self,Intertype,para)
     implicit none
     class(EDPhyCal),intent(inout)::self
     class(*),intent(inout)::Intertype
     integer,intent(in)::para(8)
     !-----------------------------------
     class(FermOper),pointer::f
     TYPE(Ham)::Htemp
     integer::jc

     select type(Intertype)
     Type is (Integer)
       allocate(f)
       call f%Initialization(self%ed%get_ns(),Intertype,para,self%print)
       GetTraceValue = self%ed%get_trace_value(f)
       deallocate(f)
     Type is (character(*))
       call Htemp%Initialization( self%ed%get_ns() , print_=self%print )
       call Htemp%StartAppendingInteraction()
       call Htemp%AppendingInteraction(Intertype,para,(1._8,0._8))
       call Htemp%EndAppendingInteraction()
       GetTraceValue = (0._8,0._8)
       do jc = 1 , Htemp%GetOptN()
          f => Htemp%GetOptact(jc)
          GetTraceValue = GetTraceValue + self%ed%get_trace_value(f)
       enddo
     endselect
   endfunction







endmodule
