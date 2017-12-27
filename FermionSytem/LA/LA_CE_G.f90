
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  :  MODULE
! NAME  :  LA_CE_G
! OBJECT:  TYPE(LACEG)
! USED  :  LA_Subspace,LA_CE_G_PQ
! DATE  :  2017-12-05
! AUTHOR:  hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            Calculate the single particle Green's function of Canonical Emsenbel
!
! STANDARD:
!            *CALL I
!
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  [sub] Initialization(LaSub,Nopt ,M_,OTH_,BZERO_,PRINT_,SHOW_)
!                        class(LASubSpace)::LaSub
!                        integer::Nopt(:,:)     Nopt(1,:) = site   Nopt(2,:) = spin
!                        integer::M_    = 90
!                        logical::OTH_  = true
!                        real*8::bzero_ = 1.e-7
!                        integer::PRINT_ = 6
!
!                  [sub] SynchronizeQP()
!                        check and Synchronize P&Q matrix
!---------------------------------------------------
! avalable gets:
!                   [sub] GetGmatrix(Nomega,Omega,G)
!                         integer::Nomega
!                         complex*16::Omega(Nomega),G(Nomega,ni,ni)
!                         NOTICE!!!:  The index of matrix G is not site index.
!                                     The are the same as the input Nopt.
!                   ![fun] GetSubDe()
!                         get the degeneracy of subspace.
! avalable is :
!                  ![fun] i
! others      :
!                  [sub] S
!
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




module LA_CE_G
  use LA_CE_G_PQ
  implicit none

   type::LACEG
     private
     logical::Initiated = .false.
     integer::M         = 90
     logical::oth       = .true.
     real*8 ::bzero     = 1.e-7!1.e-7
     TYPE(LAPQ)::P
     TYPE(LAPQ)::Q
     !--------------
     integer::Ni
     integer,allocatable::Nopt(:,:)
     !--------------
     integer::print
     INTEGER::SHOW = 0

   contains
     procedure,pass::Initialization
     final::Finalization


     procedure,pass::SynchronizeQP
     procedure,pass::GetGmatrix

   endtype



   private::Initialization,UnInitialization,Finalization

   private::SynchronizeQP,GetGmatrix

 contains


   subroutine Initialization(self,LaSub,Nopt,M_,OTH_,BZERO_,PRINT_,SHOW_)
     use basic_math_functions
     implicit none
     class(LACEG),intent(inout)             :: self
     class(LASubSpace),intent(inout),target :: LaSub
     integer,intent(in)                     :: Nopt(:,:)
     integer,intent(in),optional            :: M_
     logical,intent(in),optional            :: OTH_
     real*8,intent(in),optional             :: bzero_
     integer,intent(in),optional            :: PRINT_,SHOW_
     !------------------------------------------------------
     integer::Ni,spin,jc,para(8),ns
     type(FermOper),allocatable::fp(:),fq(:)
     class(table),pointer::Ta
     integer::sym,psubid,qsubid,subid,markp(2),markq(2),marker(2)
     TYPE(bmathf)::f

     call UnInitialization(self)  ;  self%initiated = .true.

     if(present(M_    )) self%m     = M_
     if(present(OTH_  )) self%oth   = oth_
     if(present(bzero_)) self%bzero = bzero_
     if(present(print_)) self%print = print_
     if(present(SHOW_ )) self%SHOW  = SHOW_

     Ni = size(Nopt(1,:))            ; self%ni = ni    ;allocate(self%Nopt,source=Nopt)

    !  write(*,*)nopt;stop


     Ta => LaSub%GetTablePoinnter()
     ns = ta%get_ns()
     allocate(  fp(Ni) ,fq(Ni) )

     do jc = 1 , Ni
       para    = 0
       para(1) = self%Nopt(1,jc)
       para(2) = self%Nopt(2,jc)
       call fp(jc)%Initialization(ns,2,para,self%print)
       call fq(jc)%Initialization(ns,1,para,self%print)
     enddo


     !--------------select subid------------------------------------------------
     sym = ta%get_symmetry()
     subid = LaSub%GetSubId()
     call ta%get_sub_mark_value(subid,marker)
     select case(sym)
     case(0)
       psubid = 1
       qsubid = 1
     case(1)
       markp(1) = marker(1) + 1 ; markp(2) = -1
       markq(1) = marker(1) - 1 ; markq(1) = -1
       psubid   = ta%get_subid_from_mark(markp)
       qsubid   = ta%get_subid_from_mark(markq)
     case(2)
       if ( f%IsIntegerArrayTheSame(Nopt(2,:)) )then
         spin =  Nopt(2,1)
         markp(spin+1) = marker(spin+1) + 1 ; markp(2-spin) = marker(2-spin)
         markq(spin+1) = marker(spin+1) - 1 ; markq(2-spin) = marker(2-spin)
         psubid   = ta%get_subid_from_mark(markp)    !;write(*,*)ta%get_subspace_marker(14);stop
         qsubid   = ta%get_subid_from_mark(markq)
       else
         write(self%print,*)"The input spin are different while symmetry = 2";stop
       endif
     endselect
     !--------------------------------------------------------------------------

     call self%P%Initialization(LaSub=LaSub,Ni=ni,OPT=fp,PQSubid=psubid,Asign=1,&
           M_=self%m,OTH_=self%oth,BZERO_=self%bzero,PRINT_=self%print,SHOW_=SELF%SHOW)
     call self%Q%Initialization(LaSub=LaSub,Ni=ni,OPT=fq,PQSubid=qsubid,Asign=-1,&
           M_=self%m,OTH_=self%oth,BZERO_=self%bzero,PRINT_=self%print,SHOW_=SELF%SHOW)

   endsubroutine

   subroutine UnInitialization(self)
     implicit none
     class(LACEG),intent(inout)          :: self
     !------------------------------------------------------
     if (self%initiated)then ;  self%initiated = .false.
        if (allocated(self%Nopt))  deallocate(self%Nopt)
     endif
   endsubroutine

   impure elemental subroutine Finalization(self)
     implicit none
     type(LACEG),intent(inout)          :: self
     !------------------------------------------------------
     call UnInitialization(self)
   endsubroutine


   subroutine SynchronizeQP(self)
     implicit none
     class(LACEG),intent(inout)          :: self
     !------------------------------------------------------
     call self%p%SynchronizeQP()
     call self%q%SynchronizeQP()
   endsubroutine


   subroutine GetGmatrix(self,Nomega,Omega,G)
     implicit none
     class(LACEG),intent(inout)          :: self
     integer,intent(in)::Nomega
     complex*16,intent(in)::Omega(Nomega)
     complex*16,intent(out)::G(Nomega,self%ni,self%Ni)
     !------------------------------------------------------
     complex*16::Gt(Nomega,self%ni,self%Ni)
     call SynchronizeQP(self)

     call self%p%GetGpq(Tran=.false.,Nomega=Nomega,Omega=Omega,G=G )
     call self%q%GetGpq(Tran=.true. ,Nomega=Nomega,Omega=Omega,G=Gt)

     G = g + gt
   endsubroutine
!GetGpq(Tran,Nomega,Omega,G)





endmodule
