

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : ED_GCE_G
! OBJECT: TYPE(EDGCEGreenf)
! USED  : ED_WholeSpace
! DATE  : 2017-12-01
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            consider a subspace in ED. This subspace is identified by a subid which is correpsonding to that in table.
!
! STANDARD:
!            *CALL Initialization(ED)
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  [sub] Initialization(ED,*PRINT)
!
!
! avalable gets:
!
!                  [sub] GetG(sitei,spini,sitej,spinj,Nomega,Omega,G)
!                        integer::  sitei,spini,sitej,spinj,Nomega
!                        complex*16::Omega(Nomega),G(Nomega)
!                        for an input Omega, get the dynamic Green's function G
!                        In this subroutine, a Synchronize of ED with H will be done.
!
!
! avalable is :
!                  [fun] i
! others      :
!                  [sub] p
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module ED_GCE_G
  use ED_WholeSpace
   implicit  none




     type:: EDGCEGreenf
          private

          logical::Initiated        = .false.
          integer::EigenId          = -1
          class(ED_GCE),pointer::ED =>null()
          class(table),pointer ::Ta => null()

          !----------------
          integer::print = 6


        contains
          procedure,pass::Initialization
          final::Finalization

          procedure,pass::GetG
     endtype



     private::Initialization,UnInitialization,Finalization

     private::GetG,GetG_sy0,GetG_sy1,GetG_sy2,GetG_Submn_Component

   contains


     subroutine Initialization(self,ED,print_)
       implicit none
       class(EDGCEGreenf),intent(inout)::self
       class(ED_GCE),target::ED
       integer,intent(in),optional::print_
       !------------------------------------------
       call UnInitialization(self)
       self%initiated = .true.
       self%ED => ED
       self%Ta => ED%get_Table_pointer()
       if (present(print_)) self%print = print_
     endsubroutine


     subroutine UnInitialization(self)
       implicit none
       class(EDGCEGreenf),intent(inout)::self
       !------------------------------------------
       if (self%initiated)then
          self%initiated = .false.
          self%ED => null()
          self%Ta => null()
       endif
     endsubroutine

     subroutine Finalization(self)
       implicit none
       type(EDGCEGreenf),intent(inout)::self
       !------------------------------------------
       call UnInitialization(self)
     endsubroutine



    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !                         core sub
     subroutine GetG(self,sitei,spini,sitej,spinj,Nomega,Omega,G)
       use FermionOperators
       implicit none
       class(EDGCEGreenf),intent(inout)::self
       integer,intent(in)::sitei,spini,sitej,spinj,Nomega
       complex*16,intent(in)::Omega(Nomega)
       complex*16,intent(out)::G(Nomega)
       !------------------------------------------
       TYPE(FermOper)::ci,cj
       integer::para(8)
       Logical::Diagonal
       Diagonal = ( sitei .eq.sitej   )  .and.  (spini.eq.spinj)
       para(1) = sitei ; para(2) = spini ; call ci%Initialization(self%Ta%get_ns(),1,para)
       para(1) = sitej ; para(2) = spinj ; call cj%Initialization(self%Ta%get_ns(),1,para)

       call self%ed%SynchronizeWithHamiltonian()
       select case(self%ta%get_symmetry())
       case(0)
         call GetG_sy0(self,ci,cj,Diagonal,Nomega,Omega,G)
       case(1)
         call GetG_sy1(self,ci,cj,Diagonal,Nomega,Omega,G)
       case(2)
         if (spini.ne.spinj)then
           G = (0._8,0._8)
         else
           call GetG_sy2(self,ci,cj,Diagonal,Nomega,Omega,G)
         endif
       case default
         write(self%print,*)"ERROR: symmetry =",self%ta%get_symmetry(),"is unknow value in GetG"
         write(self%print,*)"This is very strange.";stop
       endselect
     endsubroutine

     subroutine GetG_sy0(self,ci,cj,Diagonal,Nomega,Omega,G)
       use FermionOperators
       implicit none
       class(EDGCEGreenf),intent(inout)::self
       class(FermOper),intent(inout)::ci,cj
       logical,intent(in)::Diagonal
       integer,intent(in)::Nomega
       complex*16,intent(in)::Omega(Nomega)
       complex*16,intent(out)::G(Nomega)
       !------------------------------------------
       call GetG_Submn_Component(self,Ci,Cj,Diagonal,1,1,Nomega,Omega,G)
     endsubroutine

     subroutine GetG_sy1(self,ci,cj,Diagonal,Nomega,Omega,G)
       use FermionOperators
       implicit none
       class(EDGCEGreenf),intent(inout)::self
       class(FermOper),intent(inout)::ci,cj
       logical,intent(in)::Diagonal
       integer,intent(in)::Nomega
       complex*16,intent(in)::Omega(Nomega)
       complex*16,intent(out)::G(Nomega)
       !------------------------------------------
       complex*16::DeltaG(Nomega)
       Integer::subidm,subidn
       integer::nsm,markerm(2),markern(2)

       G = (0._8,0._8)
       do nsm = 1 , 2*self%ta%get_ns()
          markerm(1) = nsm     ; markerm(2) = -1
          markern(1) = nsm - 1 ; markern(2) = -1
          subidm = self%ta%get_subid_from_mark(markerm)
          subidn = self%ta%get_subid_from_mark(markern)
          call GetG_Submn_Component(self,Ci,Cj,Diagonal,Subidn,Subidm,Nomega,Omega,DeltaG)
          G = G + DeltaG
       enddo
     endsubroutine

     subroutine GetG_sy2(self,ci,cj,Diagonal,Nomega,Omega,G)
       use FermionOperators
       implicit none
       class(EDGCEGreenf),intent(inout)::self
       class(FermOper),intent(inout)::ci,cj
       logical,intent(in)::Diagonal
       integer,intent(in)::Nomega
       complex*16,intent(in)::Omega(Nomega)
       complex*16,intent(out)::G(Nomega)
       !------------------------------------------
       complex*16::DeltaG(Nomega)
       Integer::subidm,subidn,para(8)
       integer::nsm,sigma,sigmabar,NspinB,markerm(2),markern(2) ,ns,spin

       para = ci%get_para()  ; spin = para(2)
       G = (0._8,0._8)
       ns = self%ta%get_ns()
       sigma = spin  ;   sigmabar = 1 - spin
       do nsm = 1 , ns  ; Do NspinB = 0 , ns
          markerm(sigma+1) = nsm     ; markerm(sigmabar+1) = NspinB
          markern(sigma+1) = nsm - 1 ; markern(sigmabar+1) = NspinB
          subidm = self%ta%get_subid_from_mark(markerm)
          subidn = self%ta%get_subid_from_mark(markern)
          call GetG_Submn_Component(self,Ci,Cj,Diagonal,Subidn,Subidm,Nomega,Omega,DeltaG)
          G = G + DeltaG
       enddo      ;enddo
999    continue
     endsubroutine




    ! The conponent of G would comes from    <n|ci|m> x <n|cj|m>^*
    !                                     where n \in subn and m \in subm
    !                                     subm and subn can be the same(Symmetry=0).
     subroutine GetG_Submn_Component(self,OptCi,OptCj,IsDiagonal,Subidn,Subidm,Nomega,Omega,DeltaG)
       use FermionOperators
       implicit none
       class(EDGCEGreenf),intent(inout)::self
       class(FermOper),intent(inout)::OptCi,OptCj
       logical,intent(in)::IsDiagonal
       integer,intent(in)::Subidn,Subidm,Nomega
       complex*16,intent(in)::Omega(Nomega)
       complex*16,intent(out)::DeltaG(Nomega)
       !------------------------------------------
       integer::dm,dn,jm,jn
       complex*16,allocatable::si(:),sj(:)
       complex*16::upperi,upperj,upper
       Real*8::Em,En,Z,T,ExpMTm,Emn,ExpMTmn
       logical::Diagonal

       DeltaG   = (0._8,0._8)
       Z        = self%ed%GetRz()
       T        = self%ed%get_T()
       Diagonal = IsDiagonal
       dm       = self%ta%get_sub_d(subidm)
       dn       = self%ta%get_sub_d(subidn)
       allocate(si(dn),sj(dn))
       do jm = 1 , dm
         Em     = self%ed%get_E(subidm,jm)
         ExpMTm = dexp( -Em/T  )
         do jn = 1, dn
           En = self%ed%get_E(subidn,jn)
           ExpMTmn = dexp( -En/T  ) + ExpMTm
           Emn    = Em -En
           !---------------------------------
           si = self%ed%get_operator_act_on_state(Optci,subidm,jm,dn)
           upperi = self%ed%get_State_product(subidn,jn,dn,si)
           if (Diagonal)then
              upperj = upperi
           else
              sj = self%ed%get_operator_act_on_state(Optcj,subidm,jm,dn)
              upperj = self%ed%get_State_product(subidn,jn,dn,sj)
           endif
           upper  = upperi * conjg(upperj) * ExpMTmn
           DeltaG = DeltaG + upper / ( Omega - Emn)
         enddo
       enddo
       deallocate(si,sj)
       DeltaG = DeltaG / Z
     endsubroutine
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






endmodule
