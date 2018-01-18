

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  :  MODULE
! NAME  :  ED_WholeSpace
! OBJECT:  TYPE(ED_GCE)
! USED  :  ED_Subspace
! DATE  :  2017-11-28
! AUTHOR:  hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            All the subspaces in ED. This is a Grand Canonical Ensemble calculation
!
! STANDARD:
!            *CALL Initialization(T,Ta,CH)
!            *call diagonalization()
!            *call Normolize_Energy_And_GetZ
!
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  [sub] Initialization(T,Ta,CH,IsReal,print_)
!                        real(8)::T        ! temperature of the system
!                        type(table)::Ta    ! must be initiated first
!                        logical::Real
!                        type(Ham)::CH  !
!                  [sub] Normolize_Energy_And_GetZ()
!
! avalable gets:
!                  [fun] get_product(subi,i,subj,j,oper)
!                        integer(8):: subi,i,subj,j
!                        type(FermOper)::oper
!                        complex*16::get_product
!                        returnV =  <subi,i|oper|subj,j>
!
!                  [fun] get_State_product(subid,index,d,state)
!                        integer::subid,index,d
!                        complex*16::state(d)
!                        return <subid,index|state>
!                  [fun] get_operator_act_on_state(oper,subid,index,d)
!                        integer,intent(in)::subid,index,d
!                        complex*16,intent(out)::get_operator_act_on_state(d)
!                        for a given operator oper, output   |state>= oper |subid,index>
!                        state should be priviously allocated (as length of d).
!                  [fun] get_trace_value(opert)
!                        type(FermOper)::oper
!                        return complex*16
!                        returnV = 1/Z * \sum_{s} <s| oper |s> * Exp(-E_s/T)
!                        where s represent all the eigenstates in the Hilbert space and Z is the partition.
!                  [fun] get_E(subid,i)   ! get the normolized eigen energy.
!                        real(8)::get_E
!                        integer(8)::subid,i
!                  [fun] get_Vec(subid,i,d) ! get the eigen states
!                        integer(8)::subid
!                        integer(8)::i
!                        complex(8)::get_Vec(d) ! where d is the dimension of this subspace.
!                  [fun] get_Eg()
!                        get ground state energy
!                  [fun] get_Nsub()
!
!                 ! [fun] GetRz()
!                        The same as get_Z(Z), a function version
!                  [fun] get_T()
!                        real(8)::T  the temperature
!                  [fun] Get_Ns()
!                        return integer
!
!                  [sub] diagonalization()
!                        setting matrix elements is included.
!                  [fun] get_Table_pointer()
!                        class(table),pointer::get_Table_pointer
!                  [fun] get_EigenId()
!                        integer::get_EigenId
! avalable is :
!                  ![fun] i
! others      :
!                  [sub] SynchronizeWithHamiltonian()
!                        Check the EigenSytem of ED. If Hamiltonian is changed, renew EigenSytem
!                        Here renew means diagonalization +  Normolize_Energy_And_GetZ
!
!                  [sub] printE()
!                        print all eigenenergy
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module ED_WholeSpace
  use ED_Subspace
  implicit none

  TYPE::ED_GCE
    private
    logical                      :: initiated  = .false.
    logical                      :: Normolized = .false.
    logical                      :: IsReal     = .false.
    integer                      :: Nsub
    integer                      :: EigenId    = -1_8
    class(table),pointer         :: Ta         => null()
    class(Ham),pointer           :: CH         => null()
    type(EDSubSpace),ALLOCATABLE :: Hspace(:)
    real(8)                      :: T   ! temperature
    real(8)                      :: Z   ! normolized partition function
    real(8)                      :: Eg


    integer:: print = 6
  contains
    procedure,pass::Initialization
    final::Finalization

    procedure,pass::diagonalization
    procedure,pass::Normolize_Energy_And_GetZ
    procedure,pass::get_operator_act_on_state
    procedure,pass::get_product
    procedure,pass::get_State_product
    procedure,pass::get_trace_value
    !
    procedure,pass::printE
    procedure,pass::get_Nsub
    procedure,pass::get_Table_pointer
    procedure,pass::get_EigenId
    procedure,pass::get_E
    procedure,pass::get_Eg
    procedure,pass::GetRz
    procedure,pass::get_T
    procedure,pass::Get_Ns
    procedure,pass::get_Vec
    procedure,pass::SynchronizeWithHamiltonian
  endtype


  private::Initialization,UnInitialization,Finalization
  !
  private::get_operator_act_on_state
  private::diagonalization,Normolize_Energy_And_GetZ,get_product,get_trace_value
  private::get_State_product
  !
  !
  private::printE
  private::get_EigenId,get_T,Get_Ns,get_Nsub,get_Eg,get_Table_pointer,get_E,get_Vec,GetRz
  private::SynchronizeWithHamiltonian


contains

  subroutine Initialization(self,T,Ta,CH,IsReal,print_)
             implicit none
             class(ED_GCE),intent(inout)::self
             real(8),intent(in)            ::T
             class(table),target           ::Ta
             class(Ham),target             ::CH
             logical,intent(in),optional   ::IsReal
             integer,intent(in),optional   ::print_
             !-----------------------------------------
             integer::jc

             call UnInitialization(self)
             self%initiated = .true.
             self%T     = T
             self%Ta    => Ta
             if (present(IsReal)) self%IsReal=IsReal
             self%CH    => CH
             if (present(print_)) self%print = print_

             !---check T----
             if (T.le.0.0_8)then
                write(self%print,*)"Finite temperature T should >0. Now T=",T
                stop
             endif



             self%Nsub = self%Ta%get_nsub()

             allocate( self%Hspace(self%Nsub)   )
             do jc = 1, self%Nsub
                call self%Hspace(jc)%Initialization(self%Ta,self%IsReal,self%CH,jc,self%print)
             enddo

  endsubroutine


  subroutine UnInitialization(self)
              implicit none
              class(ED_GCE),intent(inout)::self
              !------------------------------------
              if (self%initiated)then
                 self%initiated = .false.
                !  do jc = 1, self%Nsub
                !     call self%Hspace(jc)%UnInitialization()
                !  enddo
                 deallocate(  self%Hspace )
                 self%Ta => null()
                 self%CH => null()
              endif
  endsubroutine

  subroutine Finalization(self)
             implicit none
             type(ED_GCE),intent(inout)::self
             !------------------------------------
             call UnInitialization(self)
  endsubroutine





  subroutine diagonalization(self)
             implicit none
             class(ED_GCE),intent(inout)::self
             !------------------------------------
             integer(8)::jc
             do jc = 1_8 ,self%Nsub
                call self%Hspace(jc)%diagonalization()
             enddo
             self%Normolized = .false.
             self%EigenId = self%CH%GetEigenId()
  endsubroutine


  subroutine Normolize_Energy_And_GetZ(self)
            implicit none
            class(ED_GCE),intent(inout)::self
            !------------------------------------
            real(8)::Eg,T
            integer::jc,jcs
            if (.not.self%Normolized)then
               !---------------------find out groundstate energy---
               self%Eg = self%Hspace(1)%get_E(1)
               do jc = 2, self%Nsub
                  Eg = self%Hspace(jc)%get_E(1)
                  if (Eg.le.self%Eg)  self%Eg = Eg
               enddo
              !----------------------normolize energy--------------
              do jc = 1, self%Nsub
                 call self%Hspace(jc)%NormolizedEnergy(self%Eg)
              enddo
              !----------------------get Z-------------------------
              T = self%T
              self%Z = 0._8
              do jc = 1,self%Nsub
                 do jcs = 1 , self%Ta%get_sub_d(jc)
                    Eg = self%Hspace(jc)%get_E(jcs)
                    self%Z = self%Z + dexp( - Eg / T   )
                 enddo
              enddo
              !------------
              self%Normolized = .true.
            endif
  endsubroutine


  integer function get_EigenId(self)
              implicit none
              class(ED_GCE),intent(inout)::self
              !-----------------------------------
              get_EigenId = self%EigenId
            endfunction


  function get_operator_act_on_state(self,oper,subid,index,d) result(state)
    implicit none
    class(ED_GCE),intent(inout)::self
    class(FermOper),intent(inout)::oper
    integer,intent(in)::subid,index,d
    complex*16::state(d)
    !---------------------------------------------
    state = self%Hspace(subid)%get_act(oper,index,d)
  endfunction

  function get_product(self,subi,i,subj,j,oper) result(returnV)
             implicit none
             class(ED_GCE),intent(inout)::self
             integer,intent(in)::subi,i,subj,j
             class(FermOper)::oper
             complex*16::returnV
             !------------------------------------------
             integer::subiD
             complex*16,allocatable::s(:)
             subiD = self%Ta%get_sub_d(subi)
             allocate(s(subiD))
             s = self%Hspace(subj)%get_act(oper,j,subiD)
             returnV = self%Hspace(subi)%get_product(s,i)
             deallocate( s  )
  endfunction

  complex*16 function get_State_product(self,subi,i,d,state)
    implicit none   !< subi,i | state  >
    class(ED_GCE),intent(inout)::self
    integer,intent(in)::subi,i,d
    complex*16,intent(in)::state(d)
    !-----------------------------------------
    get_State_product = self%Hspace(subi)%get_product(state,i)
  endfunction






  function get_trace_value(self,opert) result(returnV)
             implicit none
             class(ED_GCE),intent(inout)::self
             class(FermOper),intent(inout)::opert
             complex*16::returnV
             !---------------------------------------
             integer::jcsub,jcs
             complex*16::temp
             real(8)::E
             if (self%initiated .and.  self%Normolized)then
               returnV = (0._8,0._8)
               do jcsub = 1, self%Nsub
                  do jcs = 1, self%Ta%get_sub_d(jcsub)
                     temp = get_product(self,jcsub,jcs,jcsub,jcs,opert)
                     E = self%Hspace(jcsub)%get_E(jcs)
                     returnV = returnV + temp * dexp( - E / self%T  )
                  enddo
               enddo
               returnV = returnV / SELF%Z
             else
                write(self%print,*)"get_trace_value is called while the EigenValue is not Normolized."
                stop
             endif
  endfunction

  function get_T(self) result(T)
    implicit none
    class(ED_GCE),intent(inout)::self
    real(8)::T
    !---------------------------------------
    T = self%T
  endfunction

  integer function Get_Ns(self)
    implicit none
    class(ED_GCE),intent(inout)::self
    !---------------------------------------
    Get_Ns= self%ta%get_ns()
  endfunction


  real*8 function get_Eg(self)
    implicit none
    class(ED_GCE),intent(inout)::self
    !---------------------------------------
    get_Eg = self%Eg
  endfunction

  integer function get_Nsub(self)
    implicit none
    class(ED_GCE),intent(inout)::self
    !---------------------------------------
    get_Nsub = self%nsub
  endfunction



  subroutine printE(self)
             implicit none
             class(ED_GCE),intent(inout)::self
             !------------------------------------
             integer(8)::jc
            ! WRITE(*,*)"HERE"
            write(self%print,*)"Eg=",self%Eg
             do jc = 1_8, self%Nsub
                call self%Hspace(jc)%printE()
             enddo
  endsubroutine


  function get_Table_pointer(self) RESULT(r)
          implicit none
          class(ED_GCE),intent(inout)::self
          class(table),pointer::r
          !------------------------------------
          r => self%Ta
        endfunction


  function get_E(self,subid,i) result(E)
          implicit none
          class(ED_GCE),intent(inout)::self
          integer,intent(in)::subid
          integer,intent(in)::i
          real(8)::E
          !---------------------------
          E =  self%Hspace(subid)%get_E(i)
  endfunction



  function get_Vec(self,subid,i,d) result(V)
    implicit none
    class(ED_GCE),intent(inout)::self
    integer,intent(in)::subid
    integer,intent(in)::i,d
    complex*16::V(d)
    !---------------------------
    v = self%Hspace(subid)%get_Vec(i)
endfunction



  real*8 function GetRz(self)
    implicit none
    class(ED_GCE),intent(inout)::self
    !--------------------------------------
    if (self%initiated .and.  self%Normolized)then
      GetRz = self%Z
    else
      write(self%print,*)"get_Z is called while the eigensystem is not normolized."
      stop
    endif
  endfunction

  subroutine SynchronizeWithHamiltonian(self)
    implicit none
    class(ED_GCE),intent(inout)::self
    !--------------------------------------
    if (  self%EigenId == self%CH%GetEigenId()  )then
      continue
    else
      call diagonalization(self)
    endif
    if (.not.self%Normolized)then
      call Normolize_Energy_And_GetZ(self)
    endif
  endsubroutine















endmodule
