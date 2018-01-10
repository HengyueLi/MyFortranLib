

!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : ED_Subspace
! OBJECT: TYPE(EDSubSpace)
! USED  : fermion_table,FermionHamiltonian,class_numerical_method
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
!                  [sub] Initialization(Ta,IsReal,H,subid,*PRINT)
!                        type(table)::Ta  (will be pointed and be read only)
!                        logical    ::IsReal ! if H is a real problem.
!                        type(Ham)  ::H , will be  pointed.
!                        integer    ::subid
!                        Ta shoule be initiated before this subroutine since some information would be used.
!                  [sub] NormolizedEnergy(Eg)
!                        real(8)::Eg
!                        E -> E -Eg
!                  [sub] diagonalization()
!
! avalable gets:
!
!                  [fun] get_E(i)
!                        integer::i   ;   real(8)::get_E
!                  [fun] get_Vec(i)
!                        integer(8)::i     ;complex*16::get_Vec(d)  !where d is the dimension of the subspace.
!                  [fun] get_act(opt,i,d)
!                        type(FermOper)::opt
!                        integer(8)::i
!                        integer(8)::d
!                        complex*16::get_act(d)
!                        for a input operator opt%act
!                        opt%act |i> = S(d)    ,  thus the dimension of outputstate should be input.
!                  [fun] get_product(S,i)
!                        complex*16::S(self%d),get_product
!                        integer::i to identify the i-th eignstate in this subspace.
!                        returnV = <i|S> .      The input array S should be the same dimension of this subspace.
!
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

module ED_Subspace
   use fermion_table
   use FermionHamiltonian
   use class_numerical_method
   implicit none

   TYPE::EDSubSpace
       private
       logical                ::initiated  = .false.
       logical                ::Normolized = .false.
       logical                ::IsReal     = .false.
       integer                ::subid      = -1_8
       integer                ::d          = -1_8
       class(table),pointer   ::T          => null()
       class(Ham),pointer     ::Hi         => NULL()
       complex*16,allocatable ::Hs(:,:)
       real*8    ,allocatable ::E(:)
       real*8::Eg ! after normolized, this will save the ground state energy.
       integer::print = 6
   contains
       procedure,pass::Initialization
       procedure,pass::UnInitialization
       final::Finalization
       !
       procedure,pass::diagonalization
       procedure,pass::NormolizedEnergy
       procedure,pass::get_E
       procedure,pass::get_act
       procedure,pass::get_product
       procedure,pass::get_Vec
       !
       !
       procedure,pass::printE
   end type
   !
   private::Initialization,UnInitialization,Finalization
   !
   private::set_Hamiltonian_matrix_element,diagonalization,NormolizedEnergy,get_E
   private::get_act,get_product,printE,get_Vec
   !

 contains

  subroutine Initialization(self,T,IsReal,H,subid,print_)
    implicit none
    class(EDSubSpace),intent(inout)::self
    class(table),target,intent(inout)::T
    logical,intent(in)::IsReal
    class(Ham),target,intent(inout)::H
    integer,intent(in)::subid
    integer,intent(in),optional::print_
    !---------------------------------------
    call UnInitialization(self)
    self%initiated = .true.
    !----- check if T is initiated
    if (.not. ( T%Is_initiated()   ) )then
       write(self%print,*)"When initialize EDSubSpace, the input Table should be Initiated."
       stop
    endif
    self%T     => T
    self%IsReal= IsReal
    self%Hi    => H
    self%subid = subid
    self%d     = self%T%get_sub_d(self%subid)

    if (present(print_)) self%print = print_
    allocate(  self%Hs(self%d,self%d)   )
    allocate(  self%E(self%d)           )

  endsubroutine


  subroutine UnInitialization(self)
    implicit none
    class(EDSubSpace),intent(inout)::self
    !---------------------------------------
    if (self%initiated)then
       self%Initiated = .false.
       self%subid     = -1_8
       self%d         = -1_8
       self%T         => null()
       self%Hi        => NULL()
       deallocate( self%Hs  )
       deallocate( self%E   )
    endif
  endsubroutine

  subroutine Finalization(self)
    implicit none
    type(EDSubSpace),intent(inout)::self
    !---------------------------------------
    call UnInitialization(self)
  endsubroutine







  subroutine set_Hamiltonian_matrix_element(self)  ! only upper part of H will be setted.
    implicit none
    class(EDSubSpace),intent(inout)::self
    !---------------------------------------
    integer*8::basisIn,basisOut
    integer::jcopra,jci,jcj,Sign_,subid
    complex*16::v
    class(FermOper),pointer::opr
    self%Hs = (0._8,0._8)
    subid   = self%subid
    do jcopra = 1_8 , self%hi%GetOptN()
       opr    => self%hi%GetOptact(jcopra)
       v      =  self%Hi%GetOptV(jcopra)


                                        ! write(*,"(  (I2),',',(I2) )")jcopra,self%hi%GetOptN()
                                        ! write(*,"(I2)")opr%get_optid()
                                        ! write(*,*)opr%get_para()
                                        ! write(*,"(f6.3)")real(v)
                                        ! write(*,*)"-------------------"



       !------------------------------
       !    for ( opr  ,   v   )
       do jcj = 1 , self%d
          !-------- for subbasi  |jcs>---
          basisIn = self%T%get_subindex_to_basis(subid,jcj)
          call opr%act(basisIn,basisOut,Sign_)
          if (basisOut.ne.-1_8)then
             jci = self%T%get_basis_to_sub_index(basisOut)
             if ( jci.le.jcj  )then
               self%Hs(jci,jcj) = self%Hs(jci,jcj) + V * Sign_
             endif
          endif
       enddo
       !---------------------------------
   enddo
  endsubroutine


  subroutine diagonalization(self)
    implicit none
    class(EDSubSpace),intent(inout)::self
    !---------------------------------------
    TYPE(nummethod)::num
    real*8,allocatable::Hr(:,:)

    call set_Hamiltonian_matrix_element(self)

    !--------------
    if (self%IsReal)then
      allocate(Hr(self%d,self%d))
      Hr = Real( self%Hs ,8 )
      call num%ED_Hermitian_matrix("U",self%d,Hr     ,self%E)
      self%Hs = cmplx(Hr,0._8,8)
      deallocate( Hr  )
    else
      call num%ED_Hermitian_matrix("U",self%d,self%Hs,self%E)
    endif
    !--------------
    self%Normolized = .false.                           

  endsubroutine


  subroutine NormolizedEnergy(self,Eg)
            implicit none
            class(EDSubSpace),intent(inout)::self
            real(8),intent(in)::Eg
            !---------------------------------------
            integer(8)::jc
            if (.Not.self%normolized)then
               do jc = 1_8 , self%d
                  self%E(jc) = self%E(jc) -Eg
               enddo
               self%Eg = Eg
               self%normolized = .true.
            endif
  endsubroutine


  function get_E(self,i) result(E)
            implicit none
            class(EDSubSpace),intent(inout)::self
            integer,intent(in)::i
            real*8::E
            !---------------------------------------
            E = self%E(i)
  endfunction


  function get_Vec(self,i) result(v)
    implicit none
    class(EDSubSpace),intent(inout)::self
    integer,intent(in)::i
    complex*16::V(self%d)
    !---------------------------------------
    V = self%Hs(:,i)
endfunction



  function get_act(self,opt,i,d)  result(s)
            implicit none
            class(EDSubSpace),intent(inout)::self
            class(FermOper),intent(inout)::opt
            integer,intent(in)::i,d
            complex*16::s(d)
            !---------------------------------------
            integer::jc,sign_,id,subid
            integer*8::bi,bo
            subid = self%subid
            s = (0._8,0._8)
            do jc = 1 ,self%d
               bi = self%T%get_subindex_to_basis(subid,jc)
               call opt%act(bi,bo,sign_ )
               if ( bo.ne.-1_8 )then
                  id =  self%T%get_basis_to_sub_index(bo)
                  s(id) = s(id) + sign_ * self%Hs(jc,i)
               endif
            enddo
  endfunction


  function get_product(self,S,i) result(returnV)
            implicit none
            class(EDSubSpace),intent(inout)::self
            complex*16,intent(in)::s(self%d)
            integer,intent(in)::i
            complex*16::returnV
            !---------------------------------------
            returnV = sum( conjg(self%Hs(:,i))*S )
  endfunction






  subroutine printE(self)
            implicit none
            class(EDSubSpace),intent(inout)::self
            !---------------------------------------
            integer::jc
            write(self%print,*) self%T%get_subspace_marker(self%subid)
            do jc = 1 , self%d
              write(self%print,*)jc-1,self%e(jc)
            enddo
            ! write(self%print,*) self%E
  endsubroutine


end module
