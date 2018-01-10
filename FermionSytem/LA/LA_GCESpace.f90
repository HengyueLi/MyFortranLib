
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  :  MODULE
! NAME  :  LA_GCESpace
! OBJECT:  TYPE(LA_GCE)
! USED  :  LA_Subspace,basic_math_functions
! DATE  :  2018-01-07
! AUTHOR:  hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            Grand Cononical Emsenbel system. Lanczos method
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
!                  ![sub] Initialization(Ta,H,IsReal_,PRINT_,SHOW_,Pre_,DegPre_,Bzero_,M_,oth_)
!                        type(table)::Ta    ! must be initiated first
!                        type(Ham)::H  !
!                        logical::IsReal_ = .false.
!                        integer::PRINT_,SHOW_
!                        real*8::Pre_,DegPre_,Bzero_
!                        integer::M_
!                        logical::oth_
!
!                  ![sub] Normolize_Energy_And_GetZ()
!
! avalable gets:
!                  [fun] GetSubDe(i)
!                        degeneracy in subspace
!                  [fun] GetDe()
!                         integer    degeneracy
!                  [fun] GetEg()
!                         real*8  ground state energy
!                  [fun] GetSubEg(i)
!
!                  [fun] GetNsub()
!                        return integer
!                  [fun] GetSubSpacePointer(i)
!                        return class(LASubSpace),pointer
!                  [fun] GetNsubGS()
!                        return the number of subpsaces that contains GS
!                  [fun] GetSubIdGs(i)
!                        return subid of subspace that contains GS
!
!                  [fun] GetOperateProduct(A)
!                        TYPE(FermOper)::A
!                        return complex*16::<A> = 1/D * sum_{all GS} <GS|A|GS>
! avalable is :
!                  ![fun] i
! others      :
!                  [sub] SynchronizeWithHamiltonian()
!                        Check the EigenSytem of LA. If Hamiltonian is changed, renew EigenSytem
!                        Here renew means diagonalization +  Normolize_Energy_And_GetZ
!
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$






module Data20171205
  implicit none
  type::data
     integer::subid
  endtype
endmodule
module List20171205
  use Data20171205
  INCLUDE "../../ListStructure/ListStructure.finc"
endmodule


module LA_GCESpace
  use List20171205,only:sublist => ListStru , data
  use fermion_table
  use FermionHamiltonian
  use LA_Subspace
  implicit none


  TYPE::LA_GCE
    private
    logical:: Initiated  =  .false.
    class(table),pointer::Ta => null()
    class(Ham),pointer :: H  => null()
    integer::EigenId  = -20171204
    logical::isreal   = .false.
    integer::Nsub
    type(LASubSpace),allocatable::SUB(:)
    !--------------------------------------
    real*8::Eg
    type(sublist)::GSl
    integer::De


    !----------
    real*8 ::pre      = 1.e-14
    real*8 ::bzero    = 0.000001
    integer::M        = 30
    logical::oth      = .true.
    real*8 ::DegPre   = 1.e-6
   !------------
    integer::print = 6
    integer::show  = 0

   contains
     procedure,pass::Initialization
     final::Finalization

     procedure,pass::SynchronizeWithHamiltonian
     procedure,pass::GetDe
     procedure,pass::GetSubDe
     procedure,pass::GetEg
     procedure,pass::GetSubSpacePointer
     procedure,pass::GetNsub
     procedure,pass::GetNsubGS
     procedure,pass::GetSubIdGs
     procedure,pass::GetSubEg
     procedure,pass::GetOperateProduct
  endtype



  private::Initialization,UnInitialization,Finalization

  private::SynchronizeWithHamiltonian


  private::GetDe,GetEg
  private::GetSubSpacePointer,GetSubPointer
  private::GetNsub
  private::GetNsubGS
  private::GetSubIdGs
  private::GetSubDe
  private::GetSubEg
  PRIVATE::GetOperateProduct

contains

!,PreE_,bzero_,M_,oth_

  subroutine Initialization(self,Ta,H,IsReal_,PRINT_,SHOW_,Pre_,DegPre_,Bzero_,M_,oth_)
    implicit none
    class(LA_GCE),intent(inout)    :: Self
    class(table),intent(in),target :: Ta
    class(Ham),intent(in),target   :: H
    logical,intent(in),optional    :: IsReal_
    INTEGER,intent(in),optional    :: PRINT_,SHOW_
    real*8,intent(in),optional     :: Pre_,DegPre_,Bzero_
    integer,intent(in),optional    :: M_
    logical,intent(in),optional    :: oth_
    !----------------------------------------------
    integer::jc
    !-----------
    call UnInitialization(self) ; self%initiated = .true.
    self%Ta => Ta
    self%H  => H
    if(present(IsReal_))  self%isreal = IsReal_
    if (present(PRINT_))  self%print  = print_
    if (present(SHOW_))   self%show   = SHOW_
    if (present(DegPre_)) self%DegPre = DegPre_
    if (present(Pre_))    self%Pre    = Pre_
    if (present(Bzero_))  self%bzero  = Bzero_
    if (present(M_))      self%M      = M_
    if (present(oth_))    self%oth    = oth_
    !-------------------------------------------------
    if (.not.self%Ta%Is_initiated()) then
      write(self%print,*)"Table should be initiated befora initiate LA_GCE";stop
      stop
    endif
    !-------------------------------------------------
    self%nsub = self%Ta%get_Nsub()
    allocate(self%sub(self%nsub))
    do jc = 1 , self%nsub
       call self%sub(jc)%Initialization(Ta=self%ta,IsReal=self%IsReal,&
                H=self%H,subid=jc,PRINT_=self%print,show_=self%show,  &
                PreE_=self%pre,PreDe_=self%DegPre,bzero_=self%bzero,M_=self%m,Oth_=self%oth)
    enddo
  endsubroutine

  subroutine UnInitialization(self)
    implicit none
    class(LA_GCE),intent(inout)::self
    !----------------------------------------------
    if (self%initiated)then
       self%initiated = .false.
       !--------------------------
       deallocate(self%sub)
    endif
  endsubroutine

  subroutine Finalization(self)
    implicit none
    type(LA_GCE),intent(inout)::self
    !----------------------------------------------
    call UnInitialization(self)
  endsubroutine









  subroutine SynchronizeWithHamiltonian(self)
    use basic_math_functions
    implicit none
    class(LA_GCE),intent(inout)::self
    !----------------------------------------------
    integer::jc,con
    real*8::Eg
    TYPE(bmathf)::f
    class(data),pointer::p

    if (self%EigenId == self%H%GetEigenId()) goto 999
    Eg=Huge(Eg)
    !--------------------------------------------------
    ! diagonalization
    self%EigenId = self%H%GetEigenId()
    call self%SUB%SynchronizeWithHamiltonian()
    !--------------------------------------------------
    ! find sub
    call self%GSl%Initialization(self%print)
    ! find out one of the GS space
    do jc = 1 , self%nsub                     !;write(*,*)jc,self%sub(jc)%GetEg()
       if (self%sub(jc)%GetEg().le.Eg) Eg = self%sub(jc)%GetEg()
    enddo
    !-------------
    self%eg = eg
    !-count
    self%De     = 0
    do jc = 1 , self%nsub
      if (  f%IsTwoValuePercentageTheSame(self%sub(jc)%GetEg(),Eg,self%DegPre))then
          allocate(p) ; p%subid = jc
          call self%GSl%append(p)
          self%De = self%De + self%sub(jc)%GetDe()
      endif
    enddo

999 continue
  endsubroutine












  integer function GetDe(self)
    implicit none
    class(LA_GCE),intent(inout)::self
    !----------------------------------------------
    GetDe = self%de                                !;write(*,*)self%de ,self%nsub
  endfunction

  real*8 function GetEg(self)
    implicit none
    class(LA_GCE),intent(inout)::self
    !----------------------------------------------
    GetEg = self%Eg
  endfunction





  function GetSubPointer(subin) result(r)
    implicit none
    class(LASubSpace),target,intent(inout)::subin
    class(LASubSpace),pointer::r
    !--------------------------------
    r => subin
  endfunction


  function GetSubSpacePointer(self,i) result(r)
    implicit none
    class(LA_GCE),intent(inout)::self
    integer,intent(in)::i
    class(LASubSpace),pointer::r
    !----------------------------------------------
    r => GetSubPointer(self%SUB(i))
  endfunction


  integer function GetNsub(self)
    implicit none
    class(LA_GCE),intent(inout)::self
    !-------------------------------------
    GetNsub = self%nsub
  endfunction


  integer function GetNsubGS(self)
    implicit none
    class(LA_GCE),intent(inout)::self
    !-------------------------------------
    GetNsubGS = self%GSl%getlen()
  endfunction


  integer function GetSubIdGs(self,i)
    implicit none
    class(LA_GCE),intent(inout)::self
    integer,intent(in)::i
    !-------------------------------------
    class(data),pointer::p
    p => self%GSL%GetDataPointer(i)
    GetSubIdGs = p%subid
  endfunction


  integer function GetSubDe(self,i)
    implicit none
    class(LA_GCE),intent(inout)::self
    integer,intent(in)::i
    !-------------------------------------
    GetSubDe = self%SUB(i)%GetDe()
  endfunction

  real*8 function GetSubEg(self,i)
    implicit none
    class(LA_GCE),intent(inout)::self
    integer,intent(in)::i
    !-------------------------------------
    GetSubEg = self%SUB(i)%GetEg()
  endfunction


  complex*16 function GetOperateProduct(self,opt)
    implicit none
    class(LA_GCE),intent(inout)::self
    class(FermOper),intent(inout)::opT
    !-----------------------------------------
    integer::GsSubid,jc
    GetOperateProduct = (0._8,0._8)

    do jc = 1 , self%GetNsubGS()
      GsSubid = self%GetSubIdGs(jc)
      GetOperateProduct = GetOperateProduct + &
                         self%sub(GsSubid)%GetOperateProduct(opt) * self%sub(GsSubid)%GetDe()
    enddo
    GetOperateProduct = GetOperateProduct / self%GetDe()
  endfunction






endmodule
