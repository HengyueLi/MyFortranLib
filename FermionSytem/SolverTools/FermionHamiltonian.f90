!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : FermionHamiltonian
! OBJECT: TYPE(Ham)
! USED  : FermionOperators,functionalsubs,basic_math_functions
! DATE  : 2018-01-05
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            setting Hamiltonian. Offer interacting for ED and Lanczos
! STANDARD:
!            *call Initialization( ns , *print=6 )
!            *call StartAppendingInteraction()
!            *call AppendingInteraction(InterType,InterPara,InterV)
!            *call EndAppendingInteraction()
!
!           Interaction type is not case sensetive.
!          ┏━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
!          ┃  Interaction type ┃ index parameter     ┃            DESCRIPTION                                                        ┃
!          ┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
!          ┃    "SpinOnSite"   ┃    i,i,spin         ┃  spin depedent onsite energy:   v * n_{i,spin}                                ┃
!          ┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
!          ┃    "OnSite"       ┃    i,i              ┃          onsite energy:   v *\sum_{spin} n_{i,spin}                           ┃
!          ┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
!          ┃  "GlobalOnSite"   ┃                     ┃   global onsite energy:   v *\sum_{(i=1,Np),spin} n_{i,spin}                  ┃
!          ┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
!          ┃  "DiffOnSite"     ┃    i, j             ┃ Set a different on two orbital: v * ( N_i - N_j ),N_i = \sum_spin n_{i,spin}  ┃
!          ┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
!          ┃ "SpinHopping"     ┃     i,j,spin        ┃ spin hopping: v * c^+_{i,spin}c_{j,spin} +h.c.                                ┃
!          ┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
!          ┃     "Hopping"     ┃     i,j             ┃   hopping1: v*\sum_{spin} c^+_{i,spin}c_{j,spin} +h.c.                        ┃
!          ┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
!          ┃    "OnSiteU"      ┃     i,i             ┃   onsite Comlomb:   v * n_{i,up}*n_{i,down}                                   ┃
!          ┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
!          ┃ "GlobalOnSiteU"   ┃                     ┃  global onsite Comlomb:   v *\sum_i n_{i,up}*n_{i,down}                       ┃
!          ┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
!          ┃   "InterV"        ┃     i,j             ┃   intersite Comlomb: v * \sum_{si} * n_{i,si} * \sum_{sj}n_{j,sj}             ┃
!          ┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
!          ┃   "PairHopping"   ┃     i,j             ┃   pair hopping:V * c^+_{i up} * c^+_{i do} * c^{j do} * c_{j up} +h.c.        ┃
!          ┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
!          ┃   "Hund"          ┃     i,j             ┃ Hund:sum_{spin,spin'} V * c^+_{ispin} * c^+_{jspin'} * c_{ispin'} * c_{jspin} ┃
!          ┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
!          ┃  "GeneralInter1"  ┃ i,si,j,sj,k,sk,l,sl ┃ V*c^+_{ispini} * c^+_{jspinj} * c_{kspink} * c_{lspinl}                       ┃
!          ┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
!          ┃  "GeneralInter2"  ┃ i,si,j,sj,k,sk,l,sl ┃ c^+_{ispini} * c_{jspinj} * c^+_{kspink} * c_{lspinl}                         ┃
!          ┗━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  [sub] Initialization( ns  , *print=6 )
!
!                  [sub] StartAppendingInteraction()
!
!                  [sub] AppendingInteraction(InterType,InterPara,InterV)
!                        character(16)::InterType
!                        integer      ::InterPara(8)
!                        complex*16   ::InterV
!
!                  [sub] EndAppendingInteraction()
! avalable gets:
!                  [fun] GetOptN()
!                        return integer
!                  [fun] GetOptact(i)
!                        return pointer TYPE(FermOper)
!                  [fun] GetOptV(i)
!                        complex*16
!                  [fun] GetEigenId()
!                        return integer
! avalable  IS:
!                  [fun] Is_Usable()
!                        check if all the Hamiltonian terms are set.
!
!  Update 2017-12-24 :
!                  interface of  "SpinOnSite" is changed.
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




!     Initialization   ->(1)->   startappend -> (2) -> appending ->(2)    ->   endappend           ->  (3)
!                                allocate Int                                 allocateOPT


module HamUseData
  implicit none
  type::InterApp
      character(16)::T
      integer      ::P(8)
      complex*16   ::V
  endtype
endmodule HamUseData
module HamUseList
  use HamUseData,only: data =>  InterApp
  INCLUDE "../../ListStructure/ListStructure.finc"
endmodule HamUseList





module HamUseOPT
  use FermionOperators
  implicit none
  type::optApp
      type(FermOper)::opt
      complex*16    ::V
  endtype
endmodule HamUseOPT
module HamUseOptList
  use HamUseOPT,only: data =>  optApp
  INCLUDE "../../ListStructure/ListStructure.finc"
endmodule HamUseOptList





module FermionHamiltonian
       use FermionOperators
       use HamUseList ,only: Inlist => ListStru ,idata => data
       use HamUseOptList,only:Olist => ListStru ,odata => data
       implicit none



       type::Ham
          private
          logical::initiated = .false.
          integer::ns
          !-----------------------------------------------
          type(Olist)::optl
          !-----------------------------------------------
          type(Inlist)::l
          !-----------------------------------------------
          integer::EigenId     = -1     ! Identify a Hamiltonian .
          !-----------------------------------------------
          integer::state       = 0
          real*8 ::IntZero     = 1.e-14
          integer::print       = 6

       contains
         procedure,pass::Initialization
         FINAL::Finalization

         procedure,pass::StartAppendingInteraction
         procedure,pass::AppendingInteraction
         procedure,pass::EndAppendingInteraction

         procedure,pass::Is_Usable
         procedure,pass::GetOptN
         procedure,pass::GetOptact
         procedure,pass::GetOptV
         procedure,pass::GetEigenId


         procedure :: Copy
         generic :: assignment(=) => Copy
       endtype



       private::Initialization,UnInitialization,Finalization

       private::StartAppendingInteraction,AppendingInteraction,EndAppendingInteraction

       private::AppendOPT

       private::App_SpinOnsite,App_OnSite,App_GlobalOnSite
       private::App_SpinHopping,App_Hopping
       private::App_OnSiteU,App_GlobalOnSiteU,App_InterV
       private::App_PairHopping,App_Hund

       private::App_GeneralInter1,App_GeneralInter2

       private::CheckParaTheSame

       private::Is_Usable,GetOptN,GetOptact,GetOptV
       private::GetEigenId

       private::Copy


     contains


       subroutine Initialization(self,Ns,print_)
         implicit none
         class(Ham),intent(inout)::self
         integer,intent(in)::Ns
         integer,intent(in),optional::print_
         !---------------------------------------------
         call UnInitialization(self)
         self%initiated = .true.
         self%ns = ns
         if (present(print_))    self%print    = print_
         !--------------
         self%state = 1

       endsubroutine

       subroutine UnInitialization(self)
         implicit none
         class(Ham),intent(inout)::self
         !----------------------------------
         if (self%initiated)then
            self%initiated = .false.
         endif
       endsubroutine

       subroutine Finalization(self)
         implicit none
         type(Ham),intent(inout)::self
         !----------------------------------
         call UnInitialization(self)
       endsubroutine



       subroutine StartAppendingInteraction(self)
         implicit none
         class(Ham),intent(inout)::self
         !----------------------------------
         select case(self%state)
         case(0)
           write(self%print,*)"ERROR: Ham is not Initialized when calling 'StartAppendingInteraction'"
           stop
         case(1,2,3)
         case DEFAULT
           write(self%print,*)"ERROR: Unkonw state @ Ham @ StartAppendingInteraction"
           stop
         END SELECT
         self%state = 2
         call self%l%Initialization(print_=self%print)
       endsubroutine

       subroutine AppendingInteraction(self,InterType,InterPara,InterV)
         use functionalsubs
         implicit none
         class(Ham),intent(inout)::self
         character(16),intent(in)::InterType
         integer,intent(in)::InterPara(8)
         complex*16,intent(in)::InterV
         !----------------------------------
         character(16)::InType
         class(idata),pointer::d

         TYPE(funcsubs)::f
         if (self%state.ne.2)then
           write(self%print,*)"ERROR: AppendingInteraction is not allowed."
           write(self%print,*)"Sub StartAppendingInteraction should be called first." ;stop
         endif
         !--------------------------check NonZero------------------------------------
         if (zabs(InterV).le.self%IntZero) goto 999
         !---------------------------------------------------------------------------
         InType = f%get_string_upper(InterType)
         allocate(d)
         d%T = InType
         d%p = InterPara
         d%v = InterV
         call self%l%append(d)
         !------------------------------------------------------------------------------------------------

    999  continue
       endsubroutine


       SUBROUTINE   EndAppendingInteraction(self)
         USE basic_math_functions
         implicit none
         class(Ham),intent(inout)::self
         !------------------------------------------
         TYPE(bmathf)::F
         integer::ns,i,j,k,l,spini,spinj,spink,spinl,jc
         integer::para(8),optpara(8),print
         character(16)::Intype
         complex*16::V
         class(idata),pointer::p

         ns    = self%ns
         print = self%print

         if (self%state.ne.2)then
            write(self%print,*)"ERROR: StartAppending is not called while&
                 direct calling of EndAppendingInteraction is attempted."  ;stop
         endif


         call self%optl%Initialization(print_=self%print)


          do jc = 1 , self%l%GetLen()
             p => self%l%GetDataPointer(jc)
             Intype = p%T
             para   = p%p
             V      = p%v

            !--------------------------------------
            select case(trim(adjustl(Intype)))
            case("SPINONSITE")   !----------------SpinOnSite-------------------------------
              call App_SpinOnsite(self,V,para(1),para(3))
            case("ONSITE")       !----------------Onsite    -------------------------------
              call App_OnSite(self,V,para(1))
            case("GLOBALONSITE") !----------------global Onsite---------------------------
              CALL App_GlobalOnSite(self,V)
            case("DIFFONSITE")   !----------------set a gap  -----------------------------
              call App_SetOrbitalGap(self,para(1),para(2),V)
            case("SPINHOPPING")  !----------------spin hopping ---------------------------
              call App_SpinHopping(self,V,para(1),para(2),para(3))
            case("HOPPING")      !---------------Hopping       ---------------------------
              call App_Hopping(self,V,para(1),para(2))
            case("ONSITEU")      !--------------   U          ----------------------------
              call App_OnSiteU(self,V,para(1))
            case("GLOBALONSITEU")!----------------global U    ----------------------------
              call App_GlobalOnSiteU(self,V)
            case("INTERV")       !---------------intger site V ---------------------------
              call App_InterV(self,V,para(1),para(2))
            case("PAIRHOPPING")  !--------------- pair hopping  --------------------------
              call App_PairHopping(self,V,para(1),para(2))
            case("HUND")         !--------------  Hund term ------------------------------
              call App_Hund(self,V,para(1),para(2))
            case("GENERALINTER1")!---------------GENERAL1 --------------------------------
              call App_GeneralInter1(self,V,para)
            case("GENERALINTER2")!---------------GENERAL2 --------------------------------
              call App_GeneralInter2(self,V,para)
            case default
              write(self%print,*)"ERROR: Unkonw interacting type in Ham:",trim(adjustl(Intype));stop
            endselect
            !--------------------------------------
         enddo

         self%state = 3
         self%EigenId = f%get_random_int(0,2147483646)
       endsubroutine

  !      subroutine AppendOPT(self,optid,para,V)
  !        implicit none
  !        class(Ham),intent(inout)::self
  !        integer,intent(in)::optid , para(8)
  !        complex*16,intent(in)::v
  !        !------------------------------------------
  !        class(odata),pointer::p
  !        integer::jc,opara(8)
  !        do jc = 1 , self%optl%GetLen()
  !           p => self%optl%GetDataPointer(jc)
  !           if (optid  .eq.  p%opt%get_optid()   )then
  !             opara = p%opt%get_para()
  !             if ( CheckParaTheSame(opara, para)  )then
  !                p%v = p%v + v
  !                !p => null()
  !                goto 999
  !             endif
  !           endif
  !        enddo
  !        allocate(p)
  !        call p%opt%Initialization(self%ns,optid,para,self%print)
  !        p%v = v
  !        call self%optl%append(p)
  !  999   continue
  !      endsubroutine
  subroutine AppendOPT(self,optid,para,V)
    implicit none
    class(Ham),intent(inout)::self
    integer,intent(in)::optid , para(8)
    complex*16,intent(in)::v
    !------------------------------------------
    class(odata),pointer::p
    integer::jc,opara(8)

    call self%optl%SetMark(1)
    do jc = 1 , self%optl%GetLen()
       p => self%optl%GetMarkedPointerAndNext()
       if (optid  .eq.  p%opt%get_optid()   )then
         opara = p%opt%get_para()
         if ( CheckParaTheSame(opara, para)  )then
            p%v = p%v + v
            !p => null()
            goto 999
         endif
       endif
    enddo
    allocate(p)
    call p%opt%Initialization(self%ns,optid,para,self%print)
    p%v = v
    call self%optl%append(p)
999   continue
  endsubroutine
    logical function CheckParaTheSame(para1,para2)
       implicit none
       integer,intent(in)::para1(8),para2(8)
       !-------------------------------------
       integer::jc
       if (  sum(abs(para1-para2))==0  )then
         CheckParaTheSame = .true.
       else
         CheckParaTheSame = .false.
       endif
    endfunction






      subroutine App_SpinOnsite(self,V,i,spini)
        implicit none
        class(Ham),intent(inout)::self
        complex*16,intent(in)::V
        integer,intent(in)::i,spini
        !---------------------------------
        integer::para(8) = -1
        para(1) = i
        para(2) = spini
        call AppendOPT(self,9,para,V)
      endsubroutine

      subroutine App_OnSite(self,V,i)
        implicit none
        class(Ham),intent(inout)::self
        complex*16,intent(in)::V
        integer,intent(in)::i
        !---------------------------------
        call App_SpinOnsite(self,V,i,0)
        call App_SpinOnsite(self,V,i,1)
      endsubroutine

      subroutine App_GlobalOnSite(self,V)
        implicit none
        class(Ham),intent(inout)::self
        complex*16,intent(in)::V
        !---------------------------------
        integer::jc
        do jc = 0 , self%ns - 1
          call App_OnSite(self,V,jc)
        enddo
      endsubroutine


      subroutine App_SetOrbitalGap(self,i,j,V)
        implicit none
        class(Ham),intent(inout)::self
        integer,intent(in)::i,j
        complex*16,intent(in)::V
        !---------------------------------
        call App_OnSite(self, V,i)
        call App_OnSite(self,-V,j)
      endsubroutine

      subroutine App_SpinHopping(self,V,i,j,spin)
        implicit none
        class(Ham),intent(inout)::self
        complex*16,intent(in)::V
        integer,intent(in)::i,j,spin
        !---------------------------------
        complex*16::Vd
        integer::para(8) = -1
        para(1) = i  ; para(2) = spin  ; para(3) = j ; para(4) = spin
        call AppendOPT(self,3,para,V)
        para(1) = j  ; para(2) = spin  ; para(3) = i ; para(4) = spin
        Vd = conjg(V)
        call AppendOPT(self,3,para,Vd)
      endsubroutine


      subroutine App_Hopping(self,V,i,j)
        implicit none
        class(Ham),intent(inout)::self
        complex*16,intent(in)::V
        integer,intent(in)::i,j
        !---------------------------------
        call App_SpinHopping(self,V,i,j,0)
        call App_SpinHopping(self,V,i,j,1)
      endsubroutine


      subroutine App_OnSiteU(self,V,i)
        implicit none
        class(Ham),intent(inout)::self
        complex*16,intent(in)::V
        integer,intent(in)::i
        !---------------------------------
        integer::para(8) = -1
        para(1) = i   ;  para(2) = 0  ;  para(3) = i  ; para(4) = 1
        call AppendOPT(self,8,para,V)
      endsubroutine

      subroutine App_GlobalOnSiteU(self,V)
        implicit none
        class(Ham),intent(inout)::self
        complex*16,intent(in)::v
        !---------------------
        integer::jc
        do jc = 0 , self%ns -1
          call App_OnSiteU(self,V,jc)
        enddo
      endsubroutine

      subroutine App_InterV(self,V,i,j)
        implicit none
        class(Ham),intent(inout)::self
        complex*16,intent(in)::V
        integer,intent(in)::i,j
        !---------------------------------
        integer::para(8) = -1
        integer::spini,spinj
        do spini = 0 , 1 ;  do spinj = 0 , 1
          para(1) = i ; para(2) = spini  ; para(3) = j ; para(4) = spinj
          call AppendOPT(self,8,para,V)
        enddo            ;  enddo
      endsubroutine

      subroutine App_PairHopping(self,V,i,j)
        implicit none!V * c^+_{i up} * c^+_{i do} * c^{j do} * c_{j up} +h.c.
        class(Ham),intent(inout)::self
        complex*16,intent(in)::V
        integer,intent(in)::i,j
        !---------------------------------
        integer::para(8) = -1
        complex*16::vd
        vd = conjg(v)
        para(1) = i ; para(3) = i ; para(5) = j ; para(7) = j
        para(2) = 0 ; para(4) = 1 ; para(6) = 1 ; para(8) = 0
        call AppendOPT(self,7,para,V)
        para(1) = j ; para(3) = j ; para(5) = i ; para(7) = i
        call AppendOPT(self,7,para,Vd)
      endsubroutine

      subroutine App_Hund(self,V,i,j)
        implicit none!sum_{spin,spin'} V * c^+_{ispin} * c^+_{jspin'} * c_{ispin'} * c_{jspin}
        class(Ham),intent(inout)::self
        complex*16,intent(in)::V
        integer,intent(in)::i,j
        !---------------------------------
        integer::para(8) = -1
        integer::spin1,spin2
        do spin1 = 0 ,1  ; do spin2 = 0 , 1
          para(1) = i     ; para(3) = j     ; para(5) = i     ; para(7) = j
          para(2) = spin1 ; para(4) = spin2 ; para(6) = spin2 ; para(8) = spin1
          call AppendOPT(self,7,para,V)
       enddo            ; enddo
      endsubroutine

      SUBROUTINE App_GeneralInter1(self,V,para)
        implicit none!V*c^+_{ispini} * c^+_{jspinj} * c_{kspink} * c_{lspinl}
        class(Ham),intent(inout)::self
        complex*16,intent(in)::V
        INTEGER,INTENT(IN)::PARA(8)
        !----------------------------------
        call AppendOPT(self,7,para,V)
      ENDSUBROUTINE

      SUBROUTINE App_GeneralInter2(self,V,para)
        implicit none!V*c^+_{ispini} * c_{jspinj} * c^+_{kspink} * c_{lspinl}
        class(Ham),intent(inout)::self
        complex*16,intent(in)::V
        INTEGER,INTENT(IN)::PARA(8)
        !----------------------------------
        call AppendOPT(self,6,para,V)
      ENDSUBROUTINE



    integer function GetEigenId(self)
      implicit none
      class(Ham),intent(inout)::self
      !-------------------------------
      if (.not.Is_Usable(self))then
         write(self%print,*)"WARNNING: Trying to get an eigen id of a Hamiltonian while it is not finished settting."
         write(self%print,*)"In principle, this is an error, only except you are sure what you are doing."
      endif
      GetEigenId = self%EigenId  !;write(*,*)self%EigenId
    endfunction


    logical function Is_Usable(self)
      implicit none
      class(Ham),intent(inout)::self
      !-------------------------------
      if (self%state.eq.3)then
        Is_Usable = .true.
      else
        Is_Usable = .false.
      endif
    endfunction

    integer function GetOptN(self)
      implicit none
      class(Ham),intent(inout)::self
      !-------------------------------
      if ( Is_Usable(self) )then
        GetOptN = self%optl%GetLen()
      else
        write(self%print,*)"ERROR: GetOptN is called while setting of H is not finished yet";stop
      endif
    endfunction



    function GetOptact(self,i) result(r)
      implicit none
      class(Ham),intent(inout)::self
      integer,intent(in)::i
      class(FermOper),pointer::r
      !-------------------------------
      class(odata),pointer::p
      p => self%optl%GetDataPointer(i)
      r => p%opt
    endfunction

    complex*16 function GetOptV(self,i)
      implicit none
      class(Ham),intent(inout)::self
      integer,intent(in)::i
      !-------------------------------
      class(odata),pointer::p
      p => self%optl%GetDataPointer(i)
      GetOptV = p%v
    endfunction




    subroutine Copy(a,b)
      implicit none
      class(Ham),intent(out)::a
      class(Ham),intent(in )::b
      !-------------------------------
      a%initiated = b%initiated
      a%ns        = b%ns
      a%optl      = b%optl
      a%l         = b%l
      a%EigenId   = b%EigenId
      a%state     = b%state
      a%IntZero   = b%IntZero
      a%print     = b%print

    endsubroutine

endmodule
