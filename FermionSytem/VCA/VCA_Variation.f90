
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : VCA_Variation
! OBJECT: TYPE(VCAva)
! USED  : CodeObject , VCA_DeltaH , VCA_WaldFun , Msearchsaddle , functionalsubs
! DATE  : 2018-01-07
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            Use some nummerical method to search the stationary point of the system.
!            We follow the cross over method.
!              https://arxiv.org/abs/0806.0266 "Correlated band structure of NiO, CoO,..."
!
!           CM folw:
!           -----------
!           definition :   A = \alpha  ,
!
!           input: dA0 , Am
!
!                         ┌─────────────────┐
!                         │initiate: A0 = 0 │
!                         └────────┬────────┘
!                                  ↓
!                           ┌──────┴───────┐
!                           │  dA = dA0    │
!                           └──────┬───────┘
!            ┌────────────────────→↓ ←────────────────────────┐
!            │              ┌──────┴───────┐                  │
!            │              │ A = A0 + dA  │                  │
!            │              └──────┬───────┘                  │
!            │                     ↓                          │
!            │                ┌────┴──────┐                   │
!            │                │Last = A>=1│                   │
!            │                │ if Last ? │                   │
!            │                └─┬─────┬───┘                   │
!            │               yes↓     ↓No                     │
!            │                ┌─┴─┐   │                       │
!            │                │A=1│   │                       │
!            │                └─┬─┘   │                       │
!            │                  ↓     ↓                       ↑
!            │          ┌───────┴─────┴────────┐              │
!            │          │ backup dH            │              │
!            │          │ Searching successful?│              │
!            ↑          └──────┬─────────────┬─┘              │
!            │             yes ↓             ↓ No┌───────────┐│
!            │          ┌──────┴────┐        └───┤Recover dH ││
!            │          │ if Last ? │            │Reduce  dA ││
!            │          └──┬─────┬──┘ yes        └──────────┬┘│
!            │          NO ↓     └────→─────┐           ┌───┴─┴┐
!            │          ┌──┴───────┐     ┏━━┷━━━━━━━━━┓ │dA<Am?│
!            │          │A0=A      │     ┃out,ierr = 0┃ └┬─────┘
!            │          │Recover dA│     ┗━━━━━━━━━━━━┛  ↓yes
!            │          └─┬────────┘                  ┏━━┷━━━━━━━━━┓
!            └────────────┘                           ┃out,ierr = 1┃
!                                                     ┗━━━━━━━━━━━━┛
!
! STANDARD:
!            *CALL Initialization( )
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                   [sub] Initialization(dh,wf,jobi,jobr,print_,show_  )
!                          class(VCAdH)      ,intent(inout),target :: dh
!                          class(VCA_WaldFun),intent(inout),target :: wf
!                          integer           ,intent(in)           :: jobi(:)
!                          real*8            ,intent(in)           :: jobr(:)
!                          integer           ,intent(in) ,optional :: print_,show_
!
!                            deltaX= jobr(1)
!                            prX   = jobr(2)
!                            prY   = jobr(3)
!                            prG   = jobr(4)
!                            PRE   = jobr(5)
!
!                           dAlpha = jobr(6)    ! step of alpha
!                           Alpham = jobr(7)    ! tolarence of small alpha
!
! avalable gets:
!                   [fun] G
!
! avalable is :
!                  ![fun] i
! others      :
!                  [fun] OptimisePhasePoint()
!                        Start from recent dH, searching a new stationary point dH'.
!                        The value of dH is changed.
!                        return error code.
!
!                  [fun] CrossOverSearching(ResetDH)
!                        logical::ResetDH
!                        if ResetDH : dH will be set to zero first. If not, dH will be kept.
!                        Using CM method to find the stationary point.
!                        return the error code.
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module VCA_Variation
  use CodeObject
  use VCA_DeltaH
  use VCA_WaldFun
  use Msearchsaddle
  use functionalsubs
  implicit none



  type,extends(searchsaddle),private::Vca_saddle
        class(VCAva)        ,pointer::Varant
  contains
       !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       procedure,pass:: Func  => OverrideFun
       procedure,pass:: setx  => OverrideSet
       procedure,pass:: getx  => OverrideGet
       !-----------------------------------------
       procedure,pass:: IsRational =>OverrideIsRational
       !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  end type








  integer,parameter,private::NmaxJob = 50
  type,extends(object)::VCAva
    private
    class(VCAdH),pointer :: dH     => null()
    class(waldf),pointer :: wf     => null()
    type(Vca_saddle)     :: sd

    integer::jobi(NmaxJob)
    real*8 ::jobr(NmaxJob)
  contains
    procedure,pass::Initialization
    procedure,pass::OptimisePhasePoint
    procedure,pass::CrossOverSearching
  endtype


  private::Initialization
  private::OptimisePhasePoint
  private::CrossOverSearching



  !!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !!  Override part
  private::OverrideFun,OverrideIsRational,OverrideSet,OverrideGet
contains

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
real(8) function OverrideFun(self)
        implicit none
        class(Vca_saddle),intent(inout)::self
        !------------------------
        OverrideFun = self%Varant%wf%GetLatticeOmegaPerSite()
end function

logical function OverrideIsRational(self)
        implicit none
        class(Vca_saddle),intent(inout)::self
        !--------------------------------------------
        integer:: wtp,jc
        real*8 :: Large,V
        !-------------------
          Large = 100._8
        !-------------------
        wtp = self%Varant%getprint()
        OverrideIsRational = .true.
        do jc = 1 , self%Varant%dH%GetNumOfVarialtinalTerms()
           V = self%Varant%dH%GetValueByIndex(jc)
           if (  abs(v) >Large )then
              write(wtp,*)"+++++++++++++++++++++++++++++++++++++++++++++++++"
              write(wtp,*)"variational term:",trim(adjustl(self%Varant%dH%GetDisciption(jc)))&
                           ,"is too large"
              write(wtp,*)"+++++++++++++++++++++++++++++++++++++++++++++++++"
              OverrideIsRational = .false.
           endif
        enddo
end function

subroutine OverrideSet(self,swith,xmax,xmin) ! swith = "min"/"max"
        implicit none
        class(Vca_saddle),intent(inout)::self
        character(3),intent(in)::swith
        real(8)     ,intent(inout)::xmax(self%nmax)
        real(8)     ,intent(inout)::xmin(self%nmin)
        !------------------------
        TYPE(funcsubs)::f
        integer::jc,p
        character(3)::upper
        upper = f%get_string_upper(swith)
        DO JC = 1 , self%Varant%dH%GetTotNbyHow(swith)
           p = self%Varant%dH%GetHowId(swith,jc)
           select case(upper)
           case("MAX")
             call self%Varant%dH%SetValueByIndex(p,xmax(jc))
           case("MIN")
             call self%Varant%dH%SetValueByIndex(p,xmin(jc))
           case default
             write(self%Varant%getprint(),*)"ERROR: 'SWITH' type error in OverrideSet@VCA_Variation"
             stop
           endselect
        enddo
end subroutine



subroutine OverrideGet(self,swith,xmax,xmin) ! swith = "min"/"max"
        implicit none
        class(Vca_saddle),intent(inout)::self
        character(3),intent(in)::swith
        real(8),intent(inout)::xmax(self%nmax)
        real(8),intent(inout)::xmin(self%nmin)
        !------------------------
        TYPE(funcsubs)::f
        integer::jc,p
        character(3)::upper
        upper = f%get_string_upper(swith)
        DO JC = 1 , self%Varant%dH%GetTotNbyHow(swith)
           p = self%Varant%dH%GetHowId(swith,jc)
           select case(upper)
           case("MAX")
             xmax(jc) = self%Varant%dH%GetValueByIndex(p)
           case("MIN")
             xmin(jc) = self%Varant%dH%GetValueByIndex(p)
           case default
             write(self%Varant%getprint(),*)"ERROR: 'SWITH' type error in OverrideSet@VCA_Variation"
             stop
           endselect
        enddo
end subroutine


!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@





  subroutine Initialization(self,dh,wf,jobi,jobr,print_,show_ )
    implicit none
    class(VCAva),intent(out),target   :: self
    class(VCAdH),intent(inout),target :: dh
    class(waldf),intent(inout),target :: wf
    integer     ,intent(in)           :: jobi(:)
    real*8      ,intent(in)           :: jobr(:)
    integer     ,intent(in) ,optional :: print_,show_
    !-----------------------------------------------------------
    call self%SetInitiated(.true.)
    if (present(print_)) call self%setprint(print_)
    if (present(show_ )) call self%setshow(show_  )

    self%dh => dh
    self%wf => wf

    if (  max(  size(jobr) ,size(jobi) )> NmaxJob )then
      write(self%getprint(),*)"Over: NmaxJob=",NmaxJob," in VCA_Variation is too small"
      write(self%getprint(),*)"set it to :",max(  size(jobr) ,size(jobi) )
      stop
    endif

    self%jobi(1:size(jobi)) = jobi
    self%jobr(1:size(jobr)) = jobr



    !-------------------------------------------
    ! initiate sd
    self%sd%Varant => self
    call self%sd%Initialization(order=self%jobi(1) , mode= self%jobi(2) ,&
                  nmax= self%dh%GetTotN_Var_max(),nmin=self%dh%GetTotN_Var_min(),&
                  Maxitra = 2000,showdetails=self%getshow(),&
                  deltaX= self%jobr(1),&
                  prX   = self%jobr(2),&
                  prY   = self%jobr(3),&
                  prG   = self%jobr(4),&
                  PRE   = self%jobr(5)             )

  endsubroutine


   integer function OptimisePhasePoint(self)
     implicit none
     class(VCAva),intent(inout):: self
     !---------------------------------
     call self%sd%Searching()
     OptimisePhasePoint = self%sd%get_ierror()
  endfunction





  integer function CrossOverSearching(self,resetdh)
    implicit none
    class(VCAva),intent(inout)  :: self
    logical,intent(in)::resetdh
    !-------------------------------------
    real*8::A0,A,dA,Am
    logical::last
    integer::ierro

    if (resetdh ) call self%dH%ReSetToZero()

    A0 = 0._8   ;  dA = self%jobr(6)   ; Am =  self%jobr(7)

100 continue

    A = A0 + dA

    last = A>=1._8

    if (last) A = 1._8


    Call self%dH%BackUp()
    call self%dH%setalpha(A)
    ierro = OptimisePhasePoint(self)
    !-----------------

    if ( ierro ==0 ) then
      ! variation is successful!
      write(self%getprint(),*)"  Stationary point have been found successfully!."
      !  report stationary point
      Call self%dH%report(self%getprint())
      if (last) then
        CrossOverSearching = 0
        goto 999
      else
        A0 = A
        Call RecoverDA(dA,self%jobr(6))
        goto 100
      endif
    else
      write(self%getprint(),*)"=================================================="
      write(self%getprint(),*)"Step maybe too large,change it from",dA
      !------------------------------
      call self%dH%Recover()
      call ReduceDA(dA,self%jobr(6))
      !------------------------------
      write(self%getprint(),*)"to ",dA
      write(self%getprint(),*)"=================================================="



      if (dA <= Am)then
        CrossOverSearching = 1
        goto 999
      endif

      goto 100

    endif


999 continue

contains
  subroutine ReduceDA(dA,dA0)
    implicit none
    real*8,intent(inout)::dA
    real*8,intent(in)   ::dA0
    !-------------------------
    dA = dA / 2
  endsubroutine

  subroutine RecoverDA(dA,dA0)
    implicit none
    real*8,intent(inout)::dA
    real*8,intent(in)   ::dA0
    !-------------------------
    dA = dA0
  endsubroutine

  endfunction


!call self%sd%self%Searching()



endmodule
