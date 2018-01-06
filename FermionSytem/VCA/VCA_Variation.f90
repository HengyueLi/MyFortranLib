
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : VCA_Variation
! OBJECT: TYPE(VCAva)
! USED  : CodeObject , VCA_DeltaH , VCA_WaldFun , Msearchsaddle , functionalsubs
! DATE  : 2017-12-30
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            Use some nummerical method to search the stationary point of the system.
!            We follow the cross over method.
!              https://arxiv.org/abs/0806.0266 "Correlated band structure of NiO, CoO,..."
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
!
! avalable gets:
!                   [fun] G
!
! avalable is :
!                  ![fun] i
! others      :
!                  ![sub] T
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
  endtype


  private::Initialization



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

  subroutine CrossOverMethod(self)
    implicit none
    class(VCAva),intent(inout)  :: self
    !-------------------------------------

  endsubroutine

!call self%sd%self%Searching()



endmodule
