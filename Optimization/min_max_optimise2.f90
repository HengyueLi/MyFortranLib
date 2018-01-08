
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : Msearchsaddle
! OBJECT: TYPE(searchsaddle)
! USED  : Mminmaxoptimise
! DATE  : 2017-12-01
! AUTHOR: Hengyue Li
!--------------
! DESCRIPTION:
!            Find a searchsaddle point of system. The spirit is to find out minimum point and maximum point alternativily.
!
! STANDARD:
!              [sub] *call self%Initialization(mode,order,nmax,nmin,Maxitra,showdetails,deltaX,prX,prY,prG,PRE)
!                     order = 2n/2n +1 represent find minimum/maximum
!              [sub] *call self%Searching
!
!              Notice to override the functions below
!              !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!                  procedure,pass:: Func  => OverrideFun
!                  procedure,pass:: setx  => OverrideSet
!                  procedure,pass:: getx  => OverrideGet
!                  !-----------------------------------------
!                  procedure,pass:: IsRational =>OverrideIsRational
!              !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  [sub] a
! avalable gets:
!
!                  [fun] get_y()
!                        real(8)::get_y    output the final value of the stationary point.
!                  [fun] get_ierror()
!
! other :
!                  [sub] d
!
!
!
!
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module Msearchsaddle
    use Mminmaxoptimise
    implicit none





    type,extends(minmaxoptimise),private::findmax
        class(searchsaddle),pointer::SS    => null()  ! to searchsaddle
    contains
        procedure,pass::Func  => OverrideFunmax
        procedure,pass::get_x => OverrideGetmax
        procedure,pass::set_x => OverrideSetmax
        procedure,pass::Is_rational  => OverrideIsRationalmax
    end type

    type,extends(minmaxoptimise),private::findmin
        class(searchsaddle),pointer::SS    => null()  ! to searchsaddle
    contains
        procedure,pass::Func  => OverrideFunmin
        procedure,pass::get_x => OverrideGetmin
        procedure,pass::set_x => OverrideSetmin
        procedure,pass::Is_rational  => OverrideIsRationalmin
    end type



    TYPE::searchsaddle

        logical,private::initiated = .false.
        integer::order
        real(8)::PRE
        integer::nmax
        integer::nmin
        integer::n    ! n = nmin + nmax
        TYPE(findmax),private::fmax
        TYPE(findmin),private::fmin
        !type(ax)::something
        real(8)::y
        integer,private::ierr = -1
        !------------------------
        integer::print = 6
    contains
        procedure,pass::Initialization
        procedure,pass::Unitialization
        final::Finalization

        procedure,pass::searching
        procedure,pass::get_ierror
!        procedure,pass::get_y
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        procedure,pass:: Func  => OverrideFun
        procedure,pass:: setx  => OverrideSet
        procedure,pass:: getx  => OverrideGet
        !-----------------------------------------
        procedure,pass:: IsRational =>OverrideIsRational
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    END TYPE







    private::Initialization,Unitialization,Finalization
    private::searching
    private::OverrideFun   ,OverrideSet   ,OverrideGet   ,OverrideIsRational
    private::OverrideFunmax,OverrideSetmax,OverrideGetmax,OverrideIsRationalmax
    private::OverrideFunmin,OverrideSetmin,OverrideGetmin,OverrideIsRationalmin
    private::optimise_one_time
    private::get_xmax_xmin_to_x,set_x_to_xmax_xmin,set_Initiate_Y
    private::get_ierror,check_m9053


    contains

    !------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@  to be ovverided @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    real(8) function OverrideFun(self)
            implicit none
            class(searchsaddle),intent(inout)::self
            !------------------------
            ! For a given self%ax, calculate the function value
            !       OverrideFun = f(self%ax)
            !------------------------
            OverrideFun = 1._8  ;  stop "OverrideFun is not defined yet^"
    end function

    subroutine OverrideSet(self,swith,xmax,xmin) ! swith = "min"/"max"
            implicit none
            class(searchsaddle),intent(inout)::self
            character(3),intent(in)::swith
            real(8),intent(inout)::xmax(self%nmax)
            real(8),intent(inout)::xmin(self%nmin)
            !------------------------
            ! For a given self%ax, set x
            xmax = 0._8
            xmin = 0._8
            stop "OverrideSet is not defined in searchsaddle"
    end subroutine

    subroutine OverrideGet(self,swith,xmax,xmin) ! swith = "min"/"max"
            implicit none
            class(searchsaddle),intent(inout)::self
            character(3),intent(in)::swith
            real(8),intent(inout)::xmax(self%nmax)
            real(8),intent(inout)::xmin(self%nmin)
            !------------------------
            ! For a given self%ax, get x
            xmax = 0._8
            xmin = 0._8
            stop "OverrideGet is not defined in searchsaddle"
    end subroutine

    logical function OverrideIsRational(self)
            implicit none
            class(searchsaddle),intent(inout)::self
            !--------------------------------------------
            OverrideIsRational = .true.
    end function
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@







    subroutine Initialization(self,order,mode,nmax,nmin,Maxitra,showdetails,deltaX,prX,prY,prG,PRE,print_)
               implicit none
               class(searchsaddle),intent(inout),target::self
               integer::order,mode,nmax,nmin,Maxitra,showdetails
               real(8),intent(in)::deltaX,prX,prY,prG,PRE
               integer,intent(in),optional::print_
               !----------------------------------------------
               call Unitialization(self)                !;write(*,*)showdetails,6666
               if (present(print_)) self%print = print_
               self%initiated = .true.
               self%order     = order
               CALL self%fmax%Initialization(nmax,mode, .true. ,Maxitra,showdetails,deltaX,prX,prY,prG,self%print )
               CALL self%fmin%Initialization(nmin,mode, .false.,Maxitra,showdetails,deltaX,prX,prY,prG,self%print )
               self%fmax%SS => self
               self%fmin%SS => self
               self%PRE  = PRE
               self%nmax = nmax
               self%nmin = nmin
               self%n    = nmax + nmin
               if ( nmax*nmin .eq.0 )then
                   if (prx.ge.PRE)then
                       write(*,*)"if we only consider one case, prx >= PRE is not allowed." ;stop
                   end if
               end if
    end subroutine

    subroutine Unitialization(self)
               implicit none
               class(searchsaddle),intent(inout)::self
               !---------------------------------------
               if (self%initiated)then
                  self%fmax%SS     => null()
                  call self%fmax%Unitialization()
                  self%fmin%SS     => null()
                  call self%fmin%unitialization()
               end if
    end subroutine

    subroutine Finalization(self)
               implicit none
               type(searchsaddle),intent(inout)::self
               !---------------------------------------
               call Unitialization(self)
    end subroutine



    real(8) function OverrideFunmax(self)
            implicit none
            class(findmax),intent(inout)::self
            !------------------------
            OverrideFunmax = self%SS%Func()
    end function

    subroutine OverrideGetmax(self,x)
            implicit none
            class(findmax),intent(inout)::self
            real(8),intent(inout)::x(self%n)
            !------------------------
            real(8)::xmin( self%SS%nmin )
            call self%SS%getx("max",x,xmin)
    end subroutine

    subroutine OverrideSetmax(self,x)
            implicit none
            class(findmax),intent(inout)::self
            real(8),intent(inout)::x(self%n)
            !------------------------
            real(8)::xmin( self%SS%nmin )
            call self%SS%setx("max",x,xmin)
    end subroutine

    logical function OverrideIsRationalmax(self)
            implicit none
            class(findmax),intent(inout)::self
            !------------------------
            OverrideIsRationalmax = self%SS%IsRational()
    end function


    real(8) function OverrideFunmin(self)
            implicit none
            class(findmin),intent(inout)::self
            !------------------------
            OverrideFunmin = self%SS%Func()
    end function

    subroutine OverrideGetmin(self,x)
            implicit none
            class(findmin),intent(inout)::self
            real(8),intent(inout)::x(self%n)
            !------------------------
            real(8)::xmax( self%SS%nmax )
            call self%SS%getx("min",xmax,x)
    end subroutine

    subroutine OverrideSetmin(self,x)
            implicit none
            class(findmin),intent(inout)::self
            real(8),intent(inout)::x(self%n)
            !------------------------
            real(8)::xmax( self%SS%nmax )
            call self%SS%setx("min",xmax,x)
    end subroutine

    logical function OverrideIsRationalmin(self)
            implicit none
            class(findmin),intent(inout)::self
            !------------------------
            OverrideIsRationalmin = self%SS%IsRational()
    end function

!    subroutine searching(self)
!            implicit none
!            class(searchsaddle),intent(inout)::self
!            !-----------------------------------------------
!            real(8)::x1(self%n),x2(self%n),deltaX
!            integer::counting
!            logical::converged
!
!
!            counting  = mod(self%order,2)
!
!            call get_xmax_xmin_to_x(self,x1)
!            converged = .false.
!            Do while (.not.converged)
!                call optimise_one_time(self,counting)
!                call get_xmax_xmin_to_x(self,x2)
!                deltaX = dsqrt ( sum(  (x1-x2)**2 ) )
!                x1 = x2
!                counting  = counting + 1
!                converged = deltaX.le.self%PRE
!                if (.not.converged)then
!                    call set_Initiate_Y(self,counting)
!                end if
!            End Do
!    end subroutine

    subroutine searching(self)
            implicit none
            class(searchsaddle),intent(inout)::self
            !-----------------------------------------------
            real(8)::x1(self%n),x2(self%n),deltaX
            integer::counting
            logical::converged,maxfirst

            !-----------
            self%ierr = -1
            !-----------


            if (mod(self%order,2).eq.0)then
                maxfirst = .false.
            else
                maxfirst = .true.
            end if

            counting  = 0!

            call get_xmax_xmin_to_x(self,x1)
            converged = .false.
            Do while (.not.converged)
                if (maxfirst)then
                    call self%fmax%searching()
                    if (  check_m9053(self,self%fmax%get_ierr() )  ) goto 999
                    self%y = self%fmax%get_y()
                    call self%fmin%set_initiate_y(self%y)
                    call self%fmin%searching()
                    if (  check_m9053(self,self%fmin%get_ierr() )  ) goto 999
                    self%y = self%fmin%get_y()
                else
                    call self%fmin%searching()
                    if (  check_m9053(self,self%fmin%get_ierr() )  ) goto 999
                    self%y = self%fmin%get_y()
                    call self%fmax%set_initiate_y(self%y)
                    call self%fmax%searching()
                    if (  check_m9053(self,self%fmax%get_ierr() )  ) goto 999
                    self%y = self%fmax%get_y()
                end if
                call get_xmax_xmin_to_x(self,x2)
                deltaX = dsqrt ( sum(  (x1-x2)**2 ) )
                x1 = x2
                counting  = counting + 1
                converged = deltaX.le.self%PRE
                if (.not.converged)then
                    if (maxfirst)then
                       call self%fmax%set_initiate_y(self%y)
                    else
                       call self%fmin%set_initiate_y(self%y)
                    end if
                end if
            End Do

!            write(*,*)"---------------converged,x="
!            write(*,*)x1
!            write(*,*)"---------------------------"

       999 continue
    end subroutine

    logical function check_m9053(self,erro)
            implicit none
            class(searchsaddle),intent(inout)::self
            integer::erro
            !-----------------------------------------------
            if (erro.eq.-9053)then
                self%ierr = -9053
                check_m9053 = .true.
            else
                check_m9053 = .false.
            end if
    end function




    ! if maxormin =2n   , we searching min point
    ! if maxormin =2n+1 , we searching max point
    subroutine optimise_one_time(self,maxormin)
            implicit none
            class(searchsaddle),intent(inout)::self
            integer,intent(in)::maxormin
            !-----------------------------------------------
            if (  mod(maxormin,2).eq.0  )then
                call self%fmin%searching()
                self%y = self%fmin%get_y()
            else
                call self%fmax%searching()
                self%y = self%fmax%get_y()
            end if
    end subroutine

    subroutine set_Initiate_Y(self,maxormin)
            implicit none
            class(searchsaddle),intent(inout)::self
            integer,intent(in)::maxormin
            !-----------------------------------------------
            if (  mod(maxormin,2).eq.0  )then
                 call self%fmin%set_initiate_y(self%y)
            else
                 call self%fmax%set_initiate_y(self%y)
            end if
    end subroutine

    subroutine get_xmax_xmin_to_x(self,x)
            implicit none
            class(searchsaddle),intent(inout)::self
            real(8),intent(inout)::x(self%n)
            !-----------------------------------------------
            call self%getx( "max", x(1:self%nmax),x(self%nmax+1:self%n) )
            call self%getx( "min", x(1:self%nmax),x(self%nmax+1:self%n) )
    end subroutine

    subroutine set_x_to_xmax_xmin(self,x)
            implicit none
            class(searchsaddle),intent(inout)::self
            real(8),intent(inout)::x(self%n)
            !-----------------------------------------------
            call self%setx( "max", x(1:self%nmax),x(self%nmax+1:self%n) )
            call self%setx( "min", x(1:self%nmax),x(self%nmax+1:self%n) )
    end subroutine


    integer function get_ierror(self)
            implicit none
            class(searchsaddle),intent(inout)::self
            !-----------------------------------------------
            get_ierror = self%ierr
    end function


end module
