






!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : Mminmaxoptimise
! OBJECT: TYPE(minmaxoptimise)
! USED  : basic_math_functions
! DATE  : 2017-11-5
! AUTHOR: Hengyue Li
!--------------
! DESCRIPTION:
!            Find a maxmer or minmer point of a function.
!
! STANDARD:
!              [sub] *CALL self%Initialization(n,mode,Maxmode,MItrat,showdetails,deltaX,prX,prY,prG)
!                    integer::n      is the dimension of the problem.
!
!                    integer::mode
!                             mode = 1 for quasi newton method (copy from Naomichi Sato )
!                             mode = 2 for quasi newton method (code from book)
!                             mode = 3 searching in dimension one by one (self built, slow, unstable)
!                             mode = 4 Downhill Simplex Method (working perfect, extremely slow)
!                                      in sub DownHill_Check_Converged()  we can change the convergecy criterion.
!                    logical::Maxmode to check if we consider a maxmer problem
!                    integer::MItrat   max interaction
!                    integer::showdetails
!                    class(aux)::contains the stationary point.
!              [SUB] *CALL self%searching()
!
!              Notice that 3 process should be overrided
!                 procedure,pass::Func  => OverrideFun
!                 procedure,pass::get_x => OverrideGet
!                 procedure,pass::set_x => OverrideSet
!              And we have a abalable override function
!                 procedure,pass::Is_rational  => OverrideIsRational
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  [sub] set_initiate_y(y) .
!                        before calling self%searching(), if we set initiate value of y, we can save time.
!                        The input y is the real value of function.(not depend on Max mode or not.)
! avalable gets:
!
!                  [fun] get_n()
!                        get the dimension of the system
!                  [fun] get_get_ierr()
!                        return -9053 if Is_rational return false
!                  [fun] get_y()
!                        real(8)
! avalable IS  :
!                  [fun] Is_rational()
!                        check the recent state is rational or not.
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


module Mminmaxoptimise
    implicit none




    !-------------book Numerical Recipes in Fortran77 ....  use -------
    integer,parameter,PRIVATE::BNMAX = 50
    INTEGER,PRIVATE::ncom
    REAL*8,PRIVATE::pcom(BNMAX),xicom(BNMAX)
    !------------------------------------------------------------------



    type::minmaxoptimise


          logical,private::Initiated       = .false.
          integer,private::mode            = -1
          logical,private::Maxmode         = .false.
          integer::N
          !━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
                                                                     !┃
          !---------- input/output   -------------------             !┃
!          class(*),pointer :: ax => null()                          !┃
          !---------- final position -------------------             !┃
          real(8),allocatable,private::xout(:)                       !┃
          !---------- final value    -------------------             !┃
          real(8),private::tempY     ! temporarily save y value      !┃
          !---------- final gradient -------------------             !┃
          real(8),allocatable,private::tempG(:)                      !┃
          !━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
          integer,private::MaxItration     = 10000
          integer,private::showdetails     = 0        ! input 0,1,2,... to show details
          real(8),private::deltaX          = 1.0e-8   ! Use to calculate the differential of the funcion.
          real(8),private::preG            = 1.0e-6   ! convergecy of G:| y'- yc'|<preG
          real(8),private::preY            = 1.0e-8   ! convergecy of Y:| y - yc |<preY
          real(8),private::preX            = 1.0e-6   ! convergecy of X:| x - xc |<preX
          integer::ierr
          logical,private::useInitiate_y   = .false.

          integer::wtp = 6

    contains
          procedure,pass::Initialization
          procedure,pass::Unitialization
          final         ::Finalization

          procedure,pass::searching
          procedure,pass::set_initiate_y
          procedure,pass::get_n
          procedure,pass::get_y
          procedure,pass::get_g
          procedure,pass::get_ierr
          !@@@@@@@@@@@ need to be ovveride  @@@@@@@@@
          procedure,pass::Func  => OverrideFun
          procedure,pass::get_x => OverrideGet
          procedure,pass::set_x => OverrideSet
          !-----------avalabel
          procedure,pass::Is_rational  => OverrideIsRational
          !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    end type

    private::Initialization,Unitialization,Finalization

    private::fval,gval,fval0,gval0,Vatiation,ReportStat,lsearc,BFGS_from_Sato
    private::searching,get_n
    private::OverrideFun,OverrideGet,OverrideSet,OverrideIsRational

    private::set_initiate_y
    private::get_ierr,get_y,get_g

    private::CGmethod_my_frame,minimize
    private::check_Is_min



    private::OneDimensionalLocalMin
    private::OneDimLoc_one_direct_onestop_for_all_dim
    private::OneDimLoc_one_direct_convergy
    private::OneDimLoc_one_direct_one_step


    private::DownHill,amoeba,amotry,DownHill_Check_Converged
    private::frprmn,linmin,f1dim,brent,mnbrak
    private::lnsrch,dfpmin



    contains

!@@@@@@@@@@@@@@@@@@  to be ovverided @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    real(8) function OverrideFun(self)
            implicit none
            class(minmaxoptimise),intent(inout)::self
            !------------------------
            ! For a given self%ax, calculate the function value
            !       OverrideFun = f(self%ax)
            !------------------------
            OverrideFun = 1._8  ;  stop "OverrideFun is not defined yet^"
    end function

    subroutine OverrideGet(self,x)
            implicit none
            class(minmaxoptimise),intent(inout)::self
            real(8),intent(inout)::x(self%n)
            !------------------------
            !  for a given self%ax,  abstract the posistion x
            !       x = f(self%ax)
            x = 1._8            ;  stop "OverrideGet is not defined yet^"
    end subroutine

    subroutine OverrideSet(self,x)
            implicit none
            class(minmaxoptimise),intent(inout)::self
            real(8),intent(inout)::x(self%n)
            !------------------------
            !  for a given x,  set x to ax
            !       ax = f(x)
                                ;  stop "OverrideSet is not defined yet^"
    end subroutine

    logical function OverrideIsRational(self)
            implicit none
            class(minmaxoptimise),intent(inout)::self
            !------------------------
            OverrideIsRational = .true.
    end function
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@






    subroutine Initialization(self,N,mode,Maxmode,MItrat,showdetails,deltaX,prX,prY,prG,print_)
               implicit none
               class(minmaxoptimise),intent(inout)::self
               integer::N
               integer,intent(in)::mode
               logical,intent(in)::Maxmode
               integer,intent(in)::MItrat
               integer,intent(in)::showdetails
               real(8),intent(in)::deltaX,prX,prY,prG
               integer,intent(in),optional::print_

               !----------------------------------------
               call Unitialization(self)
               if (present(print_)) self%wtp = print_

               self%N           = N
               self%Initiated   = .true.
               self%mode        = mode
               self%Maxmode     = Maxmode
               self%MaxItration = MItrat
               self%showdetails = showdetails
               self%deltaX      = deltaX
               self%preX        = prX
               self%preY        = prY
               self%preG        = prG

               !---------------------------
               allocate( self%xout ( self%N )   )
               allocate( self%tempG( self%N )   )
    end subroutine


    subroutine Unitialization(self)
               implicit none
               class(minmaxoptimise),intent(inout)::self
               !-------------------------------------------
               if ( self%Initiated ) then
                    self%Initiated = .false.
                    deallocate(self%xout )
                    deallocate(self%tempG)
               end if
    end subroutine

    subroutine Finalization(self)
               implicit none
               type(minmaxoptimise),intent(inout)::self
               !-------------------------------------------
               call Unitialization(self)
    end subroutine


    integer function get_ierr(self)
               implicit none
               class(minmaxoptimise),intent(inout)::self
               !-------------------------------------------
               get_ierr = self%ierr
    end function

    integer function get_n(self)
               implicit none
               class(minmaxoptimise),intent(inout)::self
               !-------------------------------------------
               get_n = self%N
    end function

    subroutine get_G(self,G)
               class(minmaxoptimise),intent(inout)::self
               real(8),intent(out)::G(self%N)
               !-------------------------------------------
               if (self%Maxmode)then
                   G = -self%tempG
               else
                   G =  self%tempG
               end if
    end subroutine


    real(8) function get_y(self)
               implicit none
               class(minmaxoptimise),intent(inout)::self
               !-------------------------------------------
               !get_y = self%yout
               if (self%Maxmode)then
                   get_y = -self%tempY
               else
                   get_y =  self%tempY
               end if
    end function

    subroutine set_initiate_y(self,y)
               implicit none
               class(minmaxoptimise),intent(inout)::self
               real(8),intent(in)::y
               !-------------------------------------------
               !self%yout = y
               if (self%Maxmode)then
                   self%tempY = -y
               else
                   self%tempY =  y
               end if
               self%useInitiate_y = .true.
    end subroutine



    !===============================================================================
      subroutine fval(self, X, OBJFUN)
         implicit none
    !*******************************************************************************
    !  get function value
    !  N : dimension
    !  X : position
    !  OBJFUN : return value
    !*******************************************************************************
         class(minmaxoptimise),intent(inout)::self
         real(8) X(self%N), OBJFUN
         call self%set_x(x)

         OBJFUN = self%Func()

         IF ( self%Maxmode  )THEN
            OBJFUN=-OBJFUN
         END IF
         !------------------------------------------------------!
         !       checking the recent state                      !
         if (  .not. self%Is_rational() )  self%ierr = -9053    !
         !------------------------------------------------------!
         end subroutine fval

      subroutine fval0(self, X, OBJFUN)
         implicit none
         class(minmaxoptimise),intent(inout)::self
         real(8) X(self%N), OBJFUN
         OBJFUN =  self%tempY!self%yout
         end subroutine fval0
    !===============================================================================


    !===============================================================================
      subroutine gval(self, X, OBJGRA)
        implicit none
        class(minmaxoptimise),intent(inout)::self
        real(8),INTENT(IN)::X(self%n)
        REAL(8)::OBJGRA(self%N)
        real(8)::TEMP1,TEMP2,XT(self%N)
        INTEGER::JC

        if (self%ierr.eq.-9053) goto 999

        CALL fval(self, X, TEMP1)

        if (self%ierr.eq.-9053) goto 999

        DO JC=1,self%N
            XT         = X
            XT(JC)     = XT(JC) + self%deltaX
            CALL fval(self, XT, TEMP2)
            OBJGRA(JC) =(TEMP2-TEMP1)/self%deltaX
        END DO

    999 continue
        end subroutine gval

      subroutine gval0(self, X, OBJGRA)
        implicit none
        class(minmaxoptimise),intent(inout)::self
        real(8),INTENT(IN)::X(self%n)
        REAL(8)::OBJGRA(self%N)
        real(8)::TEMP1,TEMP2,XT(self%N)
        INTEGER::JC

!        TEMP1 =  self%yout
        TEMP1 = self%tempY

        DO JC=1,self%N
            XT         = X
            XT(JC)     = XT(JC) + self%deltaX
            CALL fval(self, XT, TEMP2)
            if (self%ierr.eq.-9053) goto 999
            OBJGRA(JC) =(TEMP2-TEMP1)/self%deltaX
        END DO
    999 continue
        end subroutine gval0

    !===============================================================================

     subroutine searching(self)
                implicit none
                class(minmaxoptimise),intent(inout)::self
                !-----------------------------------------
                if (self%n.eq.0) goto 999
                !----------------
                self%ierr = -1
                !----------------

                call self%get_x(self%xout)                              !;self%useInitiate_y = .false.

                !----------------------------------------------------
                if (self%useInitiate_y)then
                    call fval0(self, self%xout, self%tempY)
                    call gval0(self, self%xout, self%tempG)
                else
                    call fval(self, self%xout, self%tempY)
                    call gval(self, self%xout, self%tempG)
                end if
                self%useInitiate_y = .false.

                if (self%ierr.eq.-9053) goto 999
                !----------------------------------------------------
                call Vatiation(self, self%xout, self%tempY, self%tempG, self%ierr)
                call self%set_x(self%xout)
            999 continue
     end subroutine

    subroutine Vatiation(self, X, OBJFUN, OBJGRA, IND)
               implicit none
               class(minmaxoptimise),intent(inout)::self
               real(8),intent(inout)::x(self%N)
               real(8),intent(inout)::OBJFUN
               real(8),intent(inout)::OBJGRA(self%N)
               integer,intent(inout)::ind
               !----------------------------------------

               select case(self%mode)
                      case(1)
                          call BFGS_from_Sato        (self, X, OBJFUN, OBJGRA, IND) !;write(*,*)8888,IND,check_Is_min(self,x,OBJFUN)
                      case(2)
                          call dfpmin                (self, X, OBJFUN, OBJGRA, IND)
                      ! case(2)
                      !     !call CGmethod_my_frame     (self, X, OBJFUN, OBJGRA, IND)
                      !     call frprmn                 (self, X, OBJFUN, OBJGRA, IND)
                      case(3)
                          call OneDimensionalLocalMin(self, X, OBJFUN, OBJGRA, IND)
                      case(4)
                          call DownHill              (self, X, OBJFUN, OBJGRA, IND)
                      case default
                          write(*,*)"unknow input mode in sub Vatiation";stop
               end select
    end subroutine



logical function check_Is_min(self,x,y)
       implicit none
       class(minmaxoptimise),intent(inout)::self
       real(8),intent(inout)::x(self%N)
       real(8),intent(in)::y
       !-------------------------------------
       real(8)::ym,yp,x1(self%n),x2(self%n)
       integer::jc
       check_Is_min = .false.

       do jc = 1, int(self%n)
          x1 = x
          x1(jc) = x1(jc) + self%deltaX
          call fval(SELF, x1, yp)
          x1 = x
          x1(jc) = x1(jc) - self%deltaX
          call fval(SELF, x1, ym)
          !--------------------
          if ( .NOT.( ym + yp .gt. 2*y  ) ) goto 999
          !--------------------
       enddo
       check_Is_min = .true.
  999  continue
end function




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!in the following we only to used this!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            << call gval(self, X, OBJGRA) >>
!            << call fval(self, X, OBJFUN) >>
!
!           class(minmaxoptimise),intent(inout)::self
!           real*8::x(self%n),objgra(self%n),objfun
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








!===============================================================================
    subroutine BFGS_from_Sato(self, X, OBJFUN, OBJGRA, IND)
!===============================================================================

    implicit none
    class(minmaxoptimise),intent(inout)::self

    integer i, j

!       ind= error
! input level= 0,1,2
! ITMAX  max interation
! level  show details
! IT
!*******************************************************************************
!*******************************************************************************

!-------------------------------------------------------------------------------
   ! INTEGER,PARAMETER::NMAX = 100000
    integer IND, LEVEL, IT, N, ITMAX
    real(8)::X(SELF%N), OBJFUN, OBJGRA(SELF%N), INVHES(SELF%N, SELF%N), DIREC(SELF%N), NEWFUN, OLDFUN,&
            NEWX(SELF%N), OLDX(SELF%N), OLDGRA(SELF%N), GNORM, SNORM, YNORM, S(SELF%N), Y(SELF%N), W(SELF%N),&
            SY, YHY, C1, C2, EPS1, EPS2, EPS3, BOUND, ETA

!-------------------------------------------------------------------------------
!    EPS1 = 1.0d-3  !|grandient|
!    EPS2 = 1.0d-6  !|delta y|
!    EPS3 = 1.0d-8  !|delta x|
    eps1  = self%preG
    eps2  = self%preY
    eps3  = self%preX
    N     = self%N
    LEVEL = self%showdetails
    ITMAX = self%MaxItration

    BOUND = 1.0d+20
    ETA = 1.0d-4
!-------------------------------------------------------------------------------

!-------set the initial inverse Hessian approximation to be the identity matrix

    do i = 1, N  !å˜ä½è¡Œåˆ—ç”Ÿæˆ
       INVHES(i, i) = 1.0d0
       do j = 1, i-1
          INVHES(i, j) = 0.0d0
          INVHES(j, i) = 0.0d0
       end do
    end do

    IND = 0
    IT = 0        !;write(*,*)LEVEL,6666
    if(LEVEL > 0) call ReportStat(self,self%showdetails, X, OBJFUN, OBJGRA, IT)

!-------iteration

!-------compute search direction

    3 do i = 1, N
         DIREC(i) = 0.0d0 !åˆæœŸåŒ–
         do j = 1, N
            DIREC(i) = DIREC(i)-INVHES(i, j)*OBJGRA(j)
         end do
      end do

!-------perform line search

      ! NEWXã¨NEWFUNã®å€¤ã‚’å¾—ã‚‹
      call lsearc(SELF,N, X, NEWX, OBJFUN, NEWFUN, OBJGRA, DIREC, EPS3, IND)

      if (self%ierr.eq.-9053) goto 999

      if(IND .eq. 5) return

!-------evaluate gradient at the new iterate

      OLDFUN = OBJFUN
      OBJFUN = NEWFUN

      do i = 1, N  ! åº§æ¨™ã®æ›´æ–°
         OLDX(i) = X(i)
         X(i) = NEWX(i)
         S(i) = X(i) - OLDX(i) ! æ–°ã—ã„åº§æ¨™ã¨å¤ã„åº§æ¨™ã®å·®
         OLDGRA(i) = OBJGRA(i)
      end do


         !call gval(N, X, OBJGRA)
         CALL gval(self,X,OBJGRA)

         if (self%ierr.eq.-9053) goto 999


      do i = 1, N
         Y(i) = OBJGRA(i)-OLDGRA(i) ! æ›´æ–°ã«ã‚ˆã‚‹å‹¾é…ãƒ™ã‚¯ãƒˆãƒ«ã®å·®
      end do

!-------check termination criteria.

      do i = 1, N
         if(abs(X(i)) > BOUND) then
            IND = 4
            return
         end if
      end do

      GNORM = 0.0d0
      do i = 1, N
         GNORM = GNORM + OBJGRA(i)*OBJGRA(i) !å‹¾é…ãƒ™ã‚¯ãƒˆãƒ«ã®å†…ç©
      end do
      GNORM = sqrt(GNORM)

      if (OLDFUN-OBJFUN < EPS2*(1.0+abs(OLDFUN)) .and. GNORM < EPS1) return

      if (IT > ITMAX) then
         if (GNORM < EPS1) then
            IND = 1
         else if (OLDFUN-OBJFUN < EPS2*(1.0+abs(OLDFUN))) then
            IND = 2
         else
            IND = 3
         end if

         return
      end if
      !**************** I put on this myself. It looks like a bug if we leak this.
!      author : hengyue li
      if ( dsqrt(sum(y**2)).le. eps1 ) return

!-------------------------------------------------------------------------------

!-------update the inverse Hessian approximation

      SNORM = 0.0
      YNORM = 0.0
      SY = 0.0

      do i = 1, N
         SNORM = SNORM + S(i)*S(i)
         YNORM = YNORM + Y(i)*Y(i)
         SY = SY + S(i)*Y(i)
      end do

      SNORM = sqrt(SNORM)
      YNORM = sqrt(YNORM)

      if (SY < ETA*SNORM*YNORM) then ! é€†ãƒ˜ã‚·ã‚¢ãƒ³è¡Œåˆ—ã®æ­£å®šå€¤æ€§ãŒæº€ãŸã•ã‚Œãªã„å ´åˆ
         if (LEVEL .eq. 2) then
            write(*, *) 'Hessian update is skipped at iteration ', IT
         end if
         goto 12 ! é€†ã¸ã‚·ã‚¢ãƒ³ã‚’æ›´æ–°ã›ãšå†åº¦é™ä¸‹æ–¹å‘ã¨ã‚¹ãƒƒãƒ†ãƒ—å¹…ã‚’æŽ¢ã‚‹(lsearch)
      end if

      YHY = 0.0d0
      do i = 1, N
         W(i) = 0.0d0
         do j = 1, N
            W(i) = W(i) + INVHES(i, j)*Y(j) ! Hy ã®è¨ˆç®—
         end do
         YHY = YHY + Y(i)*W(i) ! yHy ã®è¨ˆç®—
      end do

      C2 = 1.0d0 / SY          ! 1/sy ã®è¨ˆç®—
      C1 = (1.0d0 + C2*YHY)*C2 ! (1 + yHy/sy)*(1/sy) ã®è¨ˆç®—

      do i = 1, N
         do j = 1, i
            INVHES(i, j) = INVHES(i, j)+C1*S(i)*S(j)-C2*(S(i)*W(j)+W(i)*S(j))
            INVHES(j, i) = INVHES(i, j) ! å¯¾ç§°è¡Œåˆ—ï¼ˆæ­£å®šå€¤è¡Œåˆ—ï¼‰
         end do
      end do

  12  IT = IT + 1
      if (LEVEL > 0) call ReportStat(self,self%showdetails, X, OBJFUN, OBJGRA, IT)

      goto 3 ! åŽæŸã™ã‚‹ã‹ITMAXã¾ã§æœ€é©è§£Xã‚’æŽ¢ã™
999 continue
    end subroutine BFGS_from_Sato
!===============================================================================


!===============================================================================
    subroutine ReportStat(self,show, X, OBJFUN, OBJGRA, IT)
      implicit none
      class(minmaxoptimise),intent(inout)::self
      integer::show,IT,n,i
      real(8) X(self%N), OBJFUN, OBJGRA(self%N)
      n = self%n
      select case(show)
      case(0)
      case(1)
        call reportLevel1(self%wtp,self%Maxmode)
      case(2)
        call reportLevel1(self%wtp,self%Maxmode)
        call reportLevel2(self%wtp,self%Maxmode)
      case(3) ! only show x, y  No G
        call reportLevel1(self%wtp,self%Maxmode)
        call NonGshow(self%wtp,self%Maxmode)
      endselect

    contains

      subroutine reportLevel1(wtp,maxmode)
                 implicit none
                 integer,intent(in)::wtp
                 logical,intent(in)::maxmode
                 !-------------------------
                 write(wtp, *)'iteration = ', IT
                 if (maxmode)then
                     write(wtp, *)'object value = ', -OBJFUN
                 else
                     write(wtp, *)'object value = ',  OBJFUN
                 end if
      endsubroutine

      subroutine reportLevel2(wtp,maxmode)
                 implicit none
                 integer,intent(in)::wtp
                 logical,intent(in)::maxmode
                 !-------------------------
                 write(wtp, *)'  variables', '                gradient'
                 if (maxmode)then
                     do i = 1, N
                        write(wtp, *) X(i),-OBJGRA(i)
                     end do
                 else
                     do i = 1, N
                        write(wtp, *) X(i), OBJGRA(i)
                     end do
                 end if
      endsubroutine

      subroutine NonGshow(wtp,maxmode)
                 implicit none
                 integer,intent(in)::wtp
                 logical,intent(in)::maxmode
                 !-------------------------
                 write(wtp, *)'  variables'
                 do i = 1, N
                    write(wtp, *) X(i)
                 end do
      endsubroutine


      end subroutine ReportStat
!===============================================================================



!===============================================================================
      subroutine lsearc(SELF,N, X, NEWX, OBJFUN, NEWFUN, OBJGRA, DIREC, EPS, IND)
!===============================================================================

        implicit none

        class(minmaxoptimise),intent(inout)::self
        integer i

!*******************************************************************************
!ç›´ç·šæŽ¢ç´¢
! X : å¤‰æ•°ã®å€¤
! DIREC : é™ä¸‹æ–¹å‘
! DIRDER : å‹¾é…ãƒ™ã‚¯ãƒˆãƒ«ã®é™ä¸‹æ–¹å‘æˆåˆ†
! ST : ã‚¹ãƒƒãƒ†ãƒ—å¹…ï¼ˆtã®å€¤ï¼‰
! NEWX : æ›´æ–°å¾Œã®å¤‰æ•°ï¼ˆXNEW = X+ST*DIRECï¼‰
! NEWFUN : æ›´æ–°å¾Œã®é–¢æ•°å€¤
! EPS : ã‚¹ãƒ†ãƒƒãƒ—å¹…(ST)ã®ä¸‹é™
! BETA : tã®åŽæŸæ€§ã‚’èª²ã™æ¡ä»¶(10d-3)
! LAMBDA, MU : æ›´æ–°ã™ã‚‹ã‚¹ãƒ†ãƒƒãƒ—å¹…tã®é¸æŠžæ¡ä»¶
!******************************************************************************

        integer N, IND
        real(8) X(N), NEWX(N), OBJFUN, NEWFUN, OBJGRA(N), DIREC(N), &
                      DIRDER, DNORM, C, BETA, LAMBDA, MU, EPS, MINST, &
                      ST, ST1, ST2

        BETA = 1.0d-3
        LAMBDA = 1.0d-1
        MU = 5.0d-1

!-------------------------------------------------------------------------------

!-------compute the directional derivative.

        DIRDER = 0.0d0
        DNORM = 0.0d0

        do i = 1, N
           DIRDER = DIRDER + OBJGRA(i)*DIREC(i)
           DNORM = DNORM + DIREC(i)*DIREC(i)   !ãƒ™ã‚¯ãƒˆãƒ«ã®å†…ç©
        end do

        DNORM = sqrt(DNORM)

        MINST = min(1.0d0, EPS/DNORM) !åŽæŸåˆ¤å®šæ¡ä»¶

!-------set the initial step-size to be unity.

        ST = 1.0d0

!-------check if the current step-size is too small.

      2 if(ST < MINST) then
           IND = 5
           return
        end if

!-------test the line search criterion.

        do i = 1, N
           NEWX(i) = X(i) + ST*DIREC(i) !å¤‰æ•°ã®æ›´æ–°
        end do

        call fval(SELF, NEWX, NEWFUN) !æ›´æ–°ã•ã‚ŒãŸé–¢æ•°å€¤(NEWFUN)ã‚’å¾—ã‚‹

        if (self%ierr.eq.-9053) goto 999

        if (NEWFUN < OBJFUN + BETA*DIRDER*ST) return !å¤§åŸŸçš„åŽæŸæ¡ä»¶(1994)

!-------calculate the new step-size.

        C = (NEWFUN-OBJFUN-DIRDER*ST)/ST**2
        ST1 = -0.5*DIRDER/C
        ST2 = min(ST1, MU*ST)
        ST = max(LAMBDA*ST, ST2) !ã‚¹ãƒ†ãƒƒãƒ—å¹…ã®æ›´æ–°

        goto 2
    999 continue
        end subroutine lsearc
!===============================================================================


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++







      subroutine  CGmethod_my_frame(self , X ,OBJFUN, OBJGRA, IND)
                  implicit none
                  class(minmaxoptimise),intent(inout)::self
                  real(8),intent(inout)::x(self%N)
                  real(8),intent(inout)::OBJFUN
                  real(8),intent(inout)::OBJGRA(self%N)
                  integer::IND
                  !-------------------------------------------------
                  integer::mode,iexit
                  real(8)::hess(self%n**2),ww(self%n**2),dfn,xprmt(self%n),g(self%n),hh,deps
                  integer::jc
                  !==============
                   iexit=0
!                   iprint=0
                   hh=1.0e-4_8
                   mode=1
                   dfn=-0.5_8
                   deps=1.0e-5_8
                   hess=0.0_8
                   g=0.0_8     ;OBJGRA = 0._8
                   ww=0.0_8
                  !==============
                  xprmt = dabs(  x ) + 1.e-15_8
                  call minimize ( self, hh,deps, mode, self%MaxItration, self%showdetails, iexit,&
                                      OBJFUN, dfn, x, OBJGRA,hess, ww, xprmt)
      end subroutine
!

!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: minimize
! TYPE   : subroutine
! PURPOSE: conjugent gradient search
! I/O    :
! VERSION: 30-Sep-95
! COMMENT: This is a most reliable conjugent gradient routine! It has
!          served us well for many years, and is capable to cope with
!          a very large number of variables. Unfortunately, we don't
!          know who wrote this routine (original name: 'va10a'), and
!          we find it very obscure.
!          Don't worry, it works just fine.
!------------- usage --------------------------------
!       funct--- external function of n variables to be minimized.
!       n ------ the number of variables of funct. (input)
!       nm ----- the space of arrays set in main(). n<=nm
!       x ------ array x(nm), the first n element are the initializing
!        		   of the n variables. The rest (nm-n) element set zero.
!                (input)
!       f ------ the final minimum value of funct. (output)
!       w ------ array w(nm), the first n element are the final values
!                of n variables.
!       xm ----- array xm(nm), determined from x(nm) in main().
!       others: g,h,dfn,hh,eps,mode,maxfn,iprint,iexit
!         ------ usage not known. set as in the example.f.
!  Taken from the http://www.lps.ens.fr/~krauth , Lisadiag.f file.
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!   Action : put on all the declaration and use the "implicit none"
!   Author : Hengyue Li
!   DATE   : 20170811
!------------------------------------------------------------------------
!noprint================================================================
!      subroutine minimize (self, x, f, g, h, w, dfn, xm, hh, eps, mode, maxfn, iprint, iexit)
      subroutine minimize (self, hh, eps, mode, maxfn, iprint, iexit, f, dfn, x, g, h, w, xm)
        implicit none
    ! call minimize ( self, x, OBJFUN, OBJGRA,hess, ww, dfn, xprmt, hh,deps, mode, self%MaxItration, self%showdetails, iexit)
!      implicit double precision (a-h,o-z)
!      dimension x(*), g(*), h(*), w(*), xm(*)
      class(minmaxoptimise),intent(inout)::self
      integer::n,mode,maxfn, iprint, iexit
      real(8)::f,dfn,hh,eps
      real(8)::x,g,h,w,xm
      dimension::x(*),g(*),h(*),w(*),xm(*)


!      external funct   x(*), g(*), h(*), w(*), xm(*)

      !----------------------
      real(8)::zero,one,half,two,sig,tot,z,zz,gys,ff,dgs,f1,aeps,alpha,df,dmin,f2,gs0
      integer::np,nn,n1,link,i,j,iu,iv,jk,ik,itn,ij,int,is,k,ib,i1,idiff,ifn
      !----------------------

      n = self%N

      zero = 0.0d0
      half = 0.5d0
      one  = 1.0d0
      two  = 2.0d0


      if (iprint .ne. 0) write (6,1000)
 1000 format (' entry into minimize')
      np = n + 1
      n1 = n - 1
      nn=(n*np)/2
      is = n
      iu = n
      iv = n + n
      ib = iv + n
      idiff = 1
      iexit = 0
      if (mode .eq. 3) go to 15
      if (mode .eq. 2) go to 10
      ij = nn + 1
      do 5 i = 1, n
      do 6 j = 1, i
      ij = ij - 1
   6  h(ij) = zero
   5  h(ij) = one
      go to 15
  10  continue
      ij = 1
      do 11 i = 2, n
      z = h(ij)
      if (z .le. zero) return
      ij = ij + 1
      i1 = ij
      do 11 j = i, n
      zz = h(ij)
      h(ij) = h(ij) / z
      jk = ij
      ik = i1
      do 12 k = i, j
      jk = jk + np - k
      h(jk) = h(jk) - h(ik) * zz
      ik = ik + 1
  12  continue
      ij = ij + 1
  11  continue
      if (h(ij) .le. zero) return
  15  continue
      ij = np
      dmin = h(1)
      do 16 i = 2, n
      if (h(ij) .ge. dmin) go to 16
      dmin = h(ij)
  16  ij = ij + np - i
      if (dmin .le. zero) return
      z = f
      itn = 0
!      call funct (n, x, f)
      call fval(self,x,f)
      if (self%ierr.eq.-9053) goto 999
      ifn = 1
      df = dfn
      if (dfn .eq. zero) df = f - z
      if (dfn .lt. zero) df = abs (df * f)
      if (df .le. zero) df = one
  17  continue
      do 19 i = 1, n
      w(i) = x(i)
  19  continue
      link = 1
      if (idiff - 1) 100, 100, 110
  18  continue
      if (ifn .ge. maxfn) go to 90
  20  continue
      if (iprint .eq. 0) go to 21
      if (mod (itn, iprint) .ne. 0) go to 21
       write (6,1001) itn, ifn
1001  format (1x,'itn = ',i5,' ifn = ',i5)
      write (6,1002) f
1002  format (1x,'f = ',e15.7)
      if (iprint .lt. 0) go to 21
      write (6,1003) (x(i), i = 1, n)
!***
!***
1003  format (1x,'x = ',4e15.7 / (5x, 4e15.7))
      write (6,1004) (g(i), i = 1, n)
1004  format (1x,'g = ',4e15.7 / (5x, 4e15.7))
!**
!***
  21  continue
      itn = itn + 1
      w(1) = -g(1)
      do 22 i = 2, n
      ij = i
      i1 = i - 1
      z = -g(i)
      do 23 j = 1, i1
      z = z - h(ij) * w(j)
      ij = ij + n - j
  23  continue
  22  w(i) = z
      w(is+n) = w(n) / h(nn)
      ij = nn
      do 25 i = 1, n1
      ij = ij - 1
      z = zero
      do 26 j = 1, i
      z = z + h(ij) * w(is+np-j)
      ij = ij - 1
  26  continue
  25  w(is+n-i) = w(n-i) / h(ij) - z
      z = zero
      gs0 = zero
      do 29 i = 1, n
      if (z * xm(i) .ge. abs (w(is+i))) go to 28
      z = abs (w(is+i)) / xm(i)
  28  gs0 = gs0 + g(i) * w(is+i)
  29  continue
      aeps = eps / z
      iexit = 2
      if (gs0 .ge. zero) go to 92
      alpha = -two * df / gs0
      if (alpha .gt. one) alpha = one
      ff = f
      tot = zero
      int = 0
      iexit = 1
  30  continue
      if (ifn .ge. maxfn) go to 90
      do 31 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  31  continue
!      call funct (n, w, f1)
      call fval(self,w,f1)
      if (self%ierr.eq.-9053) goto 999
      ifn = ifn + 1
      if (f1 .ge. f) go to 40
      f2 = f
      tot = tot + alpha
  32  continue
      do 33 i = 1, n
      x(i) = w(i)
  33  continue
      f = f1
      if (int - 1) 35, 49, 50
  35  continue
      if (ifn .ge. maxfn) go to 90
      do 34 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  34  continue
!      call funct (n, w, f1)
      call fval(self,w,f1)
      if (self%ierr.eq.-9053) goto 999
      ifn = ifn + 1
      if (f1 .ge. f) go to 50
      if ((f1 + f2 .ge. f + f) .and.  &
       (7.0d0 * f1 + 5.0d0 * f2 .gt. 12.0d0 * f)) int = 2
      tot = tot + alpha
      alpha = two * alpha
      go to 32
  40  continue
      if (alpha .lt. aeps) go to 92
      if (ifn .ge. maxfn) go to 90
      alpha = half * alpha
      do 41 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  41  continue
!      call funct (n, w, f2)
      call fval (self,w,f2)
      if (self%ierr.eq.-9053) goto 999
      ifn = ifn + 1
      if (f2 .ge. f) go to 45
      tot = tot + alpha
      f = f2
      do 42 i = 1, n
      x(i) = w(i)
  42  continue
      go to 49
  45  continue
      z = 0.1d0
      if (f1 + f .gt. f2 + f2) &
        z = one + half * (f - f1) / (f + f1 - f2 - f2)
      if (z .lt. 0.1d0) z = 0.1d0
      alpha = z * alpha
      int = 1
      go to 30
  49  continue
      if (tot .lt. aeps) go to 92
  50  continue
      alpha = tot
      do 56 i = 1, n
      w(i) = x(i)        ! ;write(*,*)ib+i,ib,i,size(w)
      w(ib+i) = g(i)
  56  continue
      link = 2
      if (idiff - 1) 100, 100, 110
  54  continue
      if (ifn .ge. maxfn) go to 90
      gys = zero
      do 55 i = 1, n
      w(i) = w(ib+i)
      gys = gys + g(i) * w(is+i)
  55  continue
      df = ff - f
      dgs = gys - gs0
      if (dgs .le. zero) go to 20
      link = 1
      if (dgs + alpha * gs0 .gt. zero) go to 52
      do 51 i = 1, n
      w(iu + i) = g(i) - w(i)
  51  continue
      sig = one / (alpha * dgs)
      go to 70
  52  continue
      zz = alpha / (dgs - alpha * gs0)
      z = dgs * zz - one
      do 53 i = 1, n
      w(iu+i) = z * w(i) + g(i)
  53  continue
      sig = one / (zz * dgs * dgs)
      go to 70
  60  continue
      link = 2
      do 61 i = 1, n
      w(iu+i) = w(i)
  61  continue
      if (dgs + alpha * gs0 .gt. zero) go to 62
      sig = one / gs0
      go to 70
  62  continue
      sig = -zz
  70  continue
      w(iv+1) = w(iu+1)
      do 71 i = 2, n
      ij = i
      i1 = i - 1
      z = w(iu+i)
      do 72 j = 1, i1
      z = z - h(ij) * w(iv+j)
      ij = ij + n - j
  72  continue
      w(iv+i) = z
  71  continue
      ij = 1
      do 75 i = 1, n
      z = h(ij) + sig * w(iv+i) * w(iv+i)
      if (z .le. zero) z = dmin
      if (z .lt. dmin) dmin = z
      h(ij) = z
      w(ib+i) = w(iv+i) * sig / z
      sig = sig - w(ib+i) * w(ib+i) * z
      ij = ij + np - i
  75  continue
      ij = 1
      do 80 i = 1, n1
      ij = ij + 1
      i1 = i + 1
      do 80 j = i1, n
      w(iu+j) = w(iu+j) - h(ij) * w(iv+i)
      h(ij) = h(ij) + w(ib+i) * w(iu+j)
      ij = ij + 1
  80  continue
      go to (60, 20), link
  90  continue
      iexit = 3
      go to 94
  92  continue
      if (idiff .eq. 2) go to 94
      idiff = 2
      go to 17
  94  continue
      if (iprint .eq. 0) return
      write (6,1005) itn, ifn, iexit
1005  format (1x,'itn = ',i5, ' ifn = ',i5,' iexit = ',i5)
      write (6,1002) f
      write (6,1003) (x(i), i = 1, n)
      write (6,1004) (g(i), i = 1, n)
      return
 100  continue
      do 101 i = 1, n
      z = hh * xm(i)
      w(i) = w(i) + z
!      call funct (n, w, f1)
      call fval(self,w,f1)
      if (self%ierr.eq.-9053) goto 999
      g(i) = (f1 - f) / z
      w(i) = w(i) - z
 101  continue
      ifn = ifn + n
      go to (18, 54), link
 110  continue
      do 111 i = 1, n
      z = hh * xm(i)
      w(i) = w(i) + z
!      call funct (n, w, f1)
      call fval (self,w,f1)
      if (self%ierr.eq.-9053) goto 999
      w(i) = w(i) - z - z
!      call funct (n, w, f2)
       call fval(self,w,f2)
       if (self%ierr.eq.-9053) goto 999
      g(i) = (f1 - f2) / (two * z)
      w(i) = w(i) + z
 111  continue
      ifn = ifn + n + n
      go to (18, 54), link

  999 continue
      end subroutine minimize




 !@ OneDimensionalLocalMin
  subroutine OneDimensionalLocalMin(self, X, OBJFUN, OBJGRA, IND)
              implicit none
              class(minmaxoptimise),intent(inout)::self
              real(8)::X(SELF%N), OBJFUN, OBJGRA(SELF%N)
              integer,intent(inout)::IND
              !---------------------------------------------------
              integer::jc,interv
              real(8)::xt(self%n)
              interv = 0
              xt = 100._8
              do while(  sum((xt - x)**2) .gt. self%preX * dsqrt(1.0_8*self%n ) )
                interv = interv + 1
                !---------------------------------------------------------------
                if ( self%showdetails.ge.1)then
                   write(*,*)"+++++++++++++++++++"
                   write(*,*)"iteration:",interv
                   write(*,*)"delta=",sum(dabs(abs(xt -x)**2)),"-->",self%preX*dsqrt(1.0_8*self%n )
                   write(*,*)"+++++++++++++++++++"
                endif
                !---------------------------------------------------------------
                 xt = x
                 call OneDimLoc_one_direct_onestop_for_all_dim(self,x,OBJFUN)
                 if (self%ierr.eq.-9053) goto 999
              enddo

         999  continue
  endsubroutine


 !@ OneDimLoc_one_direct_onestop_for_all_dim
 subroutine OneDimLoc_one_direct_onestop_for_all_dim(self,x,y)
           implicit none
           class(minmaxoptimise),intent(inout)::self
           real(8)::X(SELF%N)
           real(8),intent(inout)::y
           !------------------------------
           integer::Dimi
           do Dimi = 1,int(self%n)
             !---------------------------------------------------------------
             if ( self%showdetails.ge.1)then
                write(*,*)"@@@@@@@@@@@@@@@@@@@@@@@@@@@"
                write(*,*)"Dimi=",Dimi,"/",int(self%n)
                write(*,*)"@@@@@@@@@@@@@@@@@@@@@@@@@@@"
             endif
             !---------------------------------------------------------------
              call OneDimLoc_one_direct_convergy(self,X,Dimi,y)
              if (self%ierr.eq.-9053) goto 999
           enddo
      999  continue
    endsubroutine


 !@ OneDimLoc_one_direct_convergy
 subroutine OneDimLoc_one_direct_convergy(self,X,Dimi,y)
   ! interativly find the min in one direction
   implicit none
   class(minmaxoptimise),intent(inout)::self
   real(8)::X(SELF%N)
   integer,intent(in)::Dimi
   real(8),intent(inout)::y
   !--------------------------------------------------
   integer::NR , pos , id
   real(8)::XR , xi  ,x1,x2,r1,r2,Prei,DeltaX,delta

   !---------------
   NR    = 20
   XR    = 0.1_8
   !---------------

   xi    = x(Dimi)
   x1    = xi - xr
   x2    = xi + xr
   Prei  = self%preX*0.1





   !--------------
   id = 1
!----- start -----------
100 continue
   call OneDimLoc_one_direct_one_step(self,X,Dimi,NR,x1,x2,pos,r1,r2,y,id)

                                                                       !;write(*,*)DIMI,pos,r1,r2,888
   if (self%ierr.eq.-9053) goto 999
   if (id.eq.-1)           goto 999

   id = 0
   IF (POS.eq.0)then
      ! ----- min is found roughly
      Do while( X2-X1.gt. Prei )
         x1 = r1
         x2 = r2
         call OneDimLoc_one_direct_one_step(self,X,Dimi,NR,x1,x2,pos,r1,r2,y,id) !;write(*,*)dimi,pos,r1,r2,X2-X1, Prei
         if (self%ierr.eq.-9053) goto 999
         if (id .eq. -1) goto 999
        !  if (pos.ne.0)then
        !     write(*,*)"why? min is suddenly missing"
        !     stop
        !  endif
      enddo
   else
      call random_number(delta)
      delta =  XR / (NR + 1)*delta  ! savty shift.
      x1    =  x1 + pos*2*xr - pos*delta
      x2    =  x2 + pos*2*xr - pos*delta
      goto 100
   endif
                                              !;write(*,*)DIMI,"is ok"
999 continue
endsubroutine



 !# OneDimLoc_one_direct_one_step
  subroutine OneDimLoc_one_direct_one_step(self,X,Dimi,Ntot,x1,x2,pos,rx1,rx2,y,id)
  ! for a multiple dimensional variant X, we check its Dimi conponent
  ! we search the min in range (x1,x2) with total number Ntot.
  ! find the minimer point and renew X
  ! One should check the output pos , if pos=-1, the min is allocated at x1, means the real min <x1
  !                                   if pos= 1, the min is allocated at x2, means the real min >x2
  !                                   if pos= 0, there is a real min in (x1,x2) and it have been setted to x
  ! if pos = 0, the found min is in (rx1,rx2)
            use basic_math_functions
            implicit none
            class(minmaxoptimise),intent(inout)::self
            real(8)::X(SELF%N)
            integer,intent(in)::Dimi,Ntot
            real(8),intent(in)::x1,x2
            integer,intent(out)::pos
            real(8),intent(out)::rx1,rx2,y
            integer,intent(inout)::id
            !-------------------------------------------------------------------
            real(8)::xlist(Ntot+1),ylist(Ntot+1),xtry(self%n),y1,sumzero
            integer::jc,p,np1
            TYPE(bmathf)::m

            !sumzero = 0.0000001_8
            sumzero = self%preY


            np1 = Ntot+1
            y1  = 0._8


            do jc = 1 , np1
               xlist(jc) = x1 + (x2-x1)/Ntot * (jc - 1)
               xtry = x
               xtry(Dimi) =  xlist(jc)
               call fval(SELF, xtry, ylist(jc))
               if (jc.eq.1) y1 = ylist(jc)
               ylist(jc) = ylist(jc) - y1
            enddo

          !------> id = -1 : y is almost 0.
          !------> id = 1  : the first time check the summation.



          !  if (id.eq.1)then
              if (sum(abs(ylist))/np1.le.sumzero) then
                 id  = -1
                 goto 999
              endif
          !  endif

            p = m%get_min_array_index(np1,ylist)

            if (p.eq.1) then
               pos = -1
               rx1 = xlist(1)
               rx2 = xlist(2)
            elseif(p.eq.np1)then
               pos = 1
               rx1 = xlist(np1-1)
               rx2 = xlist(np1)
            elseif(  (p.gt.1)   .and.   (p.lt.np1)    )then
               pos = 0
               rx1 = xlist(p-1)
               rx2 = xlist(p+1)
            else
               write(*,*)"unknow error in minimer@OneDimLoc_one_direct" ;stop
            endif
            x(Dimi) = xlist(p)       !;write(*,*)555,xlist,555
            y       = ylist(p) + y1  !;write(*,*)666,ylist+y1,666
            !--------------
        999 continue
            if ( self%showdetails.ge.2)then
               write(*,*)"-----------"
               write(*,*)"id=",id
               write(*,*)"pos=",pos
               write(*,*)"in direction:",int(Dimi),"/",int(self%n),"new x is:"
               write(*,*)x
               write(*,*)"New y is:"
               write(*,*)y
               write(*,*)"-----------"
            endif
endsubroutine
























!!##############################################################################
! from
!      Numerical Recipes in Fortran 77
!      The Art of Scientific Computing
!               Second Edition
!      Volume 1 of
!      Fortran Numerical Recipes
!      William H. Press
!      Saul A. Teukolsky
!      William T. Vetterling
!      Brian P. Flannery
!
!      10.4 Downhill Simplex Method in Multidimensions


subroutine  DownHill(self, X, OBJFUN, OBJGRA, IND)
            implicit none
            class(minmaxoptimise),intent(inout)::self
            real(8)::X(SELF%N), OBJFUN, OBJGRA(SELF%N)
            integer,intent(inout)::IND
            !---------------------------------------------------
            integer::mp,np,ndim ,iter ,jc
            real(8)::basilambda,p(SELF%N+1,SELF%N) ,y(SELF%N+1),ftol



            basilambda  =   0.1_8
            ftol        =   self%preY

            ndim = int(SELF%N)
            np   = int(SELF%N)
            mp   = int(SELF%N) + 1

            p(1,:) = x
            do jc = 1 , np
              p(jc+1,:)  = x
              p(jc+1,jc) = x(jc) + basilambda
            enddo
            do jc = 1 , mp
              call fval(SELF, p(jc,:), y(jc))
            enddo
            call amoeba(self,p,y,mp,np,ndim,ftol,iter)
            !---
            x = p(1,:)

          endsubroutine





SUBROUTINE amoeba(self,p,y,mp,np,ndim,ftol,iter)
implicit none
class(minmaxoptimise),intent(inout)::self
INTEGER iter,mp,ndim,np,NMAX,ITMAX
REAL*8 ftol,p(mp,np),y(mp),TINY
PARAMETER (NMAX=25,ITMAX=5000,TINY=1.e-10)
INTEGER i,ihi,ilo,inhi,j,m,n
REAL*8 rtol,sum,swap,ysave,ytry,psum(NMAX),gg(self%n)
iter=0
1 do  n=1,ndim
      sum=0._8
      do  m=1,ndim+1
          sum=sum+p(m,n)
      enddo
      psum(n)=sum
  enddo
2 ilo=1
  if (y(1).gt.y(2)) then
    ihi=1
    inhi=2
  else
    ihi=2
    inhi=1
  endif
do  i=1,ndim+1
if(y(i).le.y(ilo)) ilo=i
if(y(i).gt.y(ihi)) then
inhi=ihi
ihi=i
else if(y(i).gt.y(inhi)) then
if(i.ne.ihi) inhi=i
endif
enddo

if ( self%showdetails.ge.1)then
  !  write(self%wtp,*)
  !  write(self%wtp,*)"recent stage :"
  !  write(self%wtp,*)"X="
  !  write(self%wtp,*)p(ilo,:)
  !  write(self%wtp,*)"y="
  !  write(self%wtp,*)y(ilo)
  !  write(self%wtp,*)"--------------------------------"
   call ReportStat(self,3, p(ilo,:), y(ilo), gg, iter)
endif



! rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
! if (rtol.lt.ftol) then
  if (DownHill_Check_Converged(self,p,y,ilo,ihi))then
    swap=y(1)
    y(1)=y(ilo)
    y(ilo)=swap
    do  n=1,ndim
        swap=p(1,n)
        p(1,n)=p(ilo,n)
        p(ilo,n)=swap
    enddo
    return
endif




if (iter.ge.ITMAX) then
  write(*,*)"ITMAX exceeded in amoeba"
  stop
endif
iter=iter+2
ytry=amotry(self,p,y,psum,mp,np,ndim,ihi,-1.0_8)
if (ytry.le.y(ilo)) then
  ytry=amotry(self,p,y,psum,mp,np,ndim,ihi,2.0_8)
  else if (ytry.ge.y(inhi)) then
    ysave=y(ihi)
ytry=amotry(self,p,y,psum,mp,np,ndim,ihi,0.5_8)
if (ytry.ge.ysave) then
do  i=1,ndim+1
if(i.ne.ilo)then
  do  j=1,ndim
psum(j)=0.5_8*(p(i,j)+p(ilo,j))
p(i,j)=psum(j)
enddo
call fval(SELF, psum, y(i))
endif
enddo
iter=iter+ndim
goto 1
endif
else
iter=iter-1
endif
goto 2
ENDsubroutine



FUNCTION amotry(self,p,y,psum,mp,np,ndim,ihi,fac)
implicit none
class(minmaxoptimise),intent(inout)::self
INTEGER ihi,mp,ndim,np,NMAX
REAL*8 amotry,fac,p(mp,np),psum(np),y(mp),gg(1:self%n)
PARAMETER (NMAX=25)
INTEGER j
REAL*8 fac1,fac2,ytry,ptry(NMAX)
fac1=(1.-fac)/ndim
fac2=fac1-fac
do  j=1,ndim
ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
enddo


if ( self%showdetails.eq.2)then
  !  write(self%wtp,*)"[evaluate :]"
  !  write(self%wtp,*)"X="
  !  write(self%wtp,*)ptry(1:self%n)
  !  write(self%wtp,*)"y="
  !  write(self%wtp,*)ytry
  !  write(self%wtp,*)"--------------------------------"
   call ReportStat(self,3, ptry(1:self%n), ytry, gg, 0)
endif
!




call fval(SELF, ptry, ytry)
if (ytry.lt.y(ihi)) then
y(ihi)=ytry
do  j=1,ndim
psum(j)=psum(j)-p(ihi,j)+ptry(j)
p(ihi,j)=ptry(j)
enddo
endif
amotry=ytry
return
ENDfunction

logical function DownHill_Check_Converged(self,p,y,idmin,idmax)
  implicit none
  class(minmaxoptimise),intent(inout)::self
  real(8),intent(in)::p(SELF%N+1,SELF%N),y(SELF%N+1)
  integer,intent(in)::idmin,idmax
  !--------------------------------
  real*8,parameter::TINY = 1.e-13
  real*8::rtol,dx(self%n)
  logical::dyConverged,dxConverged
  integer::jc


  !-------------check y ------------------------------------
  rtol=2.*abs(y(idmax)-y(idmin))/(abs(y(idmax))+abs(y(idmin))+TINY)
  if (rtol.lt.self%preY)then
     dyConverged = .True.
  else
     dyConverged = .False.
  endif
  !------------check x -------------------------------------
  dx = p(idmax,:) - p(idmin,:)
  ! dxConverged = dsqrt(sum(dx**2)/self%n).le.self%preX


  dxConverged = .True.
  do jc = 1,self%n
    dxConverged = dxConverged .and. (abs(dx(jc)).le.self%preX)
  enddo
  !---------------------------------------------------------

  !---------------------------------------------------------
  !-----  how to chose the criterion
  ! DownHill_Check_Converged = dyConverged  .or.  dxConverged
  DownHill_Check_Converged = dxConverged
  !---------------------------------------------------------
endfunction













!!##############################################################################
! from
!      Numerical Recipes in Fortran 77
!      The Art of Scientific Computing
!               Second Edition
!      Volume 1 of
!      Fortran Numerical Recipes
!      William H. Press
!      Saul A. Teukolsky
!      William T. Vetterling
!      Brian P. Flannery
!
!      10.6 Conjugate Gradient Methods in Multidimensions





SUBROUTINE mnbrak(self,ax,bx,cx,fa,fb,fc,func)
  implicit none
  class(minmaxoptimise),intent(inout)::self
REAL*8 ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
EXTERNAL func
PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)
REAL*8 dum,fu,q,r,u,ulim
fa=func(self,ax)
fb=func(self,bx)
if(fb.gt.fa)then
dum=ax;ax=bx;bx=dum;dum=fb;fb=fa;fa=dum;endif;cx=bx+GOLD*(bx-ax)
fc=func(self,cx)
1 if(fb.ge.fc)then;r=(bx-ax)*(fb-fc);q=(bx-cx)*(fb-fa)
u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r));ulim=bx+GLIMIT*(cx-bx)
if((bx-u)*(u-cx).gt.0.)then
  fu=func(self,u)
  if(fu.lt.fc)then;ax=bx;fa=fb;bx=u;fb=fu;return
else if(fu.gt.fb)then;cx=u;fc=fu;return;endif;u=cx+GOLD*(cx-bx)
  fu=func(self,u)
else if((cx-u)*(u-ulim).gt.0.)then
  fu=func(self,u)
  if(fu.lt.fc)then;bx=cx;cx=u
u=cx+GOLD*(cx-bx);fb=fc;fc=fu
fu=func(self,u)
endif;else if((u-ulim)*(ulim-cx).ge.0.)then
u=ulim
fu=func(self,u)
else;u=cx+GOLD*(cx-bx)
  fu=func(self,u)
endif;ax=bx;bx=cx;cx=u;fa=fb
fb=fc;fc=fu;goto 1;endif;return
ENDSUBROUTINE



FUNCTION brent(self,ax,bx,cx,f,tol,xmin)
implicit none
class(minmaxoptimise),intent(inout)::self
INTEGER ITMAX
REAL*8 brent,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
EXTERNAL f
PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0e-10)
INTEGER iter
REAL*8::a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
a=min(ax,cx);b=max(ax,cx);v=bx;w=v;x=v;e=0._8;fx=f(self,x);fv=fx;fw=fx
do iter=1,ITMAX;xm=0.5*(a+b);tol1=tol*abs(x)+ZEPS;tol2=2.*tol1
  if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
  if(abs(e).gt.tol1) then;r=(x-w)*(fx-fv);q=(x-v)*(fx-fw);p=(x-v)*q-(x-w)*r
    q=2.*(q-r);if(q.gt.0.) p=-p
    q=abs(q);etemp=e;e=d
    if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) goto 1
    d=p/q;u=x+d;if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
    goto 2;endif
1 if(x.ge.xm) then
e=a-x;else;e=b-x;endif;d=CGOLD*e
2 if(abs(d).ge.tol1) then
 u=x+d;else;u=x+sign(tol1,d);endif;fu=f(self,u)
if(fu.le.fx) then;if(u.ge.x) then;a=x;else;b=x;endif;v=w;fv=fw;w=x;fw=fx;x=u;fx=fu;else
if(u.lt.x) then;a=u;else;b=u;endif;if(fu.le.fw .or. w.eq.x) then;v=w;fv=fw;w=u;fw=fu
else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then;v=u;fv=fu;endif
endif;enddo;write(*,*)"brent exceed maximum iterations";stop
3 xmin=x
brent=fx;return
ENDFUNCTION




REAL*8 FUNCTION f1dim(self,x)
implicit none
class(minmaxoptimise),intent(inout)::self
REAL*8::x!,func  !!!!!!!
!------------------------------
INTEGER NMAX
PARAMETER (NMAX=BNMAX)
INTEGER j!,ncom
REAL*8::xt(self%n)

do  j=1,self%n
xt(j)=pcom(j)+x*xicom(j)
enddo
! f1dim=func(xt)
call fval(self, xt, f1dim)
ENDfunction






!            << call gval(self, X, OBJGRA) >>
!            << call fval(self, X, OBJFUN) >>
!
!           class(minmaxoptimise),intent(inout)::self
!           real*8::x(self%n),objgra(self%n),objfun


! used : brent,f1dim,mnbrak
SUBROUTINE linmin(self,p,xi,n,fret)
  implicit none
  class(minmaxoptimise),intent(inout)::self
REAL*8 p(self%n),xi(self%n),TOL
INTEGER n
real*8::fret
PARAMETER (TOL=1.e-4)
INTEGER j
REAL*8 ax,bx,fa,fb,fx,xmin,xx
ncom=n;
do  j=1,self%n
    pcom(j)=p(j)
    xicom(j)=xi(j)
enddo
ax=0._8
xx=1._8
call mnbrak(self,ax,xx,bx,fa,fx,fb,f1dim)
fret=brent(self,ax,xx,bx,f1dim,TOL,xmin)
do  j=1,n
  xi(j)=xmin*xi(j)
  p(j)=p(j)+xi(j)
enddo                                      !; WRITE(*,*)xi,999090;STOP
ENDsubroutine



!USES dfunc,func,linmin
SUBROUTINE frprmn(self,p,fret,objG, IND) !(self, X, OBJFUN, OBJGRA, IND)
!SUBROUTINE frprmn(p,n,ftol,iter,fret)
! Given a starting point p that is a vector of length n , Fletcher-Reeves-Polak-Ribiere minimiza-
! tion is performed on a function func , using its gradient as calculated by a routine dfunc .
! The convergence tolerance on the function value is input as ftol . Returned quantities are
! p (the location of the minimum), iter (the number of iterations that were performed),
! and fret (the minimum value of the function). The routine linmin is called to perform
! line minimizations.
! Parameters: NMAX is the maximum anticipated value of n ; ITMAX is the maximum allowed
! number of iterations; EPS is a small number to rectify special case of converging to exactly
  ! zero function value.
implicit none
class(minmaxoptimise),intent(inout)::self
INTEGER iter,n,NMAX,ITMAX,IND
REAL*8 fret,ftol,p(self%n),EPS,objG(self%n)
PARAMETER (NMAX=BNMAX,ITMAX=200,EPS=1.e-10)
INTEGER its,j
REAL*8 dgg,fp,gam,gg,g(bNMAX),h(bNMAX),xi(bNMAX)

n = self%n
ftol=self%preY
objG = 0._8
IND = 0

!fp=func(p)
call fval(self, p, fp)
! call dfunc(p,xi)
call gval(self, p, xi(1:SELF%N))
do  j=1,n;g(j)=-xi(j);h(j)=g(j);xi(j)=h(j);enddo
do its=1,ITMAX;iter=its
  call linmin(self,p,xi,n,fret)
  !--------------
  if (self%showdetails.ge.1)then
     write(self%wtp,*)"---------------"
     write(self%wtp,*)"recent Pos:"
     write(self%wtp,*)xi(1:self%n)
  endif
  !--------------



  if(2.*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+EPS))return;fp=fret
  ! call dfunc(p,xi)
  call gval(self, p, xi(1:SELF%N))
  gg=0._8;dgg=0._8
  do  j=1,n;gg=gg+g(j)**2
  !-----------------------------------------------------------------------
    dgg=dgg+xi(j)**2           !This statement for Fletcher-Reeves.
    ! dgg=dgg+(xi(j)+g(j))*xi(j) !This statement for Polak-Ribiere.
  !-----------------------------------------------------------------------
  enddo
  if(gg.eq.0.)return
  gam=dgg/gg
  do  j=1,n;g(j)=-xi(j);h(j)=g(j)+gam*h(j)
  xi(j)=h(j)
  enddo
enddo
write(*,*)"frprmn maximum iterations exceeded";STOP
return
endsubroutine




!

SUBROUTINE dfpmin(self,p,fret,g,idn)!,func,dfunc)
  implicit none
  class(minmaxoptimise),intent(inout)::self
INTEGER iter,n,NMAX,ITMAX
REAL*8 fret,gtol,p(self%n),EPS,STPMX,TOLX!,func
PARAMETER (NMAX=BNMAX,ITMAX=200,STPMX=100.,EPS=3.e-8,TOLX=4.*EPS)
! EXTERNAL dfunc,func
! USES dfunc,func,lnsrch
! Given a starting point p(1:n) that is a vector of length n , the Broyden-Fletcher-Goldfarb-
! Shanno variant of Davidon-Fletcher-Powell minimization is performed on a function func ,
! using its gradient as calculated by a routine dfunc . The convergence requirement on zeroing
! the gradient is input as gtol . Returned quantities are p(1:n) (the location of the mini-
! mum), iter (the number of iterations that were performed), and fret (the minimum value
! of the function). The routine lnsrch is called to perform approximate line minimizations.
! Parameters: NMAX is the maximum anticipated value of n ; ITMAX is the maximum allowed
! number of iterations; STPMX is the scaled maximum step length allowed in line searches;
! TOLX is the convergence criterion on x values.

INTEGER i,its,j,idn
LOGICAL check
REAL*8 den,fac,fad,fae,fp,stpmax,sum,sumdg,sumxi,temp,test,&
dg(NMAX),g(self%n),hdg(NMAX),hessin(NMAX,NMAX),pnew(NMAX),xi(NMAX)


! (self, X, OBJFUN, OBJGRA, IND)
!------------------
n = self%n
gtol = self%preG
idn= 0
!------------------


! fp=func(p)
call fval(self,p(1:self%n),fp)
! call dfunc(p,g)
call gval(self,p(1:self%n),g(1:self%n))


sum=0._8;do i=1,n;do j=1,n;hessin(i,j)=0._8;enddo;hessin(i,i)=1._8
xi(i)=-g(i);sum=sum+p(i)**2;enddo;stpmax=STPMX*max(sqrt(sum),float(n));do its=1,ITMAX
iter=its;call lnsrch(self,n,p,fp,g,xi,pnew,fret,stpmax,check)

call ReportStat(self,self%showdetails, p(1:self%n), fp, g(1:self%n), its)

;fp=fret;do i=1,n
xi(i)=pnew(i)-p(i);p(i)=pnew(i);enddo;test=0._8;do i=1,n;temp=abs(xi(i))/max(abs(p(i)),1.)
if(temp.gt.test)test=temp;enddo;if(test.lt.TOLX)return;do i=1,n;dg(i)=g(i);enddo
  ! call dfunc(p,g)
  call gval(self,p(1:self%n),g(1:self%n))

test=0._8;den=max(fret,1.);do i=1,n;temp=abs(g(i))*max(abs(p(i)),1.)/den
if(temp.gt.test)test=temp;enddo;if(test.lt.gtol)return;do i=1,n;dg(i)=g(i)-dg(i);enddo
do i=1,n;hdg(i)=0._8;do j=1,n;hdg(i)=hdg(i)+hessin(i,j)*dg(j);enddo;enddo;fac=0._8
fae=0._8;sumdg=0._8;sumxi=0._8;do i=1,n;fac=fac+dg(i)*xi(i);fae=fae+dg(i)*hdg(i)
sumdg=sumdg+dg(i)**2;sumxi=sumxi+xi(i)**2;enddo;if(fac.gt.sqrt(EPS*sumdg*sumxi))then
fac=1._8/fac;fad=1._8/fae;do i=1,n;dg(i)=fac*xi(i)-fad*hdg(i);enddo;do i=1,n;do j=i,n
hessin(i,j)=hessin(i,j)+fac*xi(i)*xi(j)-fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j);hessin(j,i)=hessin(i,j)
enddo;enddo;endif;do i=1,n;xi(i)=0._8;do j=1,n;xi(i)=xi(i)-hessin(i,j)*g(j)
enddo;enddo;enddo;write(*,*)"too many iterations in dfpmin";stop;return
ENDsubroutine




SUBROUTINE lnsrch(self,n,xold,fold,g,p,x,f,stpmax,check)!,func)
  implicit none
class(minmaxoptimise),intent(inout)::self
INTEGER n;LOGICAL check;REAL*8::f,fold,stpmax,g(n),p(n),x(n),xold(n),ALF,TOLX!,func
PARAMETER (ALF=1.e-4,TOLX=1.e-7)
! EXTERNAL func
! USES func
! Given an n -dimensional point xold(1:n) , the value of the function and gradient there,
! fold and g(1:n) , and a direction p(1:n) , finds a new point x(1:n) along the direction
! p from xold where the function func has decreased “sufficiently.” The new function value
! is returned in f . stpmax is an input quantity that limits the length of the steps so that you
! do not try to evaluate the function in regions where it is undefined or subject to overflow.
! p is usually the Newton direction. The output quantity check is false on a normal exit.
! It is true when x is too close to xold . In a minimization algorithm, this usually signals
! convergence and can be ignored. However, in a zero-finding algorithm the calling program
! should check whether the convergence is spurious.
! Parameters: ALF ensures sufficient decrease in function value; TOLX is the convergence
! criterion on ∆x.
INTEGER i
REAL*8::a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,test,tmplam
!---------------
!TOLX= self%preX
!---------------
check=.false.;sum=0._8;do i=1,n;sum=sum+p(i)*p(i);enddo;sum=sqrt(sum);if(sum.gt.stpmax)then
do i=1,n;p(i)=p(i)*stpmax/sum;enddo;endif;slope=0._8;do  i=1,n;slope=slope+g(i)*p(i)
enddo;if(slope.ge.0.)then;write(*,*)"roundoff problem in lnsrch";stop
 endif;test=0._8;do i=1,n;temp=abs(p(i))/max(abs(xold(i)),1.)
if(temp.gt.test)test=temp;enddo;alamin=TOLX/test;alam=1._8
1 continue
do i=1,n;x(i)=xold(i)+alam*p(i);enddo;call fval(self,X(1:self%n),f)
! f=func(x)

! !-------------------------
! !  print
! if (self%showdetails.ge.1)then
!   write(self%wtp,*)"position:"
!   write(self%wtp,*)x(1:self%n)
!   write(self%wtp,*)"y=",f
!   write(self%wtp,*)"--------------"
! endif
! !------------------------


if(alam.lt.alamin)then;do i=1,n;x(i)=xold(i);enddo
check=.true.;return;else if(f.le.fold+ALF*alam*slope)then;return;else;if(alam.eq.1.)then
tmplam=-slope/(2.*(f-fold-slope));else;rhs1=f-fold-alam*slope;rhs2=f2-fold-alam2*slope
a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2);b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
if(a.eq.0.)then;tmplam=-slope/(2.*b);else;disc=b*b-3.*a*slope;if(disc.lt.0.)then;tmplam=.5*alam
else if(b.le.0.)then;tmplam=(-b+sqrt(disc))/(3.*a);else;tmplam=-slope/(b+sqrt(disc));endif;endif
if(tmplam.gt..5*alam)tmplam=.5*alam;endif;endif;alam2=alam;f2=f;alam=max(tmplam,.1*alam)
goto 1
endsubroutine


endmodule
