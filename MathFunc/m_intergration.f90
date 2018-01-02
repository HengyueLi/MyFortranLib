



!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : class_integration
! OBJECT: TYPE(Integration)
! USED  : CodeObject , LegendreFun
! DATE  : 2018-01-02
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!        Creat list used to calculate the needed integration.
!        kind=8, input can be real(8) and complex(8)
!        For a type(Integration)::I,    first call I%Initialization(N,Itype)   ; then call I%get_xw(x1,x2,x,w) to get x,w
!        (x,w) can be real(8) or complex(8)
!        Itype=1 for Gauss_Legendre_Method,
!        Itype=0 for a returning homogeneous contribution(w(i)=1/abs(x2-x1) )
! STANDARD:
!        *CALL  I%Initialization(N,Itype)
!               Itype=1 for Gauss_Legendre_Method,
!               Itype=0 for a returning homogeneous contribution(w(i)=1/abs(x2-x1) )
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!                  [sub] Initialization(N,Itype)
!                        integer::N,Itype
!
!                  [sub] get_xw(xs,xe,x,w)
!
!
! avalable gets:
!                  [sub]
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


module class_integration
  use CodeObject
  use LegendreFun ,only: LegF => func
IMPLICIT NONE



 type,extends(object)::Integration
      integer,private::typecase         !1 for GS, 0 for homoginious contribution
      INTEGER,private::TOTN

     contains
     procedure,pass::Initialization

     generic ::get_xw=>get_xw_real,get_xw_complex
     procedure,PRIVATE,pass::get_xw_real
     procedure,PRIVATE,pass::get_xw_complex
 end type



private::get_xw_real,get_xw_complex,Initialization
private::GetOneXW

 contains




    subroutine Initialization(self,N,typecase)
            implicit none
            class(Integration),intent(inout)::self
            integer,intent(in)::n,typecase
            !---------------------------------
            call self%SetInitiated(.true.)
            self%TOTN=N
            self%typecase=typecase
    end subroutine


    SUBROUTINE get_xw_real(self,xs,xe,x,w)
        implicit none
        class(Integration),intent(inout)::self
        real(8),intent(in)::xs,xe
        real(8),intent(out)::x(self%TOTN),w(self%TOTN)
        !--------------------------------------------------
        integer::jc
        call self%CheckInitiatedOrStop()
        do jc = 1 , self%TOTN
          call GetOneXW(self,self%TOTN,jc,X(jc),W(jc))
        enddo
        do jc=1,self%TOTN
            x(jc)=x(jc)*(xe-xs)/2._8+(xs+xe)/2._8
            w(jc)=w(jc)*(xe-xs)/2._8
        end do
    END SUBROUTINE

    SUBROUTINE get_xw_complex(self,xs,xe,x,w)
        implicit none
        class(Integration),intent(inout)::self
        complex(8),intent(in)::xs,xe
        complex(8),intent(out)::x(self%TOTN)
        real(8)   ,intent(out)::w(self%TOTN)
        !--------------------------------------------------
        real(8)::xt(self%TOTN),wt(self%TOTN)
        real(8)::le
        integer::jc
        call self%CheckInitiatedOrStop()
        le=zabs(xs-xe) !;write(*,*)le
        do jc = 1 , self%TOTN
          call GetOneXW(self,self%TOTN,jc,Xt(jc),Wt(jc))
        enddo
        do jc=1,self%TOTN
            x(jc)=xt(jc)*(xe-xs)/2._8+(xe+xs)/2._8
            w(jc)=wt(jc)*le/2._8
        end do
    END SUBROUTINE



    subroutine GetOneXW(self,n,i,x,w)
      implicit none
      class(Integration),intent(inout)::self
      integer,intent(in)::n,i
      real*8,intent(out)::x,w
      !--------------------------------------------------
      TYPE(LegF)::f
      real*8::delta

      select case(self%typecase)
      case(0)  ! homoginious distribution
        delta=2._8/N
        X=(n-i)*delta + 0.5_8 * delta -1._8
        W= 2._8/N
      case(1)  ! GS distribution
        X = f%GetZerosOfLegendre(n,i)
        W=2.0_8*(1.0_8-X**2.0_8)/( f%GetLegendreP(N+1,X) *(N+1))**2.0_8
      case default
        write(self%getprint(),*)"ERROR:input type is wrong in module Integration";stop
      endselect
    endsubroutine

end module












! module class_integration
! IMPLICIT NONE
!
!
!
!  type Integration
!       logical,private::Initialized=.false.
!       integer(8),private::typecase         !1 for GS, 0 for homoginious contribution
!       INTEGER(8),private::TOTN
!
!      contains
!      procedure,pass::Initialization
!
!      generic ::get_xw=>get_xw_real,get_xw_complex
!      procedure,PRIVATE,pass::get_xw_real
!      procedure,PRIVATE,pass::get_xw_complex
!  end type
!
!
!
!
! private::PARAMETER_FOR_INTERGRAN2015_GS,PARAMETER_FOR_INTERGRAN2015_ho,LEGENDRE2015,ZERO_FOR_LEGENDRE2015,DE_FOR_LEGENDRE2015&
!          ,get_xw_real,get_xw_complex,Initialization,check_Initialized,typewrong
!
!  contains
!
!     subroutine check_Initialized(self)
!         implicit none
!         type(Integration),intent(in)::self
!         if (.not.self%initialized)then
!             write(*,*)"Integration type is not initialized yet";stop
!         end if
!     end subroutine
!
!
!
!     subroutine Initialization(self,N,typecase)
!             implicit none
!             class(Integration),intent(inout)::self
!             integer(8),intent(in)::n,typecase
!             !---------------------------------
!             self%Initialized=.true.
!             self%TOTN=N
!             self%typecase=typecase
!     end subroutine
!
!
!     SUBROUTINE get_xw_real(self,xs,xe,x,w)
!         implicit none
!         class(Integration),intent(inout)::self
!         real(8),intent(in)::xs,xe
!         real(8),intent(out)::x(self%TOTN),w(self%TOTN)
!         !--------------------------------------------------
!         real(8)::pre
!         integer(8)::jc
!         call check_Initialized(self)
!         PRE=0.000000000000001_8!PRE=10._8**(-(2*8-1))
!         select case(self%typecase)
!             case(0_8)
!                 do jc=1_8,self%TOTN
!                    call PARAMETER_FOR_INTERGRAN2015_HO(int(self%TOTN,4),int(self%TOTN,4)-int(jc,4)+1,X(jc),W(jc),pre)
!                 end do
!             case(1_8)
!                 do jc=1_8,self%TOTN
!                    call PARAMETER_FOR_INTERGRAN2015_gs(int(self%TOTN,4),int(self%TOTN,4)-int(jc,4)+1,X(jc),W(jc),pre)
!                 end do
!             case default
!                 call typewrong()
!         end select
!         do jc=1_8,self%TOTN
!             x(jc)=x(jc)*(xe-xs)/2._8+(xs+xe)/2._8
!             w(jc)=w(jc)*(xe-xs)/2._8
!         end do
!     END SUBROUTINE
!
!     SUBROUTINE get_xw_complex(self,xs,xe,x,w)
!         implicit none
!         class(Integration),intent(inout)::self
!         complex(8),intent(in)::xs,xe
!         complex(8),intent(out)::x(self%TOTN)
!         real(8)   ,intent(out)::w(self%TOTN)
!         !--------------------------------------------------
!         real(8)::xt(self%TOTN),wt(self%TOTN)
!         real(8)::pre,le
!         integer(8)::jc
!         call check_Initialized(self)
!         PRE=0.000000000000001_8    !10._8**(-(2*8-1))
!         le=zabs(xs-xe) !;write(*,*)le
!         select case(self%typecase)
!             case(0_8)
!                 do jc=1_8,self%TOTN
!                    call PARAMETER_FOR_INTERGRAN2015_HO(int(self%TOTN,4),int(self%TOTN,4)-int(jc,4)+1,Xt(jc),Wt(jc),pre)
!                 end do
!             case(1_8)
!                 do jc=1_8,self%TOTN
!                    call PARAMETER_FOR_INTERGRAN2015_GS(int(self%TOTN,4),int(self%TOTN,4)-int(jc,4)+1,Xt(jc),Wt(jc),pre)
!                 end do
!             case default
!                 call typewrong()
!         end select
!         do jc=1_8,self%TOTN
!             x(jc)=xt(jc)*(xe-xs)/2._8+(xe+xs)/2._8
!             !w(jc)=wt(jc)*(xe-xs)/2._8
!             w(jc)=wt(jc)*le/2._8
!         end do
!     END SUBROUTINE
!
!     subroutine typewrong()
!         implicit none
!         write(*,*)"input type wrong in module Integration";stop
!     end subroutine
!
!
!
! SUBROUTINE PARAMETER_FOR_INTERGRAN2015_HO(N,I,XI,WI,RE)  !
! IMPLICIT NONE
! INTEGER::N,I   ;REAL(8)::XI,WI,RE,delta
! delta=2._8/N
! XI=(n-i)*delta + 0.5_8*delta -1._8
! WI= 2._8/N
!
! ENDSUBROUTINE
!
!
! SUBROUTINE PARAMETER_FOR_INTERGRAN2015_GS(N,I,XI,WI,RE)  !
! !The details should be found in Yuan-Yao He's instruction
! !  N   IN-PUT    INTEGER     the order of LEGENDRE2015 function
! !  I   IN-PUT    INTEGER     the needed i-th zero of the correspongding LEGENDRE2015 function
! !  XI  OUT-PUT   8    the posion of the i-th zero
! !  WI  OUT-PUT   8    the corresponding weight number
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IMPLICIT NONE
! INTEGER::N,I   ;REAL(8)::XI,WI,RE
! XI=ZERO_FOR_LEGENDRE2015(N,I,RE)
! WI=2.0_8*(1.0_8-XI**2.0_8)/(LEGENDRE2015(N+1,XI)*(N+1))**2.0_8
! ENDSUBROUTINE
!
!
!
!
!
!
!
!  FUNCTION LEGENDRE2015(N,X)    !P_N(X)
!     IMPLICIT NONE
!     REAL(8)::LEGENDRE2015,X  ;    INTEGER::N
!     INTEGER::JC  ;REAL(8)::TEMPK_2,TEMPK_1
!     LEGENDRE2015=1.0_8
!     IF(N.EQ.0)THEN
!
!     ELSEIF(N.EQ.1)THEN
!         LEGENDRE2015=X
!     ELSE
!         TEMPK_2=LEGENDRE2015
!         LEGENDRE2015=X
!        DO JC=2,N
!            TEMPK_1=TEMPK_2
!            TEMPK_2=LEGENDRE2015
!            LEGENDRE2015=1.0_8/JC*( (2.0_8*JC-1)*X*TEMPK_2-TEMPK_1*(JC-1) )
!            TEMPK_1=TEMPK_2
!        ENDDO
!
!     ENDIF
!     ENDFUNCTION
!
! FUNCTION ZERO_FOR_LEGENDRE2015(N,I,RE) !THE I-TH ZERO OF THE N ODER LEGENDRE2015,RE IS THE precision
! IMPLICIT NONE
! REAL(8)::ZERO_FOR_LEGENDRE2015,RE    ;INTEGER::N,I
! REAL(8)::TEMP1
!
!     TEMP1=863.32 !RAMDON
!
! ZERO_FOR_LEGENDRE2015=COS(  (4.0_8*I-1.0_8)/(4.0_8*N+2)*3.1415926535897932384_8  )
! DO WHILE(ABS(ZERO_FOR_LEGENDRE2015-TEMP1).GT.RE)
!     TEMP1=ZERO_FOR_LEGENDRE2015
!     ZERO_FOR_LEGENDRE2015=ZERO_FOR_LEGENDRE2015-LEGENDRE2015(N,ZERO_FOR_LEGENDRE2015) &
!     /DE_FOR_LEGENDRE2015(N,ZERO_FOR_LEGENDRE2015,RE)
! ENDDO
!     ENDFUNCTION
!
! FUNCTION DE_FOR_LEGENDRE2015(N,X,RE)
!     IMPLICIT NONE
!     REAL(8) DE_FOR_LEGENDRE2015,X,RE   ;INTEGER::N
!     REAL(8)::TEMP
!                           TEMP=0.1_8*RE
!     !DE_FOR_LEGENDRE2015=(LEGENDRE2015(N,X+TEMP)-LEGENDRE2015(N,X))/TEMP
!     DE_FOR_LEGENDRE2015=n*(x*LEGENDRE2015(N,X) - LEGENDRE2015(N-1,X) )/(x**2.0_8-1._8)
!     ENDfunction
!
! end module
