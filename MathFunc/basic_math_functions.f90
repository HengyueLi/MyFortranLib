



!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : basic_math_functions
! OBJECT: TYPE(bmathf)
! USED  :
! DATE  : 2017-12-05
! AUTHOR: Hengyue Li
!--------------
! DESCRIPTION:
!             some basic math functions
! STANDARD:
!
!
!
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!
! avalable gets:
!                  [fun] get_random_int(start,endto)
!                         Get a random integer in range [start,endto]
!                         if input are integer*4, return integer*4
!                         if input are integer*8, return integer*8
!
!                  [sub] get_x_range_equal_interval_p(nx,start,endto,x):
!                             integer(8)nx; real(8)::start,endto,x(nx)
!                             creat a list of real(8) value range from start to endto
!
!                  [sub] get_random_array(A)
!                        real(8)/complex(8)::A(:)
!
!                  [sub] get_normalized_array(A)
!                        The same as get_radom_array, but keepy <A|A> = 1
!
!                  [fun] bcombination(m,n)
!                        integer(8)::m,n  /  integer::m,n
!                        integer*8::bcombination
!                        bcombination=c^m_n
!
!                  [fun] get_min_array_index(n,x)
!                        integer::get_min_array_index,n
!                        real(8)::x(n)
!                        return the position of the min point in a array.
!
! avalable Is :
!                  [fun] IsTwoValuePercentageTheSame(E1,E2,rePre)
!                        check two real*8 E1,E2 are the same or not with relavant different rePre
!
!                  [fun] IsIntegerArrayTheSame(IntArry)
!                        input an integer array  check if all elemets are the same
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

module basic_math_functions
    implicit none



    TYPE bmathf
         private


         contains
         procedure,pass::get_random_int8_p
         procedure,pass::get_random_int4_p
         generic::get_random_int => get_random_int8_p, get_random_int4_p

         procedure,pass::get_x_range_equal_interval=>get_x_range_equal_interval_p

         procedure,pass::get_min_array_index

         procedure,nopass::bcombination8
         procedure,nopass::bcombination4
         generic::bcombination => bcombination8 , bcombination4

         procedure,private,PASS::get_radom_array_r
         procedure,private,PASS::get_radom_array_C
         procedure,private,PASS::get_normalized_array_rp
         procedure,private,PASS::get_normalized_array_Cp
         generic::get_random_array=>get_radom_array_r,get_radom_array_C
         generic::get_normalized_array=>get_normalized_array_Rp,get_normalized_array_Cp

         procedure,nopass::IsTwoValuePercentageTheSame
         procedure,nopass::IsIntegerArrayTheSame

    END TYPE

   private::get_random_int8_p,get_random_int4_p,get_x_range_equal_interval_p
   PRIVATE::get_radom_array_r,get_radom_array_C,get_normalized_array_rp,get_normalized_array_Cp
   private::bcombination8,bcombination4
   private::get_min_array_index
   private::IsTwoValuePercentageTheSame
   private::IsIntegerArrayTheSame
  !-----------------------
    contains



    integer(8) function get_random_int8_p(self,start,endto)
               implicit none
               class(bmathf),intent(inout)::self
               integer(8),intent(in)::start,endto
               !------------------------------
               real(8)::r
               if ( start .le. endto    )then
                   call random_number(r)
                   get_random_int8_p = int( r * (  endto  -  start +1_8 ) ,8) + start
               else
                    write(*,*)"in get_random_int8, start <= endto should be kept,error201706051031";stop
               end if
    end function

    integer function get_random_int4_p(self,start,endto)
      implicit none
      class(bmathf),intent(inout)::self
      integer,intent(in)::start,endto
      !------------------------------
      real(8)::r
      if ( start .le. endto    )then
          call random_number(r)
          get_random_int4_p = int( r * (  endto  -  start +1 ) ,4) + start
      else
           write(*,*)"in get_random_int4, start <= endto should be kept,error201706051031";stop
      end if
end function


    subroutine get_x_range_equal_interval_p(self,nx,start,endto,x)
               implicit none
               class(bmathf),intent(inout)::self
               integer(8),intent(in)::nx
               real(8),intent(in)::start,endto
               real(8),intent(out)::x(nx)
               !----------------
               real(8)::deltax
               integer(8)::jc
               deltax = (endto - start)/nx
               do jc=1_8,nx
                  x(jc) = start + ( jc - 1_8 ) * deltax + deltax/2._8
               end do
    end subroutine

    subroutine get_radom_array_r(self,A)
               implicit none
               class(bmathf),intent(inout)::self
               real(8),intent(out)::A(:)
               !-------------------------------------
               call random_number(A)
    end subroutine

    subroutine get_radom_array_c(self,A)
               implicit none
               class(bmathf),intent(inout)::self
               COMPLEX(8),intent(out)::A(:)
               !-------------------------------------
               integer(8)::n
               real(8),allocatable::R(:),C(:)
               n = size(A)
               allocate( R(n) )
               allocate( C(n) )
               call get_radom_array_r(self,R)
               call get_radom_array_r(self,C)
               A = CMPLX( R  ,  C   ,8   )
               deallocate(R)
               deallocate(C)
    end subroutine

    subroutine get_normalized_array_rp(SELF,A)
               implicit none
               class(bmathf),intent(inout)::self
               real(8),intent(out)::A(:)
               !-------------------------------------
               call get_radom_array_r(self,A)
               A = A/dsqrt( sum(A**2) )
    end subroutine

    subroutine get_normalized_array_cp(SELF,A)
               implicit none
               class(bmathf),intent(inout)::self
               COMPLEX(8),intent(out)::A(:)
               !-------------------------------------
               call get_radom_array_c(self,A)
               A = A/dsqrt(  REAL(  sum(A*  CONJG(A) ) )   )
    end subroutine


    integer(8) function bcombination8(mu,n) !c^mu_n
           implicit none
           integer(8),intent(in)::mu,n
           !---------------------
           integer(8)::jc,m
           m = min( mu, n-mu )
           bcombination8=1_8
           do jc=n,n-m+1_8,-1_8
              bcombination8=bcombination8*jc   !;write(*,*)jc,bcombination
           end do
           do jc=1_8,m
              bcombination8=bcombination8/jc
           end do
    end function


    integer*8 function bcombination4(mu,n)
      implicit none
      integer,intent(in)::mu,n
      !---------------------
      integer*8::mu8,n8
      mu8 = mu * 1_8
      n8  = n  * 1_8
      bcombination4 = bcombination8(mu8,n8)
    endfunction



    integer function get_min_array_index(self,n,x)
            implicit none
            class(bmathf),intent(inout)::self
            integer,intent(in)::n
            real(8),intent(in)::x(n)
            !----------------------------------------
            integer::jc
            get_min_array_index = 1
            do jc = 2 , n
                if (  x(get_min_array_index) .ge.   x(jc)    )then
                    get_min_array_index = jc
                end if
            end do
    end function






    function IsTwoValuePercentageTheSame(e1,e2,Repre) result(CheckDiffSmall)
      implicit none
      real*8,intent(in)::E1,E2,rePre
      logical::CheckDiffSmall
      !----------------------------------------
      if (E1.eq.E2)then
        CheckDiffSmall = .true.
      else
        if ( abs(E1-E2) / max( abs(e1) , abs(e2)   ) .le.rePre)then
           CheckDiffSmall = .true.
        else
          CheckDiffSmall = .false.
        endif
      endif
    endfunction


    logical function IsIntegerArrayTheSame(IntArry)
      implicit none
      integer,intent(in)::IntArry(:)
      !----------------------------------------------
      integer::s ,jc
      integer,allocatable::t(:)
      s = size(IntArry)
      allocate(t(s))
      s = IntArry(1)
      if (  sum(abs(s - IntArry)).eq.0   )then
        IsIntegerArrayTheSame = .true.
      else
        IsIntegerArrayTheSame = .false.
      endif
    endfunction






end module
