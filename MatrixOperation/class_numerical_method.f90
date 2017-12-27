



!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : class_numerical_method
! OBJECT: TYPE(nummethod)
! USED  : @lapck
! DATE  : 2017-12-06
! AUTHOR: Hengyue Li
!--------------
! DESCRIPTION:
!            only method here in this module.  Therefore no Initialization will be needed.
! STANDARD:
!
!
!          [sub]  ED_Hermitian_matrix(UL,N,H,E):
!                 character::UL = “U” or “L”
!                 integer::N
!                 real(8)::E(n)
!                 INPUT/OUTPUT H, H can be real(8) or complex(8)
!
!          [sub]  ED_tridiagonal_real(N,A,B,H,E):
!                 input integer::N,  REAL(8)::A(N),B(N)
!                 OUTPUT real(8)::H(n,n) ,E(n)
!                 B(1) will not be used.
!
!          [sub] MatrixInverse(N,M)
!                integer::N
!                real*8/complex*16,intent(inout)::M(n,n)
!                calculate the inverse of the matrix M
!
!          [fun] get_determinant_of_Matrix(N,M)
!                integer::N
!                real*8/complex*16::M(N)
!                return real*8/complex*16, the determinant of matrix M(N)
!
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


module class_numerical_method
   implicit none

   type nummethod



    contains
       procedure,private,pass::ED_HERMIT_real8
       procedure,private,pass::ED_HERMIT_complex8
       generic::ED_Hermitian_matrix=>ED_HERMIT_real8,ED_HERMIT_complex8

       procedure,pass::ED_tridiagonal_real=>ED_tridiagonal_real_p

       procedure,pass::MatrixInverse_r
       procedure,pass::MatrixInverse_c
       generic::MatrixInverse=>MatrixInverse_r,MatrixInverse_c

       procedure,pass::get_determinant_of_Matrix_r
       procedure,pass::get_determinant_of_Matrix_c
       generic::get_determinant_of_Matrix=>get_determinant_of_Matrix_r,get_determinant_of_Matrix_c




   end type



private::ED_HERMIT_real8,ED_HERMIT_complex8
private::MatrixInverse_r,MatrixInverse_c
private::get_determinant_of_Matrix_c,get_determinant_of_Matrix_r

contains

    subroutine ED_HERMIT_real8(self,UL,n,H,E)
    !---------------------------------------------------------!
    !Calls the LAPACK diagonalization subroutine DSYEV        !
    !input:  h(n,n) = real symmetric matrix to be diagonalized!
    !            n  = size of a                               !
    !            UL="L"/"U" choosing lower or upper part of H.!
    !output: h(n,n) = orthonormal eigenvectors of a           !
    !        e(n) = eigenvalues of a in ascending order       !
    !---------------------------------------------------------!
        implicit NONE
        class(nummethod),intent(inout)::self
        character(1),INTENT(IN)::UL
        integer,intent(in)::n
        real(8),intent(inout)::h(n,n),e(n)
        !--------------------------------------
        INTEGER::N4,L,inf
        REAL(8)::WORK( n*(3+n/2)  )
        N4=INT(N,4)
        l=n4*(3+n4/2)
        call dsyev('V',UL,n4,H,n4,E,work,l,inf)
    endsubroutine


    subroutine ED_HERMIT_complex8(self,UL,n,H,E)
        implicit NONE
        class(nummethod),intent(inout)::self
        character(1),INTENT(IN)::UL
        integer,intent(in)::n
        complex(8),intent(inout)::h(n,n)
        real(8),intent(inout)::e(n)

        call LA_ziasym(UL,h,e,int(n,4))
    contains
             subroutine LA_ziasym(UL,a,eig,n)
             implicit none
             character(1),INTENT(IN)::UL
             integer n,l,inf
             COMPLEX(8)::a(n,n),work(n*(3+n/2))  ;REAL(8)::eig(n)   ;REAL(8),ALLOCATABLE::RWORK(:)

             l=n*(3+n/2)
             ALLOCATE(RWORK(max(1, 3*N-2)))
             call ZHEEV('V',UL,n,a,n,eig,work,l,RWORK,inf)
             DEALLOCATE(RWORK)
             end subroutine LA_ziasym
            !---------------------!
    endsubroutine



   subroutine ED_tridiagonal_real_p(self,N,A,B,H,E,ifo) ! B(1) will not be used.
        implicit NONE
        class(nummethod),intent(inout)::self
        integer,intent(in)::n
        real(8),intent(in)::a(n),b(n)
        real(8),intent(out)::h(n,n),e(n)
        integer,intent(out),optional::ifo
        !-------------------------------------
        real*8::bn_1(n-1)
        integer::LWORK,n4,LIWORK,INFO
        REAL*8,ALLOCATABLE::WORK(:)
        INTEGER,ALLOCATABLE::IWORK(:)
        n4 = int(n,4)
        E    = A
        bn_1 = b(2:n)
        LWORK  = 1 + 4*N4 + N4**2
        LIWORK = 3+5*N4
        INFO = 0
        ALLOCATE( IWORK(LIWORK)   ,WORK(LWORK))
        call DSTEVD( "V", N4, E, bn_1, H, max(1,N4), WORK, LWORK, IWORK, LIWORK, INFO )
       DEALLOCATE( IWORK ,WORK)
       if (present(ifo)) ifo = info


   end subroutine






!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   subroutine MatrixInverse_r(self,N,M)
     implicit NONE
     class(nummethod),intent(inout)::self
     INTEGER,intent(in)::N
     real*8,intent(inout)::M(N,N)
     !
     INTEGER::IPIV(N),INFO, LWORK
     real*8::WORK(N*32)
     integer::n4

     LWORK=-1
     CALL DGETRF(n , n , M, n , IPIV, INFO )
     CALL DGETRI( n, M, n, IPIV, WORK, N*32, INFO )
 ENDSUBROUTINE
   subroutine MatrixInverse_c(self,N,M)
     implicit NONE
     class(nummethod),intent(inout)::self
     INTEGER,intent(in)::N
     COMPLEX*16,intent(inout)::M(N,N)
     !
     INTEGER::IPIV(N),INFO, LWORK
     COMPLEX*16::WORK(N*32)
     integer::n4

     LWORK=-1
     CALL ZGETRF(n , n , M, n , IPIV, INFO )
     CALL ZGETRI( n, M, n, IPIV, WORK, N*32, INFO )
 ENDSUBROUTINE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 function get_determinant_of_Matrix_c(self,N,M) result(R)
 IMPLICIT none
 class(nummethod),intent(inout)::self
 INTEGER,intent(in)::N
 COMPLEX*16,intent(in)::M(N,N)
 complex*16::R
 !--------------------------------------------------------
 COMPLEX*16::A(N,N)  ;INTEGER::IPIV(N),INFO ;INTEGER::JC
 A=M
 call ZGETRF( n,	n, A, n, IPIV, INFO )
 R=(1._8,0._8)
 DO JC=1,INT(N,4)
   R=R*A(JC,JC)
   IF (.NOT. JC.EQ.IPIV(JC)) R=(-1._8)*R
 ENDDO
 END function

 function get_determinant_of_Matrix_r(self,N,M) result(R)
 IMPLICIT none
 class(nummethod),intent(inout)::self
 INTEGER,intent(in)::N
 real*8,intent(in)::M(N,N)
 real*8::R
 !--------------------------------------------------------
 real*8::A(N,N)  ;INTEGER::IPIV(N),INFO ;INTEGER::JC
 A=M
 call DGETRF( n,	n, A, n, IPIV, INFO )
 R=1._8
 DO JC=1,INT(N,4)
   R=R*A(JC,JC)
   IF (.NOT. JC.EQ.IPIV(JC)) R=(-1._8)*R
 ENDDO
 END function
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


end module
