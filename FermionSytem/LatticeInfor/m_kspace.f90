



!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : class_creat_k_space
! OBJECT: TYPE(kspace)
! USED  : class_integration
! DATE  : 2017-05-10
! AUTHOR: Hengyue Li
!--------------
! DESCRIPTION:
!            creat the reciprocal space of a lattice system.
!            step1 *call self%Initialization(Dim,n1,n2,n3) :  if Dim<3, say Dim=1, n2, n3 can be any value
!            setp2 *call self%set_ai(i,v(3)) where i=1,2,3 to set the position of lattice vector
!            setp3 call self%creatk(meshtype) to creat the Brillouin zone
!                  meshtype=0 :creat equal distant k points
!                  meshtype=1 :creat GS k points
!                  meshtype=3 :consider the symmetry and still use GS
!avalable sets:
!               set_print(integer(8)::whereto)
!
!avalable gets:
!               fun getk(jc) = real(8)::k(3)     return jc-th k point
!               fun getw(jc) = real(8)::w        return the corresponding weight
!               fun getnk() =integer(8)::nk      return the total number of k points
!               fun geta(i) = real(8)::a(3)      return vectors
!               fun getb(jc=1,2,3)=real(8)::a(3) return the reciprocal vector
!               fun getvolr=real(8)::a           return the volue of unit cell
!               fun getvolb=real(8)::a           return the volue of Brillouin zone
!               fun getd
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
module class_creat_k_space

     implicit none
     !------------------------

     !-------------
     real(8),parameter,private::pi=3.141592653589793238462643383279_8
     !-------------
     type kspace
        logical,private::Initialized=.false.
        integer(8),private::d  ! dimension of lattice
        real(8),private::ai(3,3) ! vector in real space  ai=a(:,i)
        logical,private::aisetted(3)=.false.
        logical,private::is_k_created=.false.
        integer(8),private::ni(3) !The number of point on each direction
        integer(8),private::meshtype !to be setted
        !------------------------------------------------
        integer(8),private::nk       !\prod_i n_i : the total number of k point
        real(8),allocatable,private::kpoint(:,:)
        real(8),allocatable,private::kweight(:)
        real(8),private::b(3,3)      !reprocal vectors   bi=b(:,i)
        !------------------------------------------------
        integer(8),private::whereprint=6_8
     contains
        procedure,pass::Initialization
        !-------------------------
        procedure,pass::set_ai
        procedure,pass::set_print
        !-------------------------
        procedure,pass::getd
        procedure,pass::getk
        procedure,pass::getw
        procedure,pass::getnk
        procedure,pass::geta
        procedure,pass::getb   !reciprocal vector
        procedure,pass::getvolr  !volue of unit cell
        procedure,pass::getvolb  !volue of Brillouin zone
        !-------------------------
        procedure,pass::creatk
        Final::finalization
     end type


     private::Initialization,UnInitialization,finalization
     private::check_all_needed_ai_are_setted,vector3D_multiply,creatk_12,cala_volue
     private::getk,getw,getnk,geta,getb,getvolr,getvolb,getd
     contains

        subroutine Initialization(self,D,n1,n2,n3)
            implicit none
            class(kspace),intent(inout)::self
            integer(8),intent(in)::D,n1,n2,n3
            !-----------------------------------
            integer(8)::ni(3),jc
            call UnInitialization(self)
            self%Initialized=.True.
            if (.NOT.(  (D.eq.1) .or. (D.eq.2) .or. (D.eq.3)))then
                write(self%whereprint,*)"input D=", D ," in K space module is illegal"
                stop
            end if
            self%d=D
            self%ni(1)=1_8;self%ni(2)=1_8;self%ni(3)=1_8
            ni(1)=n1  ;ni(2)=n2  ;ni(3)=n3
            do jc=1_8,D
                call check_ni(jc,ni(jc))
                self%ni(jc)=ni(jc)
            end do
            self%nk=1_8
            do jc=1_8,D
               self%nk=self%nk*self%ni(jc)
            end do
            allocate(self%kpoint(3,self%nk))
            allocate(self%kweight(self%nk))
            self%ai(1,1)=1._8 ;self%ai(2,1)=0._8  ;self%ai(3,1)=0._8
            self%ai(1,2)=0._8 ;self%ai(2,2)=1._8  ;self%ai(3,2)=0._8
            self%ai(1,3)=0._8 ;self%ai(2,3)=0._8  ;self%ai(3,3)=1._8
            contains
            subroutine check_ni(i,ni)
                implicit none
                integer(8)::i,ni
                !--------------
                if (ni.le.0)then
                    write(self%whereprint,*)"the input ni(i=)",i , "is illegal in k module setting";stop
                end if
            endsubroutine
        end subroutine

        subroutine UnInitialization(self)
            implicit none
            class(kspace),intent(inout)::self
            !--------------------------------------
            if (self%Initialized)then
              if ( allocated(self%kpoint)  )   deallocate(self%kpoint)
              if ( allocated(self%kweight)  )  deallocate(self%kweight)
              self%Initialized=.false.
            endif
        end subroutine

        subroutine finalization(self)
            implicit none
            type(kspace),intent(inout)::self
            !--------------------------------------
            call UnInitialization(self)
        end subroutine


        integer function getd(self)
          implicit none
          class(kspace),intent(inout)::self
          !----------------------------------
          getd = int(self%d)
        endfunction



        subroutine set_ai(self,i,ai)
            implicit none
            class(kspace),intent(inout)::self
            integer(8),intent(in)::i
            real(8),intent(in)::ai(3)
            if (i>self%d)then
                write(self%whereprint,*)"for that i>D, ai is not nessisary to be set. This operation is ignored"
            elseif(i<=0)then
                write(self%whereprint,*)"Error:illegal i in set_ai!!!";stop
            else
                self%ai(:,i)=ai
                self%aisetted(i)=.true.
            end if
        end subroutine

        subroutine set_print(self,whereto)
            implicit none
            class(kspace),intent(inout)::self
            integer(8),intent(in)::whereto
            !-----------------------------------
            self%whereprint=whereto
        end subroutine
        logical function check_all_needed_ai_are_setted(self)
            implicit none
            type(kspace),intent(in)::self
            !--------------------------------
            integer(8)::jc
            do jc=1_8,self%d
               if (.not.(self%aisetted(jc))) then
                    check_all_needed_ai_are_setted=.false.
                    goto 999
               end if
            end do
            check_all_needed_ai_are_setted=.true.
        999 continue
        end function

        subroutine creatk(self,meshtype)
            implicit none
            class(kspace),intent(inout)::self
            integer(8),intent(in)::meshtype
            !----------------------------------
            select case(meshtype)
               case(0,1)
                 self%meshtype=meshtype
                 call creatk_12(self)
                 self%is_k_created=.true.
               case default
                 write(self%whereprint,*)"meshtype error!!";stop
            end select
        end subroutine


       subroutine creatk_12(self)
           use class_integration
           implicit none
           type(kspace),intent(inout)::self
           !--------------------------------
           real(8)::vol,tempV(3),tempV2(3),mid(3),G(3,3),J,lenth
           real(8),allocatable::x(:,:),w(:,:)
           type(Integration)::Intlist
           integer(8)::jc,nmax,c(3),jcd
           !--------------------------------
           call vector3D_multiply(self%ai(:,1),self%ai(:,2),tempV)
           Vol=DABS(SUM(tempV*self%ai(:,3)))
           ! calculate reciprocal vector G1,G2 and G3
           call vector3D_multiply(self%ai(:,2),self%ai(:,3),G(:,1)) ;G(:,1)=G(:,1)*Pi*2._8/Vol  ;self%b(:,1)=G(:,1)
           call vector3D_multiply(self%ai(:,3),self%ai(:,1),G(:,2)) ;G(:,2)=G(:,2)*Pi*2._8/Vol  ;self%b(:,2)=G(:,2)
           call vector3D_multiply(self%ai(:,1),self%ai(:,2),G(:,3)) ;G(:,3)=G(:,3)*Pi*2._8/Vol  ;self%b(:,3)=G(:,3)
           !-----------------------
           !--calculate  J
           tempV = self%b(:,1)/dsqrt(sum(self%b(:,1)**2._8))
           tempV2= self%b(:,2)/dsqrt(sum(self%b(:,2)**2._8))
           call vector3D_multiply(tempV,tempV2,mid)
           J = sum (  mid *   self%b(:,3)/dsqrt(sum(self%b(:,3)**2._8))  )
           !-----------------------
           nmax=max(self%ni(1),self%ni(2),self%ni(3))
           allocate(x(nmax,3),w(nmax,3))
           do jc=1_8,3_8
              call Intlist%Initialization(self%ni(jc),  self%meshtype )
              call Intlist%get_xw(-1._8,1._8,x(1_8:self%ni(jc),jc),w(1_8:self%ni(jc),jc))
           end do
           !-------------------
           c(1)=0_8;c(2)=1_8;c(3)=1_8
           DO JC=1_8,self%nk
            !----------------------
               c(1)=c(1)+1_8
               if (c(1).gt.self%ni(1))then
                 c(2)=c(2)+1_8
                 if (c(2).gt.self%ni(2))then
                    c(3)=c(3)+1_8
                    c(2)=1_8
                 end if
                 c(1)=1_8
               end if
            !---------(c1,c2,c3)----------
            self%kweight(  jc)=1._8
            self%kpoint (:,jc)=0._8
            do jcd=1_8,self%d
                lenth=dsqrt(sum(G(:,JCd)**2))
                self%kpoint(:,jc)=self%kpoint(:,jc)+x(c(jcd),jcd)*G(:,JCd)/2._8
                self%kweight( jc)=self%kweight( jc)*w(c(jcd),jcd)*lenth   /2._8
            end do
            self%kweight(  jc) = self%kweight(  jc) * J
           END DO
           self%kweight=self%kweight/(2._8*pi)**real(self%d,8)*Vol
           deallocate(x,w)
       end subroutine

        subroutine vector3D_multiply(a,b,c)
            implicit none
            real(8)::a(3),b(3),c(3)
            c(1)=a(2)*b(3)-a(3)*b(2)
            c(2)=a(3)*b(1)-a(1)*b(3)
            c(3)=a(1)*b(2)-a(2)*b(1)                    ! ;write(*,*)a  ;write(*,*)b ;write(*,*)c ;write(*,*)"---------"
        end subroutine

        function getk(self,jc) result(res)
               implicit none
               class(kspace),intent(inout)::self
               integer(8),intent(in)::jc
               real(8)::res(3)
               !------------------
               if (self%is_k_created) then
                   if ( (jc>0)  .and.  (jc<=self%nk)  )then
                        res=self%kpoint(:,jc)
                    else
                        write(self%whereprint,*)"illeagal index input in getk"
                   end if
               else
                   write(self%whereprint,*)"k space is not created!";stop
               end if
        end function

        real(8) function getw(self,jc)
               implicit none
               class(kspace),intent(inout)::self
               integer(8),intent(in)::jc
               !------------------
               if (self%is_k_created) then
                   if ( (jc>0)  .and.  (jc<=self%nk)  )then
                        getw=self%kweight(jc)
                    else
                        write(self%whereprint,*)"illeagal index input in getk";stop
                   end if
               else
                   write(self%whereprint,*)"k space is not created!";stop
               end if
        end function

        integer(8) function getnk(self)
               implicit none
               class(kspace),intent(inout)::self
               !-------------------------------------
               if (self%Initialized)then
                  getnk=self%nk
               else
                  write(self%whereprint,*)"kspace is not created, get nk?";stop
               end if
        end function

        function geta(self,jc) result(res)
               implicit none
               class(kspace),intent(inout)::self
               integer(8),intent(in)::jc
               real(8)::res(3)
               !------------------
               if (self%is_k_created) then
                   if ( (jc>0)  .and.  (jc<=self%d)  )then
                        res=self%ai(:,jc)
                    elseif(jc<=3)then
                        write(self%whereprint,*)&
                        "illeagal index input in geta. Notice the dimension of the system is:"&
                        ,self%d
                        res=self%ai(:,jc)
                    else
                        write(self%whereprint,*)"illeagal index input in geta"
                        stop
                   end if
               else
                   write(self%whereprint,*)"k space is not created!";stop
               end if
        end function


        function getb(self,jc) result(res)
               implicit none
               class(kspace),intent(inout)::self
               integer(8),intent(in)::jc
               real(8)::res(3)
               !------------------
               if (self%is_k_created) then
                   if ( (jc>0)  .and.  (jc<=self%d)  )then
                        res=self%b(:,jc)
                    elseif(jc<=3)then
                        write(self%whereprint,*)&
                        "illeagal index input in getb. Notice the dimension of the system is:"&
                        ,self%d
                        res=self%b(:,jc)
                    else
                        write(self%whereprint,*)"illeagal index input in getb";stop
                   end if
               else
                   write(self%whereprint,*)"k space is not created!";stop
               end if
        end function

       real(8) function cala_volue(a1,a2,a3)
               implicit none
               real(8),intent(in)::a1(3),a2(3),a3(3)
               !--------------------------
               real(8)::temp(3)
               call vector3D_multiply(a1,a2,temp)
               cala_volue=dabs(sum( a3*temp   ))
       end function



        real(8) function getvolr(self)
               implicit none
               class(kspace),intent(inout)::self
               !--------------------------------------
               getvolr=cala_volue(self%ai(:,1),self%ai(:,2),self%ai(:,3))
        end function

        real(8) function getvolb(self)
               implicit none
               class(kspace),intent(inout)::self
               !--------------------------------------
               integer(8)::jc
               getvolb=cala_volue(self%b(:,1),self%b(:,2),self%b(:,3))
               do jc = 1_8,3_8-self%d
                  getvolb=getvolb/(2._8*pi)
               end do
        end function



end module
