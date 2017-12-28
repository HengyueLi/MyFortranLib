



!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : fermion_table
! OBJECT: TYPE(table)
! USED  : basic_math_functions
! DATE  : 2017-11-26
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            table used for all the exact diagonalization method.( Both for ED and Lanczos mothord) . Fermion system only.
!            [[[   keep  ns<=15 ]]]
! STANDARD:
!            *call Initialization( ns ,  symmetry ,*print):
!                            integer::ns is the total sites of the system
!                            integer::symmetry  denote the different symmetry in the system.
!                                        symmetry = 0  :   we can only diagonalize it dirrectly, there is only one subspace
!                                        symmetry = 1  :   particle number N is a good qunatiy. The subspace is indexed by N
!                                        symmetry = 2  :   N and spin are good quantities. subspace is indexed by ( nup,ndown)
! USING LIST:
!            :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! avalable sets:
!
! avalable gets:
!                  [fun] get_ns():
!                        integer::get_ns
!                  [fun] get_symmetry()
!                        integer::get_symmetry
!                  [fun] get_nsub():
!                                   integer::get_nsub the number of subspaces.
!                  [fun] get_sub_d(i):
!                                   integer::get_sub_d,i
!                                   return the dimension of the i-th subspace.
!                  [fun] get_basis_to_sub_index(basis):
!                                   integer::get_basis_to_sub_index
!                                   integer*8::basis
!                  [fun] get_subid_from_basis(basis):
!                                   return integer
!                                   integer*8::basis
!                  [fun] get_subindex_to_basis(subid,index):
!                                   return integer*8
!                                   integer::subid,index
!                  [fun] get_subspace_marker(subid):
!                                   character(32)::get_subspace_marker    description of the subspace
!                  [sub] get_sub_mark_value(subid,returnv):
!                                   integer::returnv(2)
!                                   for (Nup,Ndown) case , return (nup,ndown)
!                                   for (N) case, return (N, -1)
!                                   for other case , return (-1,-1)
!                  [fun] get_subid_from_mark(mark):
!                                   integer::mark(2),get_subid_from_mark
!
! avalable  IS:
!                  [fun] Is_initiated()
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module fermion_table
   implicit none

   type,private::subsapce
       logical   :: initiated = .false.
       integer   :: mark(2)   ! mark good quantities
       integer   :: d
       integer*8,allocatable::bais(:) ! bais(d)  return the nature basis.
   contains
       procedure,pass::Initialization_subsapce
       procedure,pass::UnInitialization_subsapce
       final::Final_sub
   end type

   !
   !
   type::table
       private
       logical::Initiated=.false.
       integer::symmetry
       integer::ns
       integer::nsub
       type(subsapce),allocatable::st(:)
       integer  ::HD   ! DIMENSION OF (original) H
       integer  ,allocatable::TS(:)        !TS(0:HD-1) , to sub
       integer,allocatable::TSsubid(:)   !only save subid
       !---------------
       integer::print = 6
   contains
       procedure,pass::Initialization
       final::Finalization
   !
       procedure,pass::get_ns
       procedure,pass::get_symmetry
       procedure,pass::get_nsub
       procedure,pass::get_sub_d
       procedure,pass::get_subid_from_basis
       procedure,pass::get_basis_to_sub_index
       procedure,pass::get_subindex_to_basis
       procedure,pass::get_sub_mark_value
       procedure,pass::get_subid_from_mark
      procedure,pass::get_subspace_marker
   !
       procedure,pass::Is_initiated
   end type
   !
   !
   !
   private::Initialization_subsapce,UnInitialization_subsapce,Final_sub
   private::Initialization,UnInitialization,Finalization
   private::set_table_elements
   private::get_subid_from_mark
   private::get_ns,get_symmetry,get_nsub,get_sub_d,get_basis_to_sub_index
   private::get_subindex_to_basis,get_sub_mark_value,get_subspace_marker
   private::get_subid_from_basis

   private::checkHowmany1

   contains




   subroutine Initialization_subsapce(self,d)
              implicit none
              class(subsapce),intent(inout)::self
              integer,intent(in)::d
              !--------------------------
              call UnInitialization_subsapce(self)
              self%initiated = .true.
              self%d = d
              allocate(   self%bais(d)  )
   end subroutine
   subroutine UNInitialization_subsapce(self)
              implicit none
              class(subsapce),intent(inout)::self
              !--------------------------
              if (self%initiated)then
                deallocate( self%bais  )
                 self%initiated = .false.
              endif
   end subroutine
   subroutine Final_sub(self)
     implicit none
     type(subsapce),intent(inout)::self
     !--------------------------
     call UNInitialization_subsapce(self)
end subroutine



   subroutine Initialization(self,ns,symmetry,print)
     implicit none
     class(table),intent(inout)::self
     integer,intent(in)::ns,symmetry
     integer,intent(in),optional::print
     !-------------------------------------
     call UnInitialization(self)
     self%initiated = .true.
     self%ns        = ns
     self%symmetry  = symmetry
     self%HD        = 4**ns
     if (present(print)) self%print = print

     SELECT CASE(symmetry)
     CASE(0)
       self%nsub = 1
     CASE(1)
       self%nsub = 2 * ns + 1
     CASE(2)
       self%nsub = (ns + 1)**2
     CASE DEFAULT
       write(self%print,*)"ERROR: Unknow symmetry @ table";stop
     ENDSELECT
     ALLOCATE(  self%st( self%nsub          )  )
     allocate(  self%TS( 0 : SELF%HD - 1    )  )
     allocate(  SELF%TSsubid(0 : SELF%HD - 1)  )
     call set_table_elements(self)
   endsubroutine

   subroutine UNInitialization(self)
     implicit none
     class(table),intent(inout)::self
     !-------------------------------------
     if (self%initiated)then
       deallocate(   self%st      )
       deallocate(   self%TS      )
       deallocate(   self%TSsubid )
       self%initiated = .false.
     endif
   endsubroutine

   impure elemental SUBROUTINE Finalization(SELF)
     implicit none
     TYPE(table),intent(inout)::self
     !-------------------------------------
     CALL UNInitialization(self)
   ENDSUBROUTINE



  integer function checkHowmany1(n,from,to)
     implicit none
     integer*8,intent(in)::n
     integer,intent(in)::from,to
     !----------------
     integer::jc
     checkHowmany1 = 0
     do jc = from,to
        if (btest(n,jc)) checkHowmany1 = checkHowmany1 + 1
     enddo
  endfunction


   subroutine set_table_elements(self)
     use basic_math_functions
     implicit none
     class(table),intent(inout)::self
     !-------------------------------------
     integer*8::jcbasis,jb1,jb2
     integer::jc1,jc2,d,n1,ns,counting(self%nsub),con2,marker(2)
     type(bmathf)::f

     ns = self%ns

     counting = 0_8
     select case(self%symmetry)
     case(0)
       self%TSsubid    =  1
       self%st(1)%mark = -1
       call self%st(1)%Initialization_subsapce(self%HD)
       do jcbasis = 0_8 , self%HD - 1_8
          self%TS(jcbasis)    =  jcbasis + 1_8
          self%st(1)%bais(jcbasis + 1_8) = jcbasis
       enddo
     case(1)
       do jc1 = 0 , 2 * self%ns
          jb1 = 2_8 * self%ns
          jb2 = jc1 * 1_8
          d = int(f%bcombination(jb2,jb1))
          call self%st(jc1+1)%Initialization_subsapce(d)
          self%st(jc1+1)%mark(1) = jc1
          self%st(jc1+1)%mark(2) = -1
       enddo
       do jcbasis = 0_8 , self%HD - 1_8
          n1                               = checkHowmany1(jcbasis,0,2*ns-1) + 1!subid
          counting(n1)                     = counting(n1) + 1
          self%TS( jcbasis )               = counting(n1)
          self%TSsubid(jcbasis)            = n1
          self%st(n1)%bais( counting(n1) ) = jcbasis
       enddo
     case(2)
       con2 = 0
       do jc1  = 0, ns ;  do jc2 = 0 , ns
          d    = f%bcombination(jc1,ns) * f%bcombination(jc2,ns)
          con2 = con2 + 1
          call self%st(con2)%Initialization_subsapce(d)
          self%st(con2)%mark(1) = jc1
          self%st(con2)%mark(2) = jc2
       enddo           ;  enddo
       do jcbasis = 0_8 , self%HD - 1_8
          marker(1) = checkHowmany1(jcbasis,0 ,  ns-1 )
          marker(2) = checkHowmany1(jcbasis,ns,2*ns-1 )
          n1 = get_subid_from_mark(self,marker)
          counting(n1) = counting(n1) + 1
          self%TS( jcbasis )               = counting(n1)
          self%st(n1)%bais( counting(n1) ) = jcbasis
          self%TSsubid(jcbasis)            = n1
       enddo
     endselect
   endsubroutine


  !  integer function get_subid_from_mark(self,marker)
  !    implicit none
  !    class(table),intent(inout)::self
  !    integer,intent(in)::marker(2)
  !    !-------------------------------------
  !    select case(self%symmetry)
  !    case(0)
  !      get_subid_from_mark = 1
  !    case(1)
  !      get_subid_from_mark = marker(1) + 1
  !    case(2)
  !      get_subid_from_mark = marker(2) * (self%ns + 1) + marker(1) + 1
  !    endselect
  !  endfunction
  integer function get_subid_from_mark(self,marker)
    implicit none
    class(table),intent(inout)::self
    integer,intent(in)::marker(2)
    !-------------------------------------
    integer::jc

    get_subid_from_mark = -1
    do jc = 1 , self%nsub
       if ( (self%st(jc)%mark(1).eq.marker(1))  .and. (self%st(jc)%mark(2).eq.marker(2)))then
         get_subid_from_mark = jc
         goto 999
       endif
    enddo
999 continue
  endfunction







    integer    function get_ns(self)
               implicit none
               class(table),intent(inout)::self
               !-----------------------------------------
               get_ns = self%ns
    end function


    integer    function get_symmetry(self)
                implicit none
                class(table),intent(inout)::self
                !-----------------------------------------
                get_symmetry = self%symmetry
           end function

   integer function get_nsub(self)
              implicit none
              class(table),intent(inout)::self
              !-----------------------------------------
              get_nsub = self%nsub
   end function
   !
   integer function get_sub_d(self,i)
              implicit none
              class(table),intent(inout)::self
              integer,intent(in)::i
              !-----------------------------------------
              get_sub_d = self%st(i)%d
   end function

   integer function get_basis_to_sub_index(self,basis)
              implicit none
              class(table),intent(inout)::self
              integer*8,intent(in)::basis
              !-----------------------------------------
              get_basis_to_sub_index = self%TS(basis)
   end function
   !
   integer*8 function  get_subindex_to_basis(self,subid,jci)
              implicit none
              class(table),intent(inout)::self
              integer,intent(in)::subid,jci
              !-----------------------------------------
              get_subindex_to_basis = self%st(subid)%bais(jci)
   end function

   integer function get_subid_from_basis(self,basis)
       implicit none
       class(table),intent(inout)::self
       integer*8,intent(in)::basis
       !-----------------------------------------
       get_subid_from_basis = self%TSsubid(basis)
   endfunction

   subroutine get_sub_mark_value(self,subid,res) !result(res)
              implicit none
              class(table),intent(inout)::self
              integer,intent(in)::subid
              integer,intent(out)::res(2)
              !-----------------------------------------
              res = self%st(subid)%mark
   end subroutine

   character(32) function get_subspace_marker(self,subid)
              implicit none
              class(table),intent(inout)::self
              integer,intent(in)::subid
              !-----------------------------------------
              character(32)::v1,v2
              select case(self%symmetry)
                     case(0)
                        get_subspace_marker = "Original space"
                     case(1)
                        write(v1,*)self%st(subid)%mark(1)
                        get_subspace_marker = "---subspace:  N="//trim(adjustl(v1))//" ---"
                     case(2)
                        write(v1,*)self%st(subid)%mark(1)
                        write(v2,*)self%st(subid)%mark(2)
                        get_subspace_marker = "---subspace:  Nup= "//trim(adjustl(v1))//", Ndo= "&
                        //trim(adjustl(v2))//"-----"
                     case default
                        write(*,*)"error201706292147";stop
              end select
   end function


   logical function Is_initiated(self)
     implicit none
     class(table),intent(inout)::self
     !------------------------------------
     Is_initiated = self%initiated
   endfunction


end module
