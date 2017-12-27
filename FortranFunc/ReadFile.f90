





!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! TYPE  : MODULE
! NAME  : ReadFile
! OBJECT: TYPE(readf)
! USED  :
! DATE  : 2017-11-4
! AUTHOR: hengyueli@gmail.com
!--------------
! Open-Source : No
!------------------
! DESCRIPTION:
!            action with file and data.
!
! STANDARD:
!
! USING LIST:
!
!
!
! avalable sets:
!                  [sub] Initialization(filepath)
!                        character(*)::filepath
!
!
! avalable gets:
!                  [sub] GetOneDataFromLine(linenumber,data)
!                        integer::linenumber
!                        integer/real*8/logical/characer(*)::data
!                        real a data from a line in the file
!
!
! avalable others:
!
!
!
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



module ReadFile
       implicit none



       type::readf
             private

             logical       ::Initiated = .false.
             character(len=:), allocatable::filepath

       contains
       procedure,pass::Initialization
       final::Finalization


       procedure,pass::GetOnelineString
       procedure,pass::GetLogicalFromLine
       procedure,pass::GetRealFromLine
       procedure,pass::GetIntegerFromLine
       generic::GetOneDataFromLine =>GetOnelineString,GetLogicalFromLine,GetRealFromLine,GetIntegerFromLine
       endtype


       private::Initialization,UnInitialization,Finalization

       private::GetOnelineString,GetLogicalFromLine,GetRealFromLine,GetIntegerFromLine


     contains

       subroutine Initialization(self,filepath)
                  implicit none
                  class(readf),intent(inout)::self
                  character(*),intent(in)::filepath
                  !--------------------
                  integer::l

                  call UnInitialization(self)
                  self%Initiated = .true.

                  l = len(trim(adjustl(filepath)))
                  allocate( character(len=l)::self%filepath )
                  self%filepath = trim(adjustl(filepath))

       endsubroutine


       subroutine UnInitialization(self)
         implicit none
         class(readf),intent(inout)::self
         !---------------------
         if (self%Initiated)then
           deallocate(self%filepath)
           self%Initiated = .false.
         endif
       endsubroutine


       subroutine Finalization(self)
         implicit none
         type(readf),intent(inout)::self
         call UnInitialization(self)
       endsubroutine



       !------------------------------------------------------------------------------------




       subroutine GetOnelineString(self,linenumber,string)
       implicit none
       class(readf),intent(inout)::self
       integer,intent(in)::linenumber
       character(*),intent(out)::string
       !----------------------------------
       character(30)::for   ;integer::jc
       write(for,*)len(string)
       for=trim("(A")//trim(ADJUSTL(for))//trim(")")
       open(9090,file=self%filepath,ACTION="READ")
        do jc=1,linenumber-1
           read(9090,*)
        end do
        read(9090,for)string
        close(9090)
       end subroutine


       subroutine GetIntegerFromLine(self,linenumber,returnInt)
           implicit none
           class(readf),intent(inout)::self
           integer,intent(in)::linenumber
           integer,intent(out)::returnInt
           !------------------------------
           character(256)::temp
           call GetOnelineString(self,linenumber,temp)
           read(temp,*)returnInt
       end subroutine

       subroutine GetRealFromLine(self,linenumber,returnReal)
         implicit none
         class(readf),intent(inout)::self
         integer,intent(in)::linenumber
         real*8,intent(out)::returnReal
         !------------------------------
         character(256)::temp
         call GetOnelineString(self,linenumber,temp)
         read(temp,*)returnReal
       end subroutine

       subroutine GetLogicalFromLine(self,linenumber,returnLogical)
         implicit none
         class(readf),intent(inout)::self
         integer,intent(in)::linenumber
         logical,intent(out)::returnLogical
         !------------------------------
         character(256)::temp
         call GetOnelineString(self,linenumber,temp)
         read(temp,*)returnLogical
       end subroutine







endmodule
