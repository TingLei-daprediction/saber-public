!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        submodule(mg_parameter)  mg_mppstuff
!***********************************************************************
!                                                                      !
!    Everything related to mpi communication                           !
!                                                                      !
! Library: mpi                                                         !
! Modules: kinds, mg_parameter                                         !
!                                                     M. Rancic (2020) !
!***********************************************************************
use kinds, only: i_kind
implicit none


!keep_for_now integer(i_kind):: ns,ms,ninc,minc,ninc2,minc2



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module subroutine init_mg_MPI(this)
!***********************************************************************
!                                                                      !
!     Initialize mpi                                                   !
!     Create group for filter grid                                     !
!                                                                      !
!***********************************************************************
use mpi


implicit none
class (mg_parameter_type),target:: this
integer(i_kind):: g,m
integer(i_kind), dimension(this%npes_filt):: out_ranks
integer(i_kind):: nf
integer(i_kind)::ierr
integer(i_kind):: color
include  "type_parameter_locpointer.inc"
include  "type_parameter_point2this.inc"
!-----------------------------------------------------------------------

!cltorg           mpi_comm_comp=MPI_COMM_WORLD
!***
!***  Initial MPI calls
!***
!cltorg      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,npes,ierr)
!      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      ! Create a new communicator with MPI_Comm_split
      color=1  ! just create an communicator now for the whole processes
      write(6,*)'thinkdebmype is ',mype
      call MPI_Comm_split(MPI_COMM_WORLD, color, mype, mpi_comm_comp, ierr)
      call MPI_COMM_SIZE(mpi_comm_comp,npes,ierr)





      rTYPE = MPI_REAL
      dTYPE = MPI_DOUBLE
      iTYPE = MPI_INTEGER


!***
!*** Analysis grid
!***

    nx = mod(mype,nxm)+1
    my = (mype/nxm)+1

!    if(nx==1) then
!       ns=0
!       ninc=1
!       ninc2=2
!    else 
!       ns=1
!       ninc=0
!       ninc2=1
!    endif
!
!    if(my==1) then
!       ms=0
!       minc=1
!       minc2=2
!    else 
!       ms=1
!       minc=0
!       minc2=1
!    endif


!***
!***  Define PEs that handle high generations
!***

   
      mype_hgen=-1
      my_hgen=-1

      if( mype < maxpe_filt-nxy(1)) then
        mype_hgen=mype+nxy(1)
      endif
      do g=1,gm
        if(maxpe_fgen(g-1)<= mype_hgen .and. mype_hgen< maxpe_fgen(g)) then
            my_hgen=g
         endif
      enddo
      l_hgen = mype_hgen >-1

!TEST
!      write(300+mype,*)'mype,my_hgen,l_gen,mype_hgen=',mype,my_hgen,l_hgen,mype_hgen
!TEST

!***
!***  Chars
!***
      write(c_mype,1000) mype
 1000 format(i5.5)


!-----------------------------------------------------------------------
!
      call MPI_BARRIER(mpi_comm_comp,ierr)
!
!-----------------------------------------------------------------------
!***
!***  Define group communicator for higher generations
!***
!
!  Associate a group with communicator this@mpi_comm_comp
!
      call MPI_COMM_GROUP(mpi_comm_comp,group_world,ierr)
!
!  Create a new group out of exising group
!
     do nf = 1,npes_filt
       out_ranks(nf)=nf-1
     enddo 

     call MPI_GROUP_INCL(group_world,npes_filt,out_ranks,group_work,ierr)
!
!  Now create a new communicator associated with new group
!    
     call MPI_COMM_CREATE(mpi_comm_comp, group_work, mpi_comm_work, ierr)

    if( mype < npes_filt) then


      call MPI_COMM_RANK(mpi_comm_work,mype_gr,ierr)
      call MPI_COMM_SIZE(mpi_comm_work,npes_gr,ierr)

   else
       
      mype_gr= -1
      npes_gr= npes_filt
 
   endif 

!TEST
!     write(mype+100,*) 'mype, mype_gr=',mype, mype_gr
!     print *, 'mype, mype_gr=',mype, mype_gr
!     call MPI_FINALIZE(this@mpi_comm_comp)
!     stop
!TEST
    
     

!-----------------------------------------------------------------------
!
      call MPI_BARRIER(mpi_comm_comp,ierr)
!
!-----------------------------------------------------------------------
                        endsubroutine init_mg_MPI

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                  module      subroutine barrierMPI(this)
!***********************************************************************
!                                                                      !
!     Call barrier for all                                             !
!                                                                      !
!***********************************************************************
use mpi

implicit none
        class(mg_parameter_type),target::this
integer(i_kind):: ierr
include  "type_parameter_locpointer.inc"
include  "type_parameter_point2this.inc"
!-----------------------------------------------------------------------

      call MPI_BARRIER(mpi_comm_comp,ierr)

!-----------------------------------------------------------------------
                        endsubroutine barrierMPI

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                       module  subroutine finishMPI(this)
!***********************************************************************
!                                                                      !
!     Finalize MPI                                                     !
!                                                                      !
!***********************************************************************
use mpi

implicit none
        class(mg_parameter_type),target::this
!cltthinkdeb don't need mpi_finalize if mgbf is a lib to be called from outside
!
      call MPI_FINALIZE(this%ierr)
      stop
!
!-----------------------------------------------------------------------
                        endsubroutine finishMPI

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        end submodule mg_mppstuff

