!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_mppstuff
!***********************************************************************
!                                                                      !
!    Everything related to mpi communication                           !
!                                                                      !
! Library: mpi                                                         !
! Modules: kinds, mg_parameter                                         !
!                                                     M. Rancic (2020) !
!***********************************************************************
use mpi
use kinds, only: i_kind
use mg_parameter
implicit none

character(len=5):: c_mype
integer(i_kind):: mype
integer(i_kind):: npes,iTYPE,rTYPE,dTYPE,mpi_comm_comp,ierr,ierror
integer(i_kind):: mpi_comm_work,group_world,group_work
integer(i_kind):: mype_gr,npes_gr

integer(i_kind) my_hgen
integer(i_kind) mype_hgen
logical:: l_hgen
integer(i_kind):: nx,my
!keep_for_now integer(i_kind):: ns,ms,ninc,minc,ninc2,minc2

type  mppstuff_type
contains
procedure,nopass ::  init_mg_MPI,finishMPI,barrierMPI
end type  mppstuff_type

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine init_mg_MPI
!***********************************************************************
!                                                                      !
!     Initialize mpi                                                   !
!     Create group for filter grid                                     !
!                                                                      !
!***********************************************************************
use mpi


implicit none
integer(i_kind):: g,m
integer(i_kind), dimension(npes_filt):: out_ranks
integer(i_kind):: nf
!-----------------------------------------------------------------------
!clt mgbf4jedi
           mpi_comm_comp=MPI_COMM_WORLD
!***
!***  Initial MPI calls
!***
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(mpi_comm_comp,mype,ierr)
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
!  Associate a group with communicator mpi_comm_comp
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
!     call MPI_FINALIZE(mpi_comm_comp)
!     stop
!TEST
    
     

!-----------------------------------------------------------------------
!
      call MPI_BARRIER(mpi_comm_comp,ierr)
!
!-----------------------------------------------------------------------
                        endsubroutine init_mg_MPI

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine barrierMPI
!***********************************************************************
!                                                                      !
!     Call barrier for all                                             !
!                                                                      !
!***********************************************************************
use mpi

implicit none
integer:: ierr
!-----------------------------------------------------------------------

      call MPI_BARRIER(mpi_comm_comp,ierr)

!-----------------------------------------------------------------------
                        endsubroutine barrierMPI

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine finishMPI
!***********************************************************************
!                                                                      !
!     Finalize MPI                                                     !
!                                                                      !
!***********************************************************************
use mpi

implicit none
integer:: ierr

!-----------------------------------------------------------------------
!
      call MPI_FINALIZE(ierr)
      stop
!
!-----------------------------------------------------------------------
                        endsubroutine finishMPI

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_mppstuff

