!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_transfer
!***********************************************************************
!                                                                      !
!  Transfer data between analysis and filter grid                      !
!                                                                      !
! Modules: kinds, mg_parameter, mg_intstate, mg_bocos, mg_interpolate, !
!          mg_timers, mg_mppstuff                                      !
!                                                     M. Rancic (2021) !
!***********************************************************************
use mpi
use kinds, only: r_kind,i_kind
use mg_parameter
use mg_intstate, only: VALL,WORKA
use mg_mppstuff, only:  mype,ierror,mpi_comm_world
use mg_mppstuff, only: nx,my,mpi_comm_comp

implicit none

integer(i_kind):: n,m,l,k,i,j

public anal_to_filt_all
public filt_to_anal_all

public stack_to_composite
public composite_to_stack
public 
type mg_transfer_type 
 contains
 procedure,nopass :: anal_to_filt_all
 procedure,nopass :: filt_to_anal_all 
end type mg_transfer_type 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine anal_to_filt_all
!***********************************************************************
!                                                                      !
!  Transfer data from analysis to first generaton of filter grid       !
!                                                                      !
!***********************************************************************
use mg_interpolate, only: lsqr_adjoint_offset
use mg_bocos, only:  bocoT_2d
implicit none

real(r_kind),allocatable,dimension(:,:,:):: VLOC  

!----------------------------------------------------------------------

    allocate(VLOC(km,i0-ib:im+ib,j0-jb:jm+jb))                      


!T                                                 call btim(  aintp_tim)

      VLOC=0.
         call lsqr_adjoint_offset(WORKA,VLOC,km)


!T                                                 call etim(  aintp_tim)


!***
!***  Apply adjoint lateral bc on PKF and WKF
!***
    

         call bocoT_2d(VLOC,km,im,jm,ib,jb)
 
       VALL=0.
       VALL(1:km,i0:im,j0:jm)=VLOC(1:km,i0:im,j0:jm)
      

    deallocate(VLOC)

!                                            call etim(   btrns1_tim)

!----------------------------------------------------------------------
                        endsubroutine anal_to_filt_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filt_to_anal_all
!***********************************************************************
!                                                                      !
!  Transfer data from filter to analysis grid                          !
!                                                                      !
!***********************************************************************
use mg_interpolate, only: lsqr_direct_offset
use mg_bocos, only:  boco_2d
implicit none


real(r_kind),allocatable,dimension(:,:,:):: VLOC   


!----------------------------------------------------------------------

!T                                            call btim(   btrns2_tim)

!***
!***  Define VLOC
!***

    allocate(VLOC(1:km,i0-ib:im+ib,j0-jb:jm+jb))                     

      VLOC=0.
      VLOC(1:km,i0:im,j0:jm)=VALL(1:km,i0:im,j0:jm)
        

!***
!***  Supply boundary conditions for VLOC
!***
         call boco_2d(VLOC,km,im,jm,ib,jb)


!***
!*** Interpolate to analysis grid composite variables
!***


!T                                                 call btim(   intp_tim)

         call lsqr_direct_offset(VLOC,WORKA,km)

!T                                                 call etim(   intp_tim)
    deallocate(VLOC)


!T                                                 call etim(   btrns2_tim)

!----------------------------------------------------------------------
                        endsubroutine filt_to_anal_all


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine stack_to_composite                   &
!***********************************************************************
!                                                                      !
!  Transfer data from stack to composite variables                     !
!                                                                      !
!***********************************************************************
(ARR_ALL,A2D,A3D)
!----------------------------------------------------------------------
implicit none
real(r_kind),dimension(km ,i0-hx:im+hx,j0-hy:jm+hy),   intent(in):: ARR_ALL
real(r_kind),dimension(km3,i0-hx:im+hx,j0-hy:jm+hy,lm),intent(out):: A3D
real(r_kind),dimension(km2,i0-hx:im+hx,j0-hy:jm+hy)   ,intent(out):: A2D
!----------------------------------------------------------------------
    do L=1,lm
      do j=j0-hy,jm+hy
      do i=i0-hx,im+hx
        A3D(1,i,j,L)=ARR_ALL(     L,i,j)
        A3D(2,i,j,L)=ARR_ALL(  lm+L,i,j) 
        A3D(3,i,j,L)=ARR_ALL(2*lm+L,i,j)
        A3D(4,i,j,L)=ARR_ALL(3*lm+L,i,j)
        A3D(5,i,j,L)=ARR_ALL(4*lm+L,i,j)
        A3D(6,i,j,L)=ARR_ALL(5*lm+L,i,j)
      enddo
      enddo
    enddo


    A2D(1,:,:)=ARR_ALL(6*lm+1,:,:)
    A2D(2,:,:)=ARR_ALL(6*lm+2,:,:)
    A2D(3,:,:)=ARR_ALL(6*lm+3,:,:)
    A2D(4,:,:)=ARR_ALL(6*lm+4,:,:)

!----------------------------------------------------------------------
                        endsubroutine stack_to_composite

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine composite_to_stack                   &
!***********************************************************************
!                                                                      !
!  Transfer data from composite to stack variables                     !
!                                                                      !
!***********************************************************************
(A2D,A3D,ARR_ALL)
!----------------------------------------------------------------------
implicit none
real(r_kind),dimension(km2,i0-hx:im+hx,j0-hy:jm+hy),   intent(in):: A2D
real(r_kind),dimension(km3,i0-hx:im+hx,j0-hy:jm+hy,lm),intent(in):: A3D
real(r_kind),dimension(km ,i0-hx:im+hx,j0-hy:jm+hy),   intent(out):: ARR_ALL
integer(i_kind):: i,j,L
!----------------------------------------------------------------------
    do L=1,lm
      do j=j0-hy,jm+hy
      do i=i0-hx,im+hx
        ARR_ALL(     L,i,j)= A3D(1,i,j,L)
        ARR_ALL(  lm+L,i,j)= A3D(2,i,j,L)
        ARR_ALL(2*lm+L,i,j)= A3D(3,i,j,L)
        ARR_ALL(3*lm+L,i,j)= A3D(4,i,j,L)
        ARR_ALL(4*lm+L,i,j)= A3D(5,i,j,L)
        ARR_ALL(5*lm+L,i,j)= A3D(6,i,j,L)
      enddo
      enddo
    enddo


    ARR_ALL(6*lm+1,:,:)= A2D(1,:,:)
    ARR_ALL(6*lm+2,:,:)= A2D(2,:,:)
    ARR_ALL(6*lm+3,:,:)= A2D(3,:,:)
    ARR_ALL(6*lm+4,:,:)= A2D(4,:,:)

!----------------------------------------------------------------------
                        endsubroutine composite_to_stack 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_transfer
