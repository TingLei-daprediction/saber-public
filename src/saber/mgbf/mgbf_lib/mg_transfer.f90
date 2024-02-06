!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        submodule(mg_intstate)  mg_transfer 
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
!TEST
!use mg_output, only: output_spec1_2dd
!TEST


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                       module  subroutine anal_to_filt_all(this,WORKA)
!***********************************************************************
!                                                                      !
!  Transfer data from analysis to first generaton of filter grid       !
!                                                                      !
!***********************************************************************
implicit none
class(mg_intstate_type),target::this
real (r_kind):: WORKA(this%km,this%n0:this%nm,this%m0:this%mm)

real(r_kind),allocatable,dimension(:,:,:):: VLOC  
include  "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include  "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"

!----------------------------------------------------------------------

    allocate(VLOC(km,i0-ib:im+ib,j0-jb:jm+jb))                      


!T                                                 call btim(  aintp_tim)

      VLOC=0.
         call this%lsqr_adjoint_offset(WORKA,VLOC,km)


!T                                                 call etim(  aintp_tim)


!***
!***  Apply adjoint lateral bc on PKF and WKF
!***
    

         call this%bocoT_2d(VLOC,km,im,jm,ib,jb)
 
       VALL=0.
       VALL(1:km,i0:im,j0:jm)=VLOC(1:km,i0:im,j0:jm)
      

    deallocate(VLOC)

!                                            call etim(   btrns1_tim)

!----------------------------------------------------------------------
                        endsubroutine anal_to_filt_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                       module  subroutine filt_to_anal_all(this,WORKA)
!***********************************************************************
!                                                                      !
!  Transfer data from filter to analysis grid                          !
!                                                                      !
!***********************************************************************
implicit none
class(mg_intstate_type),target::this
real (r_kind):: WORKA(this%km,this%n0:this%nm,this%m0:this%mm)


real(r_kind),allocatable,dimension(:,:,:):: VLOC   
include  "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include  "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!TEST
!real(r_kind), allocatable, dimension(:,:):: PA
!TEST

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
         call this%boco_2d(VLOC,km,im,jm,ib,jb)


!***
!*** Interpolate to analysis grid composite variables
!***
!TEST
!    allocate(PA(1:im,1:jm))
!
!     PA(1:im,1:jm)=VLOC(3*lm+lm/2,1:im,1:jm)
!
!        call output_spec1_2dd(PA,im,jm)
!
!     call finishMPI
!TEST


!T                                                 call btim(   intp_tim)

         call this%lsqr_direct_offset(VLOC,WORKA,this%km) !cltthink

!T                                                 call etim(   intp_tim)
    deallocate(VLOC)


!T                                                 call etim(   btrns2_tim)

!----------------------------------------------------------------------
                        endsubroutine filt_to_anal_all


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module subroutine stack_to_composite                   &
!***********************************************************************
!                                                                      !
!  Transfer data from stack to composite variables                     !
!                                                                      !
!***********************************************************************
(this,ARR_ALL,A2D,A3D)
!----------------------------------------------------------------------
implicit none
class(mg_intstate_type),target::this
real(r_kind),dimension(this%km ,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),   intent(in):: ARR_ALL
real(r_kind),dimension(this%km3,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy,this%lm),intent(out):: A3D
real(r_kind),dimension(this%km2,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy)   ,intent(out):: A2D
!----------------------------------------------------------------------
integer(i_kind)::i,j,k, L
include  "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include  "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
    do L=1,lm
      do j=j0-hy,jm+hy
      do i=i0-hx,im+hx
        do k=1,km3
          A3D(k,i,j,L)=ARR_ALL( (k-1)*lm+L,i,j )
        enddo
      enddo
      enddo
    enddo

        do k=1,km2
          A2D(k,:,:)=ARR_ALL(km3*lm+k,:,:)
        enddo 

!----------------------------------------------------------------------
                        endsubroutine stack_to_composite

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                       module  subroutine composite_to_stack                   &
!***********************************************************************
!                                                                      !
!  Transfer data from composite to stack variables                     !
!                                                                      !
!***********************************************************************
(this,A2D,A3D,ARR_ALL)
!----------------------------------------------------------------------
implicit none
class(mg_intstate_type),target::this
real(r_kind),dimension(this%km2,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),   intent(in):: A2D
real(r_kind),dimension(this%km3,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy,this%lm),intent(in):: A3D
real(r_kind),dimension(this%km ,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),   intent(out):: ARR_ALL
integer(i_kind):: i,j,k,L
include  "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include  "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!----------------------------------------------------------------------
    do L=1,lm
      do j=j0-hy,jm+hy
      do i=i0-hx,im+hx
        do k=1,km3
          ARR_ALL( (k-1)*lm+L,i,j )=A3D(k,i,j,L)
        enddo
      enddo
      enddo
    enddo

        do k=1,km2
          ARR_ALL(km3*lm+k,:,:)=A2D(k,:,:)
        enddo 

!----------------------------------------------------------------------
                        endsubroutine composite_to_stack 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        end submodule mg_transfer
