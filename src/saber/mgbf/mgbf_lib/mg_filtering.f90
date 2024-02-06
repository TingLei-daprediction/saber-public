!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        submodule(mg_intstate) mg_filtering
!***********************************************************************
!                                                                      !
! Contains all multigrid filtering prodecures                          ! 
!                                                                      ! 
!                                                     M. Rancic (2020) !
!***********************************************************************
use mg_timers
use kinds, only: r_kind,i_kind
!clt use jp_pbfil, only: rbeta,rbetaT
use jp_pbfil3, only: dibetat,dibeta
use mpi




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
   module                  subroutine mg_filtering_procedure(this,mg_filt) 
!***********************************************************************
!                                                                      !
! Driver for Multigrid filtering procedures with Helmholtz operator    !
!                                                                      !
!   1, 2, 3: Radial filter                                             !
!         1: 2d radial filter for all variables                        !
! ->      2: 2d radial filter with 1d in vertical for 3d variables     !
!         3: 3d radial filter for 3d variables                         !
!                                                                      !
!   4, 5, 6: Line filter                                               !
!         4: 2d line filter for all variables                          !
!         5: 2d line filter with 1d in vertical for 3d variables       !
!         6: 3d line filter for 3d variables                           !
!                                                                      !
!                                                                      !
!***********************************************************************
implicit none 
class(mg_intstate_type),target::this

integer(i_kind),intent(in):: mg_filt
include  "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include  "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------
  if(mgbf_line) then
    if(mg_filt<4) then
       print*,'("Line filters have options 4-6")'
       stop
    endif
  else
    if(mg_filt>3) then
       print*,'("Radial filters have options 1-3")'
       stop
    endif
  endif 
      select case(mg_filt)
        case(1)
          call this%mg_filtering_rad1
        case(2)
          call this%mg_filtering_rad2
        case(3)
          call this%mg_filtering_rad3
        case(4)
          call this%mg_filtering_lin1
        case(5)
          call this%mg_filtering_lin2
        case(6)
          call this%mg_filtering_lin3
        case default
          call this%mg_filtering_fast          
       end select

!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_procedure    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
           module             subroutine mg_filtering_rad1(this)
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 1:                                     !
!                                                                      !
!     - Multiple of 2D and 3D variables                                !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 2d radial filter only for all variables                        !
!                                                                      !
!***********************************************************************
implicit none
class(mg_intstate_type),target:: this

integer(i_kind) L,i,j,g
include  "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include  "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------


!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend (Step 1)
!***
     
                                                 call btim( upsend_tim)
       call this%upsending_all(VALL,HALL,lquart)
                                                 call etim( upsend_tim)
!----------------------------------------------------------------------


!----------------------------------------------------------------------

!***
!*** Apply adjoint of Beta filter at all generations 
!***
                                                 call btim(    bfiltT_tim)



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

      call this%rbetaT(km,hx,i0,im,hy,j0,jm,pasp2,ss2,VALL(:,:,:))
  if(l_hgen)  then
      call this%rbetaT(km,hx,i0,im,hy,j0,jm,pasp2,ss2,HALL(:,:,:))
  endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

        write(6,*)'thinkdeb33 1 ', km,im,jm,hx,hy 
        call this%bocoT_2d(VALL,km,im,jm,hx,hy)
        call this%bocoT_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)


                                                 call etim(    bfiltT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

      call this%weighting_all(VALL,HALL,lhelm)


                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations
!***

                                                 call btim( bfilt_tim)

      call this%boco_2d(VALL,km,im,jm,hx,hy)
      call this%boco_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!

      call this%rbeta(km,hx,i0,im,hy,j0,jm,pasp2,ss2,VALL(:,:,:))
  if(l_hgen)  then
      call this%rbeta(km,hx,i0,im,hy,j0,jm,pasp2,ss2,HALL(:,:,:))
  endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfilt_tim)

!***
!***  Downsend, interpolate and add, then zero high generations 
!***

                                                 call btim(   dnsend_tim)
       call this%downsending_all(HALL,VALL,lquart)

                                                 call etim(   dnsend_tim)


!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_rad1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
             module           subroutine mg_filtering_rad2(this)
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 2:                                     !
!                                                                      !
!     - Multiple of 2D and 3D variables                                !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 2d radial filter + 1d vertical filter                          !
!                                                                      !
!***********************************************************************
implicit none
class (mg_intstate_type),target::this

real(r_kind), allocatable, dimension(:,:,:):: VM2D
real(r_kind), allocatable, dimension(:,:,:):: HM2D
real(r_kind), allocatable, dimension(:,:,:,:):: VM3D
real(r_kind), allocatable, dimension(:,:,:,:):: HM3D

integer(i_kind) L,i,j
include  "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include  "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------

allocate(VM3D(km3,i0-hx:im+hx,j0-hy:jm+hy,lm))                  ; VM3D=0.
allocate(VM2D(km2,i0-hx:im+hx,j0-hy:jm+hy   ))                  ; VM2D=0.
allocate(HM3D(km3,i0-hx:im+hx,j0-hy:jm+hy,lm))                  ; HM3D=0.
allocate(HM2D(km2,i0-hx:im+hx,j0-hy:jm+hy   ))                  ; HM2D=0.



!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend 
!***
     
                                                 call btim( upsend_tim)
       call this%upsending_all(VALL,HALL,lquart)
                                                 call etim( upsend_tim)
!----------------------------------------------------------------------


!----------------------------------------------------------------------

!***
!*** Apply adjoint of Beta filter at all generations 
!***
                                                 call btim(    bfiltT_tim)



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

      call this%rbetaT(km,hx,i0,im,hy,j0,jm,pasp2,ss2,VALL(:,:,:))
      call this%stack_to_composite(VALL,VM2D,VM3D)

  if(l_hgen)  then
      call this%rbetaT(km,hx,i0,im,hy,j0,jm,pasp2,ss2,HALL(:,:,:))
      call this%stack_to_composite(HALL,HM2D,HM3D)
  endif

      call this%sup_vrbeta1T(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call this%composite_to_stack(VM2D,VM3D,VALL)
  if(l_hgen)  then
      call this%sup_vrbeta1T(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,HM3D)
      call this%composite_to_stack(HM2D,HM3D,HALL)
   endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


        write(6,*)'thinkdeb33 2 ', km,im,jm,hx,hy 
        call this%bocoT_2d(VALL,km,im,jm,hx,hy)
        call this%bocoT_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)


                                                 call etim(    bfiltT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

      call this%weighting_all(VALL,HALL,lhelm)


                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations (Step 7)
!***
                                                 call btim( bfilt_tim)

      call this%boco_2d(VALL,km,im,jm,hx,hy)
      call this%boco_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!

      call this%rbeta(km,hx,i0,im,hy,j0,jm,pasp2,ss2,VALL(:,:,:))
      call this%stack_to_composite(VALL,VM2D,VM3D)
  if(this%l_hgen)  then
      call this%rbeta(km,hx,i0,im,hy,j0,jm,pasp2,ss2,HALL(:,:,:))
      call this%stack_to_composite(HALL,HM2D,HM3D)
  endif

      call this%sup_vrbeta1(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call this%composite_to_stack(VM2D,VM3D,VALL)
  if(l_hgen)  then
      call this%sup_vrbeta1(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,HM3D)
      call this%composite_to_stack(HM2D,HM3D,HALL)
   endif
       call this%barrierMPI



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfilt_tim)

!***
!***  Downsend, interpolate and add (Step 4)
!***  Then zero high generations (Step 5)
!***

                                                 call btim(   dnsend_tim)
       call this%downsending_all(HALL,VALL,lquart)

                                                 call etim(   dnsend_tim)

deallocate(VM3D) 
deallocate(VM2D)
deallocate(HM3D)
deallocate(HM2D)

!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_rad2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          module              subroutine mg_filtering_rad3(this)
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 2:                                     !
!                                                                      !
!     - Multiple of 2D and 3D variables                                !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 3d radial filter 
!                                                                      !
!***********************************************************************
!-----------------------------------------------------------------------
implicit none
class (mg_intstate_type),target::this


real(r_kind), allocatable, dimension(:,:,:):: VM2D
real(r_kind), allocatable, dimension(:,:,:):: HM2D
real(r_kind), allocatable, dimension(:,:,:,:):: VM3D
real(r_kind), allocatable, dimension(:,:,:,:):: HM3D


integer(i_kind) L,i,j
include  "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include  "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"

!----------------------------------------------------------------------
allocate(VM3D(km3,i0-hx:im+hx,j0-hy:jm+hy,lm))                 ; VM3D=0.
allocate(VM2D(km2,i0-hx:im+hx,j0-hy:jm+hy   ))                 ; VM2D=0.
allocate(HM3D(km3,i0-hx:im+hx,j0-hy:jm+hy,lm))                 ; HM3D=0.
allocate(HM2D(km2,i0-hx:im+hx,j0-hy:jm+hy   ))                 ; HM2D=0.

!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend 
!***
     
                                                 call btim( upsend_tim)
       call this%upsending_all(VALL,HALL,lquart)
                                                 call etim( upsend_tim)


!***
!*** Apply adjoint of Beta filter at all generations 
!***
                                                 call btim(    bfiltT_tim)



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Adjoint filtering
!
      call this%stack_to_composite(VALL,VM2D,VM3D)
        call this%rbetaT(km2,hx,i0,im,hy,j0,jm,pasp2,ss2,VM2D)
        call this%sup_vrbeta3T(km3,hx,hy,hz,im,jm,lm,pasp3,ss3,VM3D)
      call this%composite_to_stack(VM2D,VM3D,VALL)

    if(l_hgen) then
      call this%stack_to_composite(HALL,HM2D,HM3D)
        call this%rbetaT(km2,hx,i0,im,hy,j0,jm,pasp2,ss2,HM2D)
        call this%sup_vrbeta3T(km3,hx,hy,hz,im,jm,lm,pasp3,ss3,HM3D)
      call this%composite_to_stack(HM2D,HM3D,HALL)
    endif 

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
        write(6,*)'thinkdeb33 3 ', km,im,jm,hx,hy 

        call this%bocoT_2d(VALL,km,im,jm,hx,hy)
        call this%bocoT_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)


                                                 call etim(    bfiltT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                 call btim( weight_tim)

      call this%weighting_all(VALL,HALL,lhelm)

                                                 call etim( weight_tim)


!***
!*** Apply Beta filter at all generations 
!***

                                                 call btim( bfilt_tim)

      call this%boco_2d(VALL,km,im,jm,hx,hy)
      call this%boco_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!
      call this%stack_to_composite(VALL,VM2D,VM3D)
        call this%rbeta(km2,hx,i0,im,hy,j0,jm,pasp2,ss2,VM2D(:,:,:))
        call this%sup_vrbeta3(km3,hx,hy,hz,im,jm,lm,pasp3,ss3,VM3D)
      call this%composite_to_stack(VM2D,VM3D,VALL)
  if(l_hgen)  then
      call this%stack_to_composite(HALL,HM2D,HM3D)
        call this%rbeta(km2,hx,i0,im,hy,j0,jm,pasp2,ss2,HM2D(:,:,:))
        call this%sup_vrbeta3(km3,hx,hy,hz,im,jm,lm,pasp3,ss3,HM3D)
      call this%composite_to_stack(HM2D,HM3D,HALL)
  endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfilt_tim)

!***
!***  Downsend, interpolate and add 
!***  Then zero high generations 
!***

                                                 call btim(   dnsend_tim)

       call this%downsending_all(HALL,VALL,lquart)

                                                 call etim(   dnsend_tim)
deallocate(VM3D)
deallocate(VM2D)
deallocate(HM3D)
deallocate(HM2D)


!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_rad3   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          module              subroutine mg_filtering_lin1(this)
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 4:                                     !
!                                                                      !
!     - Multiple of 2D  line filter                                    !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 2d line filter only for all variables                          !
!                                                                      !
!***********************************************************************
implicit none
class (mg_intstate_type),target::this

integer(i_kind) L,i,j
integer(i_kind) icol,iout,jout
logical:: ff
include  "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include  "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------


!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend (Step 1)
!***
     
                                                 call btim( upsend_tim)
       call this%upsending_all(VALL,HALL,lquart)
                                                 call etim( upsend_tim)
!----------------------------------------------------------------------


!----------------------------------------------------------------------

!***
!*** Apply adjoint of Beta filter at all generations 
!***
                                                 call btim(    bfiltT_tim)



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

       do icol=3,1,-1
         call dibetat(km,i0-hx,i0,im,im+hx, j0-hy,j0,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VALL, ff, iout,jout)

         call this%bocoT_2d(VALL,km,im,jm,hx,hy)
       enddo

     do icol=3,1,-1
       if(l_hgen)  then
         call dibetat(km,i0-hx,i0,im,im+hx, j0-hy,j0,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HALL, ff, iout,jout)
       endif


        write(6,*)'thinkdeb33 4 ', km,im,jm,hx,hy 
         call this%bocoT_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)
     enddo


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfiltT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

      call this%weighting_all(VALL,HALL,lhelm)


                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations
!***


                                                 call btim( bfilt_tim)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!
       do icol=1,3
         call this%boco_2d(VALL,km,im,jm,hx,hy)
         call dibeta(km,i0-hx,0,im,im+hx, j0-hy,0,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VALL, ff, iout,jout)
       enddo

     do icol=1,3
         call this%boco_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)
       if(l_hgen)  then
         call dibeta(km,i0-hx,0,im,im+hx, j0-hy,0,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HALL, ff, iout,jout)
       endif
     enddo


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfilt_tim)

!***
!***  Downsend, interpolate and add, then zero high generations 
!***

                                                 call btim(   dnsend_tim)
       call this%downsending_all(HALL,VALL,lquart)

                                                 call etim(   dnsend_tim)


!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_lin1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
             module           subroutine mg_filtering_lin2(this)
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 5:                                     !
!                                                                      !
!     - Multiple of 2D  line filter                                    !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 2d radial filter + 1d vertical filter
!                                                                      !
!***********************************************************************
implicit none
class (mg_intstate_type),target::this

integer(i_kind) L,i,j
integer(i_kind) icol,iout,jout
logical:: ff

real(r_kind), allocatable, dimension(:,:,:):: VM2D
real(r_kind), allocatable, dimension(:,:,:):: HM2D
real(r_kind), allocatable, dimension(:,:,:,:):: VM3D
real(r_kind), allocatable, dimension(:,:,:,:):: HM3D
include  "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include  "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"

!----------------------------------------------------------------------

allocate(VM3D(km3,i0-hx:im+hx,j0-hy:jm+hy,lm))                 ; VM3D=0.
allocate(VM2D(km2,i0-hx:im+hx,j0-hy:jm+hy   ))                 ; VM2D=0.
allocate(HM3D(km3,i0-hx:im+hx,j0-hy:jm+hy,lm))                 ; HM3D=0.
allocate(HM2D(km2,i0-hx:im+hx,j0-hy:jm+hy   ))                 ; HM2D=0.


!-----------------------------------------------------------------------


!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend (Step 1)
!***
     
                                                 call btim( upsend_tim)
       call this%upsending_all(VALL,HALL,lquart)
                                                 call etim( upsend_tim)
!----------------------------------------------------------------------


!----------------------------------------------------------------------

!***
!*** Apply adjoint of Beta filter at all generations 
!***
                                                 call btim(    bfiltT_tim)



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Horizontal
!

       do icol=3,1,-1
         call dibetat(km,i0-hx,i0,im,im+hx, j0-hy,j0,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VALL, ff, iout,jout)
         call this%bocoT_2d(VALL,km,im,jm,hx,hy)
       enddo

     do icol=3,1,-1
       if(l_hgen)  then
         call dibetat(km,i0-hx,i0,im,im+hx, j0-hy,j0,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HALL, ff, iout,jout)
       endif
        write(6,*)'thinkdeb33 5 ', km,im,jm,hx,hy 
         call this%bocoT_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)
     enddo
!
! Vertical
!

       call this%stack_to_composite(VALL,VM2D,VM3D)
         call this%sup_vrbeta1T(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
       call this%composite_to_stack(VM2D,VM3D,VALL)

    if(l_hgen)  then
      call this%stack_to_composite(HALL,HM2D,HM3D)
        call this%sup_vrbeta1T(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,HM3D)
      call this%composite_to_stack(HM2D,HM3D,HALL)
    endif

        write(6,*)'thinkdeb33 6 ', km,im,jm,hx,hy 
        call this%bocoT_2d(VALL,km,im,jm,hx,hy)
        call this%bocoT_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                call etim(    bfiltT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

      call this%weighting_all(VALL,HALL,lhelm)


                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations
!***


                                                 call btim( bfilt_tim)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Horizontal
!
       do icol=1,3
         call this%boco_2d(VALL,km,im,jm,hx,hy)
         call dibeta(km,i0-hx,i0,im,im+hx, j0-hy,j0,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VALL, ff, iout,jout)
       enddo

     do icol=1,3
         call this%boco_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)
       if(l_hgen)  then
         call dibeta(km,i0-hx,i0,im,im+hx, j0-hy,j0,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HALL, ff, iout,jout)
       endif
     enddo
!
! Vertical
!

      call this%boco_2d(VALL,km,im,jm,hx,hy)
      call this%boco_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)


      call this%stack_to_composite(VALL,VM2D,VM3D)
        call this%sup_vrbeta1(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call this%composite_to_stack(VM2D,VM3D,VALL)

    if(l_hgen)  then
      call this%stack_to_composite(HALL,HM2D,HM3D)
        call this%sup_vrbeta1(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,HM3D)
      call this%composite_to_stack(HM2D,HM3D,HALL)
    endif


       call this%barrierMPI
!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


                                                 call etim(    bfilt_tim)

!***
!***  Downsend, interpolate and add, then zero high generations 
!***

                                                 call btim(   dnsend_tim)
       call this%downsending_all(HALL,VALL,lquart)

                                                 call etim(   dnsend_tim)


deallocate(VM3D)
deallocate(VM2D)
deallocate(HM3D)
deallocate(HM2D)


!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_lin2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          module              subroutine mg_filtering_lin3(this)
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 6:                                     !
!                                                                      !
!     - Multiple of 2D  line filter                                    !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 3d line filter                      
!                                                                      !
!***********************************************************************
!TEST
use, intrinsic :: ieee_arithmetic
!TEST
use jp_pkind2, only: fpi
implicit none
class (mg_intstate_type),target::this

integer(i_kind) k,i,j,L
integer(i_kind) icol,iout,jout,lout
logical:: ff

real(r_kind), allocatable, dimension(:,:,:):: VM2D
real(r_kind), allocatable, dimension(:,:,:):: HM2D
real(r_kind), allocatable, dimension(:,:,:,:):: VM3D
real(r_kind), allocatable, dimension(:,:,:,:):: HM3D

real(r_kind), allocatable, dimension(:,:,:,:):: W
real(r_kind), allocatable, dimension(:,:,:,:):: H

integer(fpi), allocatable, dimension(:,:,:):: JCOL
include  "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include  "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"


allocate(VM3D(km3,i0-hx:im+hx,j0-hy:jm+hy,lm))                 ; VM3D=0.
allocate(VM2D(km2,i0-hx:im+hx,j0-hy:jm+hy   ))                 ; VM2D=0.
allocate(HM3D(km3,i0-hx:im+hx,j0-hy:jm+hy,lm))                 ; HM3D=0.
allocate(HM2D(km2,i0-hx:im+hx,j0-hy:jm+hy   ))                 ; HM2D=0.

allocate(W(km3,i0-hx:im+hx,j0-hy:jm+hy,1-hz:lm+hz))            ; W=0.
allocate(H(km3,i0-hx:im+hx,j0-hy:jm+hy,1-hz:lm+hz))            ; H=0.

allocate(JCOL(0:im,0:jm,1:Lm))                                  ; JCOL=0

!-----------------------------------------------------------------------


!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend
!***
     
                                                 call btim( upsend_tim)
       call this%upsending_all(VALL,HALL,lquart)
                                                 call etim( upsend_tim)
!----------------------------------------------------------------------


!----------------------------------------------------------------------

!***
!*** Apply adjoint of Beta filter at all generations 
!***
                                                 call btim(    bfiltT_tim)

!
! From single stack to composite variables
!

       call this%stack_to_composite(VALL,VM2D,VM3D)
     if(l_hgen)  then
       call this%stack_to_composite(HALL,HM2D,HM3D)
     endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
!  Apply adjoint filter to 2D variables first
!

       do icol=3,1,-1
         call dibetat(km2,i0-hx,i0,im,im+hx, j0-hy,j0,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VM2D, ff, iout,jout)
        write(6,*)'thinkdeb33 11.0 ', km2,im,jm,hx,hy 
         call this%bocoT_2d(VM2D,km2,im,jm,hx,hy)
       enddo
    
     do icol=3,1,-1
       if(l_hgen)  then
         call dibetat(km2,i0-hx,i0,im,im+hx, j0-hy,j0,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HM2D, ff, iout,jout)
       endif
        write(6,*)'thinkdeb33 11 ', km2,im,jm,hx,hy 
         call this%bocoT_2d(HM2D,km2,im,jm,hx,hy,Fimax,Fjmax,2,gm)
     enddo

!
! Create and apply adjoint filter to extended 3D variables
!

         W(:,:,:,1:lm)=VM3D(:,:,:,1:lm)

       do icol=7,1,-1
         do L=1,hz
           W(:,:,:,1-L )=0.
           W(:,:,:,LM+L)=0.
         end do
         call dibetat(km3,i0-hx,i0,im,im+hx, j0-hy,j0,jm,jm+hy, 1-hz,1,lm,lm+hz,icol, nfil  &
                     ,qcols,dixs3,diys3,dizs3,JCOL,vpasp3, W, ff, iout,jout,lout)
         call this%bocoT_3d(W,km3,im,jm,Lm,hx,hy,hz,Fimax,Fjmax)
       enddo

     if(l_hgen)  then
          H(:,:,:,1:lm)=HM3D(:,:,:,1:lm)
     endif

     do icol=7,1,-1
       if(l_hgen)  then
         do L=1,hz
           H(:,:,:,1-L )=0.
           H(:,:,:,LM+L)=0.
         end do

         call dibetat(km3,i0-hx,i0,im,im+hx, j0-hy,j0,jm,jm+hy, 1-hz,1,lm,lm+hz,icol, nfil  &
                     ,qcols,dixs3,diys3,dizs3,JCOL,vpasp3, H, ff, iout,jout,lout)
       endif
         call this%bocoT_3d(H,km3,im,jm,Lm,hx,hy,hz,Fimax,Fjmax,2,gm)
     enddo


!
! Go back from extended 3D variables and combine them with 2D variables in one stacked variable
!

       VM3D(:,:,:,1:lm)= W(:,:,:,1:lm)
       call this%composite_to_stack(VM2D,VM3D,VALL)

     if(l_hgen)  then
       HM3D(:,:,:,1:lm)=H(:,:,:,1:lm)
       call this%composite_to_stack(HM2D,HM3D,HALL)
     endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfiltT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

      call this%weighting_all(VALL,HALL,lhelm)


                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations
!***
!
! From single stacked to composite variables
!
                                                 call btim( bfilt_tim)

       call this%stack_to_composite(VALL,VM2D,VM3D)
     if(l_hgen)  then
       call this%stack_to_composite(HALL,HM2D,HM3D)
     endif



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
!  Apply filter to 2D variables first
!
       do icol=1,3
         call this%boco_2d(VM2D,km2,im,jm,hx,hy)
         call dibeta(km2,i0-hx,i0,im,im+hx, j0-hy,j0,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VM2D, ff, iout,jout)
       enddo

     do icol=1,3
         call this%boco_2d(HM2D,km2,im,jm,hx,hy,Fimax,Fjmax,2,gm)
       if(l_hgen)  then
         call dibeta(km2,i0-hx,i0,im,im+hx, j0-hy,j0,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HM2D, ff, iout,jout)
       endif
     enddo

!
! Create and apply filter to extended 3D variables
!

           W(:,:,:,1:lm)=VM3D(:,:,:,1:lm)
        do L=1,hz
          do j=j0-hy,jm+hy
          do i=i0-hx,im+hx
              W(:,i,j,1-L )=VM3D(:,i,j, 1+L)
              W(:,i,j,LM+L)=VM3D(:,i,j,LM-L)
          end do
          end do
        end do

       do icol=1,7
         call this%boco_3d(W,km3,im,jm,lm,hx,hy,hz,Fimax,Fjmax)
         call dibeta(km3,i0-hx,i0,im,im+hx, j0-hy,j0,jm,jm+hy, 1-hz,1,lm,lm+hz,icol, nfil  &
                    ,qcols,dixs3,diys3,dizs3,JCOL,vpasp3, W, ff, iout,jout,lout)
        enddo 

     if(l_hgen)  then
           H(:,:,:,1:lm)=HM3D(:,:,:,1:lm)
        do L=1,hz
          do j=j0-hy,jm+hy
          do i=i0-hx,im+hx
              H(:,i,j,1-L )=HM3D(:,i,j, 1+L)
              H(:,i,j,LM+L)=HM3D(:,i,j,LM-L)
          end do
          end do
        end do
     endif
       do icol=1,7
         call this%boco_3d(H,km3,im,jm,lm,hx,hy,hz,Fimax,Fjmax,2,gm)
         if(l_hgen)  then
         call dibeta(km3,i0-hx,i0,im,im+hx, j0-hy,j0,jm,jm+hy, 1-hz,1,lm,lm+hz,icol, nfil  &
                    ,qcols,dixs3,diys3,dizs3,JCOL,vpasp3, H, ff, iout,jout,lout)
         endif
       enddo

!
! Go back from extended 3D variables and combine them with 2D variables in one stacked variable
!

       VM3D(:,:,:,1:lm)= W(:,:,:,1:lm)
       call this%composite_to_stack(VM2D,VM3D,VALL)

     if(l_hgen)  then
       HM3D(:,:,:,1:lm)=H(:,:,:,1:lm)
       call this%composite_to_stack(HM2D,HM3D,HALL)
     endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfilt_tim)

!***
!***  Downsend, interpolate and add, then zero high generations 
!***

                                                 call btim(   dnsend_tim)
       call this%downsending_all(HALL,VALL,lquart)

                                                 call etim(   dnsend_tim)


!-----------------------------------------------------------------------

deallocate(VM3D)
deallocate(VM2D)
deallocate(HM3D)
deallocate(HM2D)

deallocate(W)
deallocate(H)

deallocate(JCOL)

!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_lin3

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          module              subroutine mg_filtering_fast(this)
!***********************************************************************
!                                                                      !
! Fast multigrid filtering procedure:                                  !
!                                                                      !
!     - Multiple of 2D and 3D variables                                !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 1d+1d horizontal filter + 1d vertical filter                   !
!                                                                      !
!***********************************************************************
implicit none
class (mg_intstate_type),target::this

real(r_kind), allocatable, dimension(:,:,:):: VM2D
real(r_kind), allocatable, dimension(:,:,:):: HM2D
real(r_kind), allocatable, dimension(:,:,:,:):: VM3D
real(r_kind), allocatable, dimension(:,:,:,:):: HM3D

integer(i_kind) L,i,j
include  "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include  "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------

allocate(VM3D(km3,i0-hx:im+hx,j0-hy:jm+hy,lm))                  ; VM3D=0.
allocate(VM2D(km2,i0-hx:im+hx,j0-hy:jm+hy   ))                  ; VM2D=0.
allocate(HM3D(km3,i0-hx:im+hx,j0-hy:jm+hy,lm))                  ; HM3D=0.
allocate(HM2D(km2,i0-hx:im+hx,j0-hy:jm+hy   ))                  ; HM2D=0.



!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend 
!***
     
                                                 call btim( upsend_tim)
       call this%upsending_all(VALL,HALL,lquart)
                                                 call etim( upsend_tim)
!----------------------------------------------------------------------


!----------------------------------------------------------------------

!***
!*** Apply adjoint of Beta filter at all generations 
!***
                                                 call btim(    bfiltT_tim)



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Horizontally
!

    do j=0,jm
      call this%rbetaT(km,hx,1,im,paspx,ssx,VALL(:,:,j))
    enddo
       call this%bocoTx(VALL,km,im,jm,hx,hy)

    do i=0,im
      call this%rbetaT(km,hy,1,jm,paspy,ssy,VALL(:,i,:))
    enddo
      call this%bocoTy(VALL,km,im,jm,hx,hy)

      call this%stack_to_composite(VALL,VM2D,VM3D)

  if(l_hgen)  then
    do j=0,jm
      call this%rbetaT(km,hx,i0,im,paspx,ssx,HALL(:,:,j))
    enddo
  endif
      call this%bocoTx(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)

  if(l_hgen)  then
    do i=0,im
      call this%rbetaT(km,hy,j0,jm,paspy,ssy,HALL(:,i,:))
    enddo
  endif
      call this%bocoTy(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)

!
! Vertically 
!
      call this%stack_to_composite(HALL,HM2D,HM3D)
      call this%sup_vrbeta1T(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call this%composite_to_stack(VM2D,VM3D,VALL)
  if(l_hgen)  then
      call this%sup_vrbeta1T(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,HM3D)
      call this%composite_to_stack(HM2D,HM3D,HALL)
   endif


       call this%barrierMPI
!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


                                                 call etim(    bfiltT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

      call this%weighting_all(VALL,HALL,lhelm)


                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations (Step 7)
!***
                                                 call btim( bfilt_tim)


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!
! Horizonatally

     call this%bocox(VALL,km,im,jm,hx,hy)
    do j=0,jm
      call this%rbeta(km,hx,i0,im,paspx,ssx,VALL(:,:,j))
    enddo

      call this%bocoy(VALL,km,im,jm,hx,hy)
    do i=0,im
      call this%rbeta(km,hy,j0,jm,paspy,ssy,VALL(:,i,:))
    enddo

      call this%stack_to_composite(VALL,VM2D,VM3D)

      call this%bocox(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)
  if(l_hgen)  then
    do j=0,jm
      call this%rbeta(km,hx,i0,im,paspx,ssx,HALL(:,:,j))
    enddo
  endif
      call this%bocoy(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)
  if(l_hgen)  then
    do i=0,im
      call this%rbeta(km,hy,j0,jm,paspy,ssy,HALL(:,i,:))
    enddo
  endif
  if(l_hgen)  then
    call this%stack_to_composite(HALL,HM2D,HM3D)
  endif

!
! Vertically
!

      call this%sup_vrbeta1(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call this%composite_to_stack(VM2D,VM3D,VALL)
  if(l_hgen)  then
      call this%sup_vrbeta1(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,HM3D)
      call this%composite_to_stack(HM2D,HM3D,HALL)
   endif

       call this%barrierMPI
!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfilt_tim)

!***
!***  Downsend, interpolate and add (Step 4)
!***  Then zero high generations (Step 5)
!***

                                                 call btim(   dnsend_tim)
       call this%downsending_all(HALL,VALL,lquart)

                                                 call etim(   dnsend_tim)

deallocate(VM3D) 
deallocate(VM2D)
deallocate(HM3D)
deallocate(HM2D)

!-----------------------------------------------------------------------
                        endsubroutine mg_filtering_fast

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
           module             subroutine sup_vrbeta1                        &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta1                                          *
!                                                                     *
!**********************************************************************
(this,kmax,hx,hy,hz,im,jm,lm, pasp,ss, V)
!----------------------------------------------------------------------
implicit none
        class(mg_intstate_type),target::this

integer(i_kind),intent(in):: kmax,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:kmax,this%i0-hx:im+hx,this%j0-hy:jm+hy,1:lm),intent(inout):: V
real(r_kind),dimension(1,1,1:lm), intent(in):: pasp
real(r_kind),dimension(1:lm), intent(in):: ss

real(r_kind),dimension(1:kmax,1-hz:lm+hz):: W

integer(i_kind):: i,j,L

!----------------------------------------------------------------------

        do j=this%j0,jm
        do i=this%i0,im
          do L=1,Lm
            W(:,L)=V(:,i,j,L)
          end do
          do L=1,hz
            W(:,1-L)=W(:,1+L)
            W(:,LM+L)=W(:,LM-L)
          end do
             call this%rbeta(kmax,hz,1,lm,  pasp,ss,W)
          do l=1,Lm
            V(:,i,j,L)=W(:,L)
          end do
        end do
        end do

  
!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
             module           subroutine sup_vrbeta1T                        &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta1T                                         *
!                                                                     *
!**********************************************************************
(this,kmax,hx,hy,hz,im,jm,lm,  pasp,ss, V)
!----------------------------------------------------------------------
implicit none
        class(mg_intstate_type),target::this

integer(i_kind),intent(in):: kmax,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:kmax,this%i0-hx:im+hx,this%j0-hy:jm+hy,1:lm),intent(inout):: V
real(r_kind),dimension(1,1,1:lm), intent(in):: pasp
real(r_kind),dimension(1:lm), intent(in):: ss

real(r_kind),dimension(1:kmax,1-hz:lm+hz):: W

integer(i_kind):: i,j,L

!----------------------------------------------------------------------

        do j=this%j0,jm
        do i=this%i0,im
          do L=1,Lm
            W(:,L)=V(:,i,j,L)
          end do
          do L=1,hz
            W(:,1-L )=W(:,1+L )
            W(:,LM+L)=W(:,LM-L)
          end do
             call this%rbetaT(kmax,hz,1,lm, pasp,ss,W)
!
! Apply adjoint at the edges of domain
!
          do L=1,hz
            W(:,1+L)=W(:,1+L)+W(:,1-L)
            W(:,LM-L)=W(:,LM-L)+W(:,LM+L)
          enddo
          do l=1,Lm
            V(:,i,j,L)=W(:,L)
          end do
        end do
        end do

!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta1T

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                  module      subroutine sup_vrbeta3                        &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta3                                           *
!                                                                     *
!**********************************************************************
(this,kmax,hx,hy,hz,im,jm,lm, pasp,ss, V)
!----------------------------------------------------------------------
implicit none
        class(mg_intstate_type),target::this

integer(i_kind),intent(in):: kmax,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:kmax,this%i0-hx:im+hx,this%j0-hy:jm+hy,1:lm),intent(inout):: V
real(r_kind),dimension(3,3,this%i0:im,this%j0:jm,1:lm), intent(in):: pasp
real(r_kind),dimension(this%i0:im,this%j0:jm,1:lm), intent(in):: ss

real(r_kind),dimension(1:kmax,this%i0-hx:im+hx,this%j0-hy:jm+hy,1-hz:lm+hz):: W

integer(i_kind):: i,j,L

!----------------------------------------------------------------------

          do L=1,Lm
          do j=this%j0-hy,jm+hy
          do i=this%i0-hx,im+hx
            W(:,i,j,L)=V(:,i,j,L)
          end do
          end do
          end do

        do L=1,hz
          do j=this%j0-hy,jm+hy
          do i=this%i0-hx,im+hx
              W(:,i,j,1-L )=W(:,i,j,1+L )
              W(:,i,j,LM+L)=W(:,i,j,LM-L)
          end do
          end do
        end do
    
    
           call this%rbeta(kmax,hx,this%i0,im, hy,this%j0,jm, hz,1,lm, pasp,ss,W)

  
          do l=1,Lm
          do j=this%j0,jm
          do i=this%i0,im
            V(:,i,j,L)=W(:,i,j,L)
          end do
          end do
          end do

!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta3
 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
             module           subroutine sup_vrbeta3T                       &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta3                                           *
!                                                                     *
!**********************************************************************
(this,kmax,hx,hy,hz,im,jm,lm, pasp,ss,V)
!----------------------------------------------------------------------
implicit none
        class(mg_intstate_type), target::this

integer(i_kind),intent(in):: kmax,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:kmax,this%i0-hx:im+hx,this%j0-hy:jm+hy,1:lm),intent(inout):: V
real(r_kind),dimension(3,3,this%i0:im,this%j0:jm,1:lm), intent(in):: pasp
real(r_kind),dimension(this%i0:im,this%j0:jm,1:lm), intent(in):: ss

real(r_kind),dimension(1:kmax,this%i0-hx:im+hx,this%j0-hy:jm+hy,1-hz:lm+hz):: W

integer(i_kind):: i,j,l

!----------------------------------------------------------------------

          do L=1,Lm
          do j=this%j0-hy,jm+hy
          do i=this%i0-hx,im+hx
            W(:,i,j,L)=V(:,i,j,L)
          end do
          end do
          end do

        do L=1,hz
          do j=this%j0-hy,jm+hy
          do i=this%i0-hx,im+hx
              W(:,i,j,1-L )=W(:,i,j, 1+L)
              W(:,i,j,LM+L)=W(:,i,j,LM-L)
          end do
          end do
        end do
    
    
           call this%rbetaT(kmax,hx,this%i0,im, hy,this%j0,jm, hz,1,lm, pasp,ss,W)

!
! Apply adjoint at the edges of domain
!
        do L=1,hz
          do j=this%j0-hy,jm+hy
          do i=this%i0-hx,im+hx
              W(:,i,j,1+L )=W(:,i,j, 1+L)+W(:,i,j, 1-L)
              W(:,i,j,LM-L)=W(:,i,j,LM-L)+W(:,i,j,LM+L)
          end do
          end do
         end do
  
          do l=1,lm
          do j=this%j0,jm
          do i=this%i0,im
            V(:,i,j,l)=W(:,i,j,l)
          end do
          end do
          end do

!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta3T

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        end submodule mg_filtering
