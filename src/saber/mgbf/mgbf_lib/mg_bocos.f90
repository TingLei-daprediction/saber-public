submodule(mg_intstate) mg_bocos
!$$$  submodule documentation block
!                .      .    .                                       .
! module:   mg_bocos
!   prgmmr: rancic           org: NCEP/EMC            date: 2022
!
! abstract:  Provide communication between subdomains and supply halos
!            on filter grid (offset version)
!
! module history log:
!   2023-04-19  lei     - object-oriented coding
!   2024-01-11  rancic  - optimization for ensemble localization
!   2024-02-20  yokota  - refactoring to apply for GSI
!
! Subroutines Included:
!   boco_2d_g1 -
!   boco_2d_gh -
!   bocoT_2d_g1 -
!   bocoT_2d_gh -
!   boco_3d_g1 -
!   boco_3d_gh -
!   bocoT_3d_g1 -
!   bocoT_3d_gh -
!   upsend_all_g1 -
!   upsend_all_gh -
!   downsend_all_gh -
!   downsend_all_g2 -
!   bocox_2d_g1 -
!   bocox_2d_gh -
!   bocoy_2d_g1 -
!   bocoy_2d_gh -
!   bocoTx_2d_g1 -
!   bocoTx_2d_gh -
!   bocoTy_2d_g1 -
!   bocoTy_2d_gh -
!   boco_2d_loc -
!   bocoT_2d_loc -
!   upsend_loc_g12 -
!   upsend_loc_g23 -
!   upsend_loc_g34 -
!   downsend_loc_g43 -
!   downsend_loc_g32 -
!   downsend_loc_g21 -
!
! Functions Included:
!
! remarks:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use mgbf_kinds, only: r_kind,i_kind

implicit none

interface boco_2d
  module procedure boco_2d_g1
  module procedure boco_2d_gh
endinterface

interface bocoT_2d
  module procedure bocoT_2d_g1
  module procedure bocoT_2d_gh
endinterface

interface boco_3d
  module procedure boco_3d_g1 
  module procedure boco_3d_gh 
endinterface

interface bocoT_3d
  module procedure bocoT_3d_g1 
  module procedure bocoT_3d_gh 
endinterface

interface upsend_all
  module procedure upsend_all_g1
  module procedure upsend_all_gh
endinterface

interface downsend_all     
  module procedure downsend_all_gh
  module procedure downsend_all_g2
endinterface

interface bocox
  module procedure bocox_2d_g1
  module procedure bocox_2d_gh
endinterface

interface bocoy
  module procedure bocoy_2d_g1
  module procedure bocoy_2d_gh
endinterface

interface bocoTx
  module procedure bocoTx_2d_g1
  module procedure bocoTx_2d_gh
endinterface

interface bocoTy
  module procedure bocoTy_2d_g1
  module procedure bocoTy_2d_gh
endinterface
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine boco_2d_g1 &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies (nbx,nby) lines of halos in (x,y) directions assuming       !
! mirror boundary conditions. Version for generation 1                 !
!                                                                      !
!                        - offset version -                            !
!                                                                      !
!**********************************************************************!
(this,W,km_in,im_in,jm_in,nbx,nby)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby
real(r_kind),dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
!-----------------------------------------------------------------------
integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth
integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,i,j
integer(i_kind) ndatay,ndatax
integer(i_kind) g_ind
integer(i_kind), parameter :: DIR_S=1, DIR_N=2, DIR_E=3, DIR_W=4
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain
!
          g_ind = 1

          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)

          imax = im_in
          jmax = jm_in

!-----------------------------------------------------------------------
      ndatay = km_in*imax*nby
      ndatax = km_in*(jmax+2*nby)*nbx

      sHandle(:) = MPI_REQUEST_NULL
      rHandle(:) = MPI_REQUEST_NULL

!
!  Post receives toward SOUTH and NORTH
!
      if( itarg_s >= 0 ) then
        call MPI_IRECV( W(1,1,1-nby), ndatay, dtype, itarg_s, itarg_s, &
                        mpi_comm_comp, rHandle(DIR_S), ierr )
      end if

      if( itarg_n >= 0 ) then
        call MPI_IRECV( W(1,1,jmax+1), ndatay, dtype, itarg_n, itarg_n, &
                        mpi_comm_comp, rHandle(DIR_N), ierr )
      end if

!
!  Post receives toward WEST and EAST
!
      if( itarg_w >= 0 ) then
        call MPI_IRECV( W(1,1-nbx,1-nby), ndatax, dtype, itarg_w, itarg_w, &
                        mpi_comm_comp, rHandle(DIR_W), ierr )
      end if

      if( itarg_e >= 0 ) then
        call MPI_IRECV( W(1,imax+1,1-nby), ndatax, dtype, itarg_e, itarg_e, &
                        mpi_comm_comp, rHandle(DIR_E), ierr )
      end if

!
!  SEND boundaries toward SOUTH and NORTH
!
      if( itarg_s >= 0 ) then
        call MPI_ISEND( W(1,1,1), ndatay, dtype, itarg_s, mype, &
                        mpi_comm_comp, sHandle(DIR_S), ierr )
      end if

      if( itarg_n >= 0 ) then
        call MPI_ISEND( W(1,1,jmax-nby+1), ndatay, dtype, itarg_n, mype, &
                        mpi_comm_comp, sHandle(DIR_N), ierr )
      end if

!
!  SEND extended boundaries toward WEST and EAST
!
      if( itarg_w >= 0 ) then
        call MPI_ISEND( W(1,1,1-nby), ndatax, dtype, itarg_w, mype, &
                        mpi_comm_comp, sHandle(DIR_W), ierr )
      end if

      if( itarg_e >= 0 ) then
        call MPI_ISEND( W(1,imax-nbx+1,1-nby), ndatax, dtype, itarg_e, mype, &
                        mpi_comm_comp, sHandle(DIR_E), ierr )
      end if

!
!  Complete non-blocking receives and sends
!
      if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
      end if

      if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
      end if

!
!  Fill physical boundaries when required
!
   if( lsouth ) then
     do j=1,nby
       do i=1,imax
         W(:,i,-nby+j) = W(:,i,nby+1-j)
       end do
     end do
   end if

   if( lnorth ) then
     do j=1,nby
       do i=1,imax
         W(:,i,jmax+j) = W(:,i,jmax+1-j)
       end do
     end do
   end if

   if( lwest ) then
     do j=1-nby,jmax+nby
       do i=1,nbx
         W(:,-nbx+i,j) = W(:,nbx+1-i,j)
       end do
     end do
   end if

   if( least ) then
     do j=1-nby,jmax+nby
       do i=1,nbx
         W(:,imax+i,j) = W(:,imax+1-i,j)
       end do
     end do
   end if

!-----------------------------------------------------------------------
endsubroutine boco_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine boco_2d_gh &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies (nbx,nby) lines of halos in (x,y) directions assuming       !
! mirror boundary conditions. Version for high generations             !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!**********************************************************************!
(this,W,km_in,im_in,jm_in,nbx,nby,Fimax_in,Fjmax_in,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby,mygen_min,mygen_max
real(r_kind),dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
!-----------------------------------------------------------------------
integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth
integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,i,j
integer(i_kind) ndatay,ndatax
integer(i_kind) g_ind,g
logical l_sidesend
integer(i_kind), parameter :: DIR_S=1, DIR_N=2, DIR_E=3, DIR_W=4
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!
 
       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2 
         g = my_hgen
         l_sidesend=.true.
       else
         l_sidesend=.false.
       endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 

          if(least) then
            imax = Fimax_in(g)
          else 
            imax = im_in       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax_in(g)
          else  
            jmax = jm_in
          endif

!-----------------------------------------------------------------------
      ndatay = km_in*imax*nby
      ndatax = km_in*(jmax+2*nby)*nbx

      sHandle(:) = MPI_REQUEST_NULL
      rHandle(:) = MPI_REQUEST_NULL

!
!  Post receives toward SOUTH and NORTH
!
      if( itarg_s >= 0 ) then
        call MPI_IRECV( W(1,1,1-nby), ndatay, dtype, itarg_s, itarg_s, &
                        mpi_comm_work, rHandle(DIR_S), ierr )
      end if

      if( itarg_n >= 0 ) then
        call MPI_IRECV( W(1,1,jmax+1), ndatay, dtype, itarg_n, itarg_n, &
                        mpi_comm_work, rHandle(DIR_N), ierr )
      end if

!
!  Post receives toward WEST and EAST
!
      if( itarg_w >= 0 ) then
        call MPI_IRECV( W(1,1-nbx,1-nby), ndatax, dtype, itarg_w, itarg_w, &
                        mpi_comm_work, rHandle(DIR_W), ierr )
      end if

      if( itarg_e >= 0 ) then
        call MPI_IRECV( W(1,imax+1,1-nby), ndatax, dtype, itarg_e, itarg_e, &
                        mpi_comm_work, rHandle(DIR_E), ierr )
      end if

!
!  SEND boundaries to SOUTH and NORTH
!
      if( itarg_s >= 0 ) then
        call MPI_ISEND( W(1,1,1), ndatay, dtype, itarg_s, mype, &
                        mpi_comm_work, sHandle(DIR_S), ierr )
      end if

      if( itarg_n >= 0 ) then
        call MPI_ISEND( W(1,1,jmax-nby+1), ndatay, dtype, itarg_n, mype, &
                        mpi_comm_work, sHandle(DIR_N), ierr )
      end if

!
!  SEND extended boundaries to WEST and EAST
!
      if( itarg_w >= 0 ) then
        call MPI_ISEND( W(1,1,1-nby), ndatax, dtype, itarg_w, mype, &
                        mpi_comm_work, sHandle(DIR_W), ierr )
      end if

      if( itarg_e >= 0 ) then
        call MPI_ISEND( W(1,imax-nbx+1,1-nby), ndatax, dtype, itarg_e, mype, &
                        mpi_comm_work, sHandle(DIR_E), ierr )
      end if

!
!  Complete non-blocking receives and sends
!
      if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
      end if

      if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
      end if

!
!  Fill physical boundaries when required
!
   if( lsouth ) then
     do j=1,nby
       do i=1,imax
         W(:,i,-nby+j) = W(:,i,nby+1-j)
       end do
     end do
   end if

   if( lnorth ) then
     do j=1,nby
       do i=1,imax
         W(:,i,jmax+j) = W(:,i,jmax+1-j)
       end do
     end do
   end if

   if( lwest ) then
     do j=1-nby,jmax+nby
       do i=1,nbx
         W(:,-nbx+i,j) = W(:,nbx+1-i,j)
       end do
     end do
   end if

   if( least ) then
     do j=1-nby,jmax+nby
       do i=1,nbx
         W(:,imax+i,j) = W(:,imax+1-i,j)
       end do
     end do
   end if

!-----------------------------------------------------------------------

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
endsubroutine boco_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine bocoT_2d_g1 &
!***********************************************************************
!                                                                      !
! Adjoint of side sending subroutine:                                  !
! Supplies (nbx,nby) lines of halos in (x,y) directions, including     !
! values at the edges of the subdomains and assuming mirror boundary   !
! conditions just for generation 1                                     !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(this,W,km_in,im_in,jm_in,nbx,nby)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby
real(r_kind),dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
!-----------------------------------------------------------------------
integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical lwest,least,lsouth,lnorth                                       
integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,i,j
integer(i_kind) ndatay,ndatax
integer(i_kind), parameter :: DIR_S=1, DIR_N=2, DIR_E=3, DIR_W=4
real(r_kind) :: rBuf_W(km_in,nbx,1-nby:jm_in+nby)
real(r_kind) :: rBuf_E(km_in,nbx,1-nby:jm_in+nby)
real(r_kind) :: rBuf_S(km_in,im_in,nby)
real(r_kind) :: rBuf_N(km_in,im_in,nby)
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!
! Limit comminications to selected number of generations
!

         itarg_n = Fitarg_n(1)
         itarg_s = Fitarg_s(1)
         itarg_w = Fitarg_w(1)
         itarg_e = Fitarg_e(1)

         lwest   = Flwest(1)
         least   = Fleast(1)
         lsouth  = Flsouth(1)
         lnorth  = Flnorth(1)

         imax = im_in    
         jmax = jm_in

!----------------------------------------------------------------------
      ndatax =km_in*(jmax+2*nby)*nbx
      ndatay =km_in*imax*nby

      sHandle(:) = MPI_REQUEST_NULL
      rHandle(:) = MPI_REQUEST_NULL

!
! RECEIVE extended halos from EAST and WEST
!
      if( itarg_w >= 0 ) then
        call MPI_IRECV( rBuf_W, ndatax, dtype, itarg_w, itarg_w, &
                        mpi_comm_comp, rHandle(DIR_W), ierr )
      end if

      if( itarg_e >= 0 ) then
        call MPI_IRECV( rBuf_E, ndatax, dtype, itarg_e, itarg_e, &
                        mpi_comm_comp, rHandle(DIR_E), ierr )
      end if

!
! RECEIVE boundaries from SOUTH and NORTH
!
      if( itarg_s >= 0 ) then
        call MPI_IRECV( rBuf_S, ndatay, dtype, itarg_s, itarg_s, &
                        mpi_comm_comp, rHandle(DIR_S), ierr )
      end if

      if( itarg_n >= 0 ) then
        call MPI_IRECV( rBuf_N, ndatay, dtype, itarg_n, itarg_n, &
                        mpi_comm_comp, rHandle(DIR_N), ierr )
      end if

!
! SEND extended halos toward WEST and EAST
!
      if( itarg_w >= 0) then
        call MPI_ISEND( W(1,1-nbx,1-nby), ndatax, dtype, itarg_w, mype, &
                        mpi_comm_comp, sHandle(DIR_W), ierr )
      end if

      if( itarg_e >= 0 ) then
        call MPI_ISEND( W(1,imax+1,1-nby), ndatax, dtype, itarg_e, mype, &
                        mpi_comm_comp, sHandle(DIR_E), ierr )
      end if

!
! SEND boundaries SOUTH and NORTH
!
      if( itarg_s >= 0 ) then
        call MPI_ISEND( W(1,1,1-nby), ndatay, dtype, itarg_s, mype, &
                        mpi_comm_comp, sHandle(DIR_S), ierr )
      end if

      if( itarg_n >= 0 ) then
        call MPI_ISEND( W(1,1,jmax+1), ndatay, dtype, itarg_n, mype, &
                        mpi_comm_comp, sHandle(DIR_N), ierr )
      end if

!
!  Complete non-blocking operations
!
      if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
      end if

      if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
      end if

!
! Assign received halos from WEST and EAST to interrior of domains
!
   if(lwest) then
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,i,j)= W(:,i,j)+W(:,1-i,j)
     end do
     end do
   else if( itarg_w >= 0 ) then
     do j=1-nby,jmax+nby
     do i=1,nbx
      W(:,i,j)= W(:,i,j)+rBuf_W(:,i,j)
     end do
     end do
   end if

   if(least) then
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+W(:,imax+1+nbx-i,j)
     end do
     end do
   else if( itarg_e >= 0 ) then
     do j=1-nby,jmax+nby
     do i=1,nbx  
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+rBuf_E(:,i,j)
     end do
     end do
   end if

!
! ASSIGN received values from SOUTH and NORTH
!
   if(lsouth) then
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+W(:,i,1-j)
     end do
     end do
   else if( itarg_s >= 0 ) then
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+rBuf_S(:,i,j)
     end do
     end do
   end if

   if(lnorth) then
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+W(:,i,jmax+1+nby-j)
     enddo
     enddo
   else if( itarg_n >= 0 ) then
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+rBuf_N(:,i,j)
     enddo
     enddo
   end if

!-----------------------------------------------------------------------
endsubroutine bocoT_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine bocoT_2d_gh &
!***********************************************************************
!                                                                      !
!  Supply n-lines inside of domains, including edges, with halos from  !
!  the surrounding domains.  Assume mirror boundary conditions at the  !
!  boundaries of the domain. For high multigrid generations.           !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(this,W,km_in,im_in,jm_in,nbx,nby,Fimax_in,Fjmax_in,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby,mygen_min,mygen_max
real(r_kind), dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:)::                           &
                                   rBuf_N,rBuf_E,rBuf_S,rBuf_W   
integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical lwest,least,lsouth,lnorth                                       
integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,i,j
integer(i_kind) ndatay,ndatax
integer(i_kind) g_ind,g
logical l_sidesend
integer(i_kind), parameter :: DIR_S=1, DIR_N=2, DIR_E=3, DIR_W=4
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!
! Limit comminications to selected number of generations
!


       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2
         g = my_hgen
         l_sidesend=.true.
       else 
         l_sidesend=.false.
       endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 

          if(least) then
            imax = Fimax_in(g)
          else 
            imax = im_in       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax_in(g)
          else  
            jmax = jm_in
          endif


!-----------------------------------------------------------------------
      ndatay = km_in*imax*nby
      ndatax = km_in*(jmax+2*nby)*nbx

      sHandle(:) = MPI_REQUEST_NULL
      rHandle(:) = MPI_REQUEST_NULL

!
!  Post receives from neighbours
!
      if( itarg_w >= 0 ) then
        allocate( rBuf_W(km_in,nbx,1-nby:jmax+nby), stat = ierr )
        call MPI_IRECV( rBuf_W, ndatax, dtype, itarg_w, itarg_w, &
                        mpi_comm_work, rHandle(DIR_W), ierr )
      end if

      if( itarg_e >= 0 ) then
        allocate( rBuf_E(km_in,nbx,1-nby:jmax+nby), stat = ierr )
        call MPI_IRECV( rBuf_E, ndatax, dtype, itarg_e, itarg_e, &
                        mpi_comm_work, rHandle(DIR_E), ierr )
      end if

      if( itarg_s >= 0 ) then
        allocate( rBuf_S(km_in,imax,nby), stat = ierr )
        call MPI_IRECV( rBuf_S, ndatay, dtype, itarg_s, itarg_s, &
                        mpi_comm_work, rHandle(DIR_S), ierr )
      end if

      if( itarg_n >= 0 ) then
        allocate( rBuf_N(km_in,imax,nby), stat = ierr )
        call MPI_IRECV( rBuf_N, ndatay, dtype, itarg_n, itarg_n, &
                        mpi_comm_work, rHandle(DIR_N), ierr )
      end if

!
!  SEND halos to neighbours
!
      if( itarg_w >= 0) then
        call MPI_ISEND( W(1,1-nbx,1-nby), ndatax, dtype, itarg_w, mype, &
                        mpi_comm_work, sHandle(DIR_W), ierr )
      end if

      if( itarg_e >= 0 ) then
        call MPI_ISEND( W(1,imax+1,1-nby), ndatax, dtype, itarg_e, mype, &
                        mpi_comm_work, sHandle(DIR_E), ierr )
      end if

      if( itarg_s >= 0 ) then
        call MPI_ISEND( W(1,1,1-nby), ndatay, dtype, itarg_s, mype, &
                        mpi_comm_work, sHandle(DIR_S), ierr )
      end if

      if( itarg_n >= 0 ) then
        call MPI_ISEND( W(1,1,jmax+1), ndatay, dtype, itarg_n, mype, &
                        mpi_comm_work, sHandle(DIR_N), ierr )
      end if

!
!  Complete non-blocking operations
!
      if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
      end if

      if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
      end if

!
! Assign received halos from WEST and EAST
!
   if(lwest) then

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,i,j)= W(:,i,j)+W(:,1-i,j)
     end do
     end do

   else if( itarg_w >= 0 ) then 

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,i,j)= W(:,i,j)+rBuf_W(:,i,j)
     enddo
     enddo

   end if

   if(least) then

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+W(:,imax+1+nbx-i,j)
     end do
     end do

   else if( itarg_e >= 0 ) then 

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+rBuf_E(:,i,j)
     enddo
     enddo

   end if

!
! Assign received values from SOUTH and NORTH
!
   if(lsouth) then
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+W(:,i,1-j)
     end do
     end do
   else if( itarg_s >= 0 ) then
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+rBuf_S(:,i,j)
     end do
     end do
   end if

   if(lnorth) then
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+W(:,i,jmax+1+nby-j)
     enddo
     enddo
   else if( itarg_n >= 0 ) then
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+rBuf_N(:,i,j)
     enddo
     enddo
   end if

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!
      if( allocated(rBuf_W) ) deallocate( rBuf_W, stat = ierr)
      if( allocated(rBuf_E) ) deallocate( rBuf_E, stat = ierr)
      if( allocated(rBuf_S) ) deallocate( rBuf_S, stat = ierr)
      if( allocated(rBuf_N) ) deallocate( rBuf_N, stat = ierr)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
endsubroutine bocoT_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine boco_3d_g1 &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies (nbx,nby) lines of halos in (x,y) directions assuming       !
! mirror boundary conditions. Version for generation 1                 !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!**********************************************************************!
(this,W,km_in,im_in,jm_in,Lm_in,nbx,nby,nbz,Fimax_in,Fjmax_in)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,Lm_in,nbx,nby,nbz
real(r_kind),dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby,1-nbz:Lm_in+nbz)    &
                      ,intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
!-----------------------------------------------------------------------
integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth
integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,i,j,L
integer(i_kind) ndatay,ndatax
integer(i_kind) g_ind
integer(i_kind), parameter :: DIR_S=1, DIR_N=2, DIR_E=3, DIR_W=4
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!
! Limit communications to generation one
!
          g_ind=1

          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)

          imax = im_in
          jmax = jm_in

!-----------------------------------------------------------------------
      ndatay = km_in*imax*nby*Lm_in
      ndatax = km_in*(jmax+2*nby)*nbx*Lm_in

      sHandle(:) = MPI_REQUEST_NULL
      rHandle(:) = MPI_REQUEST_NULL

!
!  Post receives toward SOUTH and NORTH
!
      if( itarg_s >= 0 ) then
        call MPI_IRECV( W(1,1,1-nby,1), ndatay, dtype, itarg_s, itarg_s, &
                        mpi_comm_comp, rHandle(DIR_S), ierr )
      end if

      if( itarg_n >= 0 ) then
        call MPI_IRECV( W(1,1,jmax+1,1), ndatay, dtype, itarg_n, itarg_n, &
                        mpi_comm_comp, rHandle(DIR_N), ierr )
      end if

!
!  Post receives toward WEST and EAST
!
      if( itarg_w >= 0 ) then
        call MPI_IRECV( W(1,1-nbx,1-nby,1), ndatax, dtype, itarg_w, itarg_w, &
                        mpi_comm_comp, rHandle(DIR_W), ierr )
      end if

      if( itarg_e >= 0 ) then
        call MPI_IRECV( W(1,imax+1,1-nby,1), ndatax, dtype, itarg_e, itarg_e, &
                        mpi_comm_comp, rHandle(DIR_E), ierr )
      end if

!
!  SEND boundaries toward SOUTH and NORTH
!
      if( itarg_s >= 0 ) then
        call MPI_ISEND( W(1,1,1,1), ndatay, dtype, itarg_s, mype, &
                        mpi_comm_comp, sHandle(DIR_S), ierr )
      end if

      if( itarg_n >= 0 ) then
        call MPI_ISEND( W(1,1,jmax-nby+1,1), ndatay, dtype, itarg_n, mype, &
                        mpi_comm_comp, sHandle(DIR_N), ierr )
      end if

!
!  SEND extended boundaries toward WEST and EAST
!
      if( itarg_w >= 0 ) then
        call MPI_ISEND( W(1,1,1-nby,1), ndatax, dtype, itarg_w, mype, &
                        mpi_comm_comp, sHandle(DIR_W), ierr )
      end if

      if( itarg_e >= 0 ) then
        call MPI_ISEND( W(1,imax-nbx+1,1-nby,1), ndatax, dtype, itarg_e, mype, &
                        mpi_comm_comp, sHandle(DIR_E), ierr )
      end if

!
!  Complete non-blocking receives and sends
!
      if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
      end if

      if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
      end if

!
!  Apply physical boundary conditions where needed
!
   if( lsouth ) then
     do L=1,Lm_in
       do j=1,nby
         do i=1,imax
           W(:,i,-nby+j,L)=W(:,i,nby+1-j,L)
         end do
       end do
     end do
   end if

   if( lnorth ) then
     do L=1,Lm_in
       do j=1,nby
         do i=1,imax
           W(:,i,jmax+j,L)=W(:,i,jmax+1-j,L)
         end do
       end do
     end do
   end if

   if( lwest ) then
     do L=1,Lm_in
       do j=1-nby,jmax+nby
         do i=1,nbx
           W(:,-nbx+i,j,L)= W(:,nbx+1-i,j,L)
         end do
       end do
     end do
   end if

   if( least ) then
     do L=1,Lm_in
       do j=1-nby,jmax+nby
         do i=1,nbx
           W(:,imax+i,j,L)=W(:,imax+1-i,j,L)
         end do
       end do
     end do
   end if

!------------------------------------------------------------------
endsubroutine boco_3d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine boco_3d_gh &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies (nbx,nby) lines of halos in (x,y) directions assuming       !
! mirror boundary conditions. Version for high generations             !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!**********************************************************************!
(this,W,km_in,im_in,jm_in,Lm_in,nbx,nby,nbz,Fimax_in,Fjmax_in,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,Lm_in,nbx,nby,nbz,mygen_min,mygen_max
real(r_kind),dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby,1-nbz:Lm_in+nbz),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
!-----------------------------------------------------------------------
integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth
integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,i,j,L
integer(i_kind) ndatay,ndatax
integer(i_kind) g_ind,g
logical l_sidesend
integer(i_kind), parameter :: DIR_S=1, DIR_N=2, DIR_E=3, DIR_W=4
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!
       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2
         g = my_hgen
         l_sidesend=.true.
       else
         l_sidesend=.false.
       endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain
!
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)

          if(least) then
            imax = Fimax_in(g)
          else
            imax = im_in
          endif
          if(lnorth) then
            jmax = Fjmax_in(g)
          else
            jmax = jm_in
          endif

!-----------------------------------------------------------------------
      ndatay = km_in*imax*nby*Lm_in
      ndatax = km_in*(jmax+2*nby)*nbx*Lm_in

      sHandle(:) = MPI_REQUEST_NULL
      rHandle(:) = MPI_REQUEST_NULL

!
!  Post receives toward SOUTH and NORTH
!
      if( itarg_s >= 0 ) then
        call MPI_IRECV( W(1,1,1-nby,1), ndatay, dtype, itarg_s, itarg_s, &
                        mpi_comm_work, rHandle(DIR_S), ierr )
      end if

      if( itarg_n >= 0 ) then
        call MPI_IRECV( W(1,1,jmax+1,1), ndatay, dtype, itarg_n, itarg_n, &
                        mpi_comm_work, rHandle(DIR_N), ierr )
      end if

!
!  Post receives toward WEST and EAST
!
      if( itarg_w >= 0 ) then
        call MPI_IRECV( W(1,1-nbx,1-nby,1), ndatax, dtype, itarg_w, itarg_w, &
                        mpi_comm_work, rHandle(DIR_W), ierr )
      end if

      if( itarg_e >= 0 ) then
        call MPI_IRECV( W(1,imax+1,1-nby,1), ndatax, dtype, itarg_e, itarg_e, &
                        mpi_comm_work, rHandle(DIR_E), ierr )
      end if

!
!  SEND boundaries to SOUTH and NORTH
!
      if( itarg_s >= 0 ) then
        call MPI_ISEND( W(1,1,1,1), ndatay, dtype, itarg_s, mype, &
                        mpi_comm_work, sHandle(DIR_S), ierr )
      end if

      if( itarg_n >= 0 ) then
        call MPI_ISEND( W(1,1,jmax-nby+1,1), ndatay, dtype, itarg_n, mype, &
                        mpi_comm_work, sHandle(DIR_N), ierr )
      end if

!
!  SEND extended boundaries to WEST and EAST
!
      if( itarg_w >= 0 ) then
        call MPI_ISEND( W(1,1,1-nby,1), ndatax, dtype, itarg_w, mype, &
                        mpi_comm_work, sHandle(DIR_W), ierr )
      end if

      if( itarg_e >= 0 ) then
        call MPI_ISEND( W(1,imax-nbx+1,1-nby,1), ndatax, dtype, itarg_e, mype, &
                        mpi_comm_work, sHandle(DIR_E), ierr )
      end if

!
!  Complete non-blocking receives and sends
!
      if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
      end if

      if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
      end if

!
!  Apply physical boundary conditions where needed
!
   if( lsouth ) then
     do L=1,Lm_in
       do j=1,nby
         do i=1,imax
           W(:,i,-nby+j,L)=W(:,i,nby+1-j,L)
         end do
       end do
     end do
   end if

   if( lnorth ) then
     do L=1,Lm_in
       do j=1,nby
         do i=1,imax
           W(:,i,jmax+j,L)=W(:,i,jmax+1-j,L)
         end do
       end do
     end do
   end if

   if( lwest ) then
     do L=1,Lm_in
       do j=1-nby,jmax+nby
         do i=1,nbx
           W(:,-nbx+i,j,L)= W(:,nbx+1-i,j,L)
         end do
       end do
     end do
   end if

   if( least ) then
     do L=1,Lm_in
       do j=1-nby,jmax+nby
         do i=1,nbx
           W(:,imax+i,j,L)=W(:,imax+1-i,j,L)
         end do
       end do
     end do
   end if

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
endsubroutine boco_3d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine bocoT_3d_g1 &
!***********************************************************************
!                                                                      !
! Adjoint of side sending subroutine:                                  !
! Supplies (nbx,nby) lines of halos in (x,y) directions, including     !
! values at the edges of the subdomains and assuming mirror boundary   !
! conditions just for generation 1                                     !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(this,W,km_in,im_in,jm_in,Lm_in,nbx,nby,nbz,Fimax_in,Fjmax_in)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,Lm_in,nbx,nby,nbz
real(r_kind), dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby,1-nbz:Lm_in+nbz),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:,:) :: rBuf_W,rBuf_E,rBuf_S,rBuf_N
integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical lwest,least,lsouth,lnorth                                       
integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,i,j,L
integer(i_kind) ndatay,ndatax
integer(i_kind) g_ind
integer(i_kind), parameter :: DIR_S=1, DIR_N=2, DIR_E=3, DIR_W=4
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!
! Limit communications to selected number of generations
!

         g_ind=1
!
! from mg_domain
!
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)

          imax = im_in    
          jmax = jm_in


!----------------------------------------------------------------------
      ndatay =km_in*imax*nby*Lm_in
      ndatax =km_in*(jmax+2*nby)*nbx*Lm_in

      sHandle(:) = MPI_REQUEST_NULL
      rHandle(:) = MPI_REQUEST_NULL

!
! RECEIVE extended halos from EAST and WEST
!
      if( itarg_w >= 0 ) then
        allocate( rBuf_W(km_in,nbx,1-nby:jmax+nby,Lm_in), stat = ierr )
        call MPI_IRECV( rBuf_W, ndatax, dtype, itarg_w, itarg_w, &
                        mpi_comm_comp, rHandle(DIR_W), ierr )
      end if

      if( itarg_e >= 0 ) then
        allocate( rBuf_E(km_in,nbx,1-nby:jmax+nby,Lm_in), stat = ierr )
        call MPI_IRECV( rBuf_E, ndatax, dtype, itarg_e, itarg_e, &
                        mpi_comm_comp, rHandle(DIR_E), ierr )
      end if

!
! RECEIVE boundaries from SOUTH and NORTH
!
      if( itarg_s >= 0 ) then
        allocate( rBuf_S(km_in,imax,nby,Lm_in), stat = ierr )
        call MPI_IRECV( rBuf_S, ndatay, dtype, itarg_s, itarg_s, &
                        mpi_comm_comp, rHandle(DIR_S), ierr )
      end if

      if( itarg_n >= 0 ) then
        allocate( rBuf_N(km_in,imax,nby,Lm_in), stat = ierr )
        call MPI_IRECV( rBuf_N, ndatay, dtype, itarg_n, itarg_n, &
                        mpi_comm_comp, rHandle(DIR_N), ierr )
      end if

!
! SEND extended halos toward WEST and EAST
!
      if( itarg_w >= 0) then
        call MPI_ISEND( W(1,1-nbx,1-nby,1), ndatax, dtype, itarg_w, mype, &
                        mpi_comm_comp, sHandle(DIR_W), ierr )
      end if

      if( itarg_e >= 0 ) then
        call MPI_ISEND( W(1,imax+1,1-nby,1), ndatax, dtype, itarg_e, mype, &
                        mpi_comm_comp, sHandle(DIR_E), ierr )
      end if

!
! SEND boundaries SOUTH and NORTH
!
      if( itarg_s >= 0 ) then
        call MPI_ISEND( W(1,1,1-nby,1), ndatay, dtype, itarg_s, mype, &
                        mpi_comm_comp, sHandle(DIR_S), ierr )
      end if

      if( itarg_n >= 0 ) then
        call MPI_ISEND( W(1,1,jmax+1,1), ndatay, dtype, itarg_n, mype, &
                        mpi_comm_comp, sHandle(DIR_N), ierr )
      end if

!
!  Complete non-blocking operations
!
      if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
      end if

      if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
      end if

!
! Assign received halos from WEST and EAST to interior of domains
!
   if(lwest) then
     do L=1,Lm_in
       do j=1-nby,jmax+nby
         do i=1,nbx
           W(:,i,j,L)= W(:,i,j,L)+W(:,1-i,j,L)
         end do
       end do
     end do
   else if( allocated(rBuf_W) ) then
     do L=1,Lm_in
       do j=1-nby,jmax+nby
         do i=1,nbx
           W(:,i,j,L)= W(:,i,j,L)+rBuf_W(:,i,j,L)
         end do
       end do
     end do
   end if

   if(least) then
     do L=1,Lm_in
       do j=1-nby,jmax+nby
         do i=1,nbx
           W(:,imax-nbx+i,j,L)= W(:,imax-nbx+i,j,L)+W(:,imax+1+nbx-i,j,L)
         end do
       end do
     end do
   else if( allocated(rBuf_E) ) then
     do L=1,Lm_in
       do j=1-nby,jmax+nby
         do i=1,nbx
           W(:,imax-nbx+i,j,L)= W(:,imax-nbx+i,j,L)+rBuf_E(:,i,j,L)
         end do
       end do
     end do
   end if

!
! ASSIGN received values from SOUTH and NORTH
!
   if(lsouth) then
     do L=1,Lm_in
       do j=1,nby
         do i=1,imax
           W(:,i,j,L)= W(:,i,j,L)+W(:,i,1-j,L)
         end do
       end do
     end do
   else if( allocated(rBuf_S) ) then
     do L=1,Lm_in
       do j=1,nby
         do i=1,imax
           W(:,i,j,L)= W(:,i,j,L)+rBuf_S(:,i,j,L)
         end do
       end do
     end do
   end if

   if(lnorth) then
     do L=1,Lm_in
       do j=1,nby
         do i=1,imax
           W(:,i,jmax-nby+j,L)= W(:,i,jmax-nby+j,L)+W(:,i,jmax+1+nby-j,L)
         end do
       end do
     end do
   else if( allocated(rBuf_N) ) then
     do L=1,Lm_in
       do j=1,nby
         do i=1,imax
           W(:,i,jmax-nby+j,L)= W(:,i,jmax-nby+j,L)+rBuf_N(:,i,j,L)
         end do
       end do
     end do
   end if

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!
      if( allocated(rBuf_W) ) deallocate( rBuf_W, stat = ierr)
      if( allocated(rBuf_E) ) deallocate( rBuf_E, stat = ierr)
      if( allocated(rBuf_S) ) deallocate( rBuf_S, stat = ierr)
      if( allocated(rBuf_N) ) deallocate( rBuf_N, stat = ierr)

!-----------------------------------------------------------------------
endsubroutine bocoT_3d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine bocoT_3d_gh &
!***********************************************************************
!                                                                      !
!  Supply n-lines inside of domains, including edges, with halos from  !
!  the surrounding domains.  Assume mirror boundary conditions at the  !
!  boundaries of the domain. For high multigrid generations.           !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(this,W,km_in,im_in,jm_in,Lm_in,nbx,nby,nbz,Fimax_in,Fjmax_in,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,Lm_in,nbx,nby,nbz,mygen_min,mygen_max
real(r_kind), dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby,1-nbz:Lm_in+nbz),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:,:):: rBuf_W,rBuf_E,rBuf_S,rBuf_N
integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical lwest,least,lsouth,lnorth                                       
integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,i,j,L
integer(i_kind) ndatay,ndatax
integer(i_kind) g_ind,g
logical l_sidesend
integer(i_kind), parameter :: DIR_S=1, DIR_N=2, DIR_E=3, DIR_W=4
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!
! Limit comminications to selected number of generations
!


       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2
         g = my_hgen
         l_sidesend=.true.
       else 
         l_sidesend=.false.
       endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 

          if(least) then
            imax = Fimax_in(g)
          else 
            imax = im_in       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax_in(g)
          else  
            jmax = jm_in
          endif


!----------------------------------------------------------------------
      ndatax =km_in*(jmax+2*nby)*nbx*Lm_in
      ndatay =km_in*imax*nby*Lm_in

      sHandle(:) = MPI_REQUEST_NULL
      rHandle(:) = MPI_REQUEST_NULL

!
! RECEIVE extended halos from EAST and WEST
!
      if( itarg_w >= 0 ) then
        allocate( rBuf_W(km_in,nbx,1-nby:jmax+nby,Lm_in), stat = ierr )
        call MPI_IRECV( rBuf_W, ndatax, dtype, itarg_w, itarg_w, &
                        mpi_comm_work, rHandle(DIR_W), ierr )
      end if
      if( itarg_e >= 0 ) then
        allocate( rBuf_E(km_in,nbx,1-nby:jmax+nby,Lm_in), stat = ierr )
        call MPI_IRECV( rBuf_E, ndatax, dtype, itarg_e, itarg_e, &
                        mpi_comm_work, rHandle(DIR_E), ierr )
      end if

!
! RECEIVE boundaries from SOUTH and NORTH
!
      if( itarg_s >= 0 ) then
        allocate( rBuf_S(km_in,imax,nby,Lm_in), stat = ierr )
        call MPI_IRECV( rBuf_S, ndatay, dtype, itarg_s, itarg_s, &
                        mpi_comm_work, rHandle(DIR_S), ierr )
      end if
      if( itarg_n >= 0 ) then
        allocate( rBuf_N(km_in,imax,nby,Lm_in), stat = ierr )
        call MPI_IRECV( rBuf_N, ndatay, dtype, itarg_n, itarg_n, &
                        mpi_comm_work, rHandle(DIR_N), ierr )
      end if

!
! SEND extended halos toward WEST and EAST
!
      if( itarg_w >= 0) then
        call MPI_ISEND( W(1,1-nbx,1-nby,1), ndatax, dtype, itarg_w, mype, &
                        mpi_comm_work, sHandle(DIR_W), ierr )
      end if
      if( itarg_e >= 0 ) then
        call MPI_ISEND( W(1,imax+1,1-nby,1), ndatax, dtype, itarg_e, mype, &
                        mpi_comm_work, sHandle(DIR_E), ierr )
      end if

!
! SEND boundaries SOUTH and NORTH
!
      if( itarg_s >= 0 ) then
        call MPI_ISEND( W(1,1,1-nby,1), ndatay, dtype, itarg_s, mype, &
                        mpi_comm_work, sHandle(DIR_S), ierr )
      end if
      if( itarg_n >= 0 ) then
        call MPI_ISEND( W(1,1,jmax+1,1), ndatay, dtype, itarg_n, mype, &
                        mpi_comm_work, sHandle(DIR_N), ierr )
      end if

!
!  Complete non-blocking operations
!
      if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
      end if

      if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
      end if

!
! Assign received halos from WEST and EAST to interior of domains
!
   if(lwest) then
     do L=1,Lm_in
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,i,j,L)= W(:,i,j,L)+W(:,1-i,j,L)
     end do
     end do
     end do
   else if( allocated(rBuf_W) ) then
     do L=1,Lm_in
     do j=1-nby,jmax+nby
     do i=1,nbx
      W(:,i,j,L)= W(:,i,j,L)+rBuf_W(:,i,j,L)
     end do
     end do
     end do
   end if

   if(least) then
     do L=1,Lm_in
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax-nbx+i,j,L)= W(:,imax-nbx+i,j,L)+W(:,imax+1+nbx-i,j,L)
     end do
     end do
     end do
   else if( allocated(rBuf_E) ) then
     do L=1,Lm_in
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax-nbx+i,j,L)= W(:,imax-nbx+i,j,L)+rBuf_E(:,i,j,L)
     end do
     end do
     end do
   end if

!
! Assign received values from SOUTH and NORTH
!
   if(lsouth) then
     do L=1,Lm_in
     do j=1,nby
     do i=1,imax
       W(:,i,j,L)= W(:,i,j,L)+W(:,i,1-j,L)
     end do
     end do
     end do
   else if( allocated(rBuf_S) ) then
     do L=1,Lm_in
     do j=1,nby
     do i=1,imax
       W(:,i,j,L)= W(:,i,j,L)+rBuf_S(:,i,j,L)
     end do
     end do
     end do
   end if

   if(lnorth) then
     do L=1,Lm_in
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j,L)= W(:,i,jmax-nby+j,L)+W(:,i,jmax+1+nby-j,L)
     end do
     end do
     end do
   else if( allocated(rBuf_N) ) then
     do L=1,Lm_in
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j,L)= W(:,i,jmax-nby+j,L)+rBuf_N(:,i,j,L)
     end do
     end do
     end do
   end if

!
!                           DEALLOCATE rBufferes
!
      if( allocated(rBuf_W) ) deallocate( rBuf_W, stat = ierr)
      if( allocated(rBuf_E) ) deallocate( rBuf_E, stat = ierr)
      if( allocated(rBuf_S) ) deallocate( rBuf_S, stat = ierr)
      if( allocated(rBuf_N) ) deallocate( rBuf_N, stat = ierr)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
endsubroutine bocoT_3d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine upsend_all_g1 &
!***********************************************************************
!                                                                      !
!         Upsend data from generation one to generation two            !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(this,Harray,Warray,km_in)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in
real(r_kind), dimension(km_in,1:this%imL,1:this%jmL),intent(in):: Harray
real(r_kind), dimension(km_in,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(out):: Warray
!-----------------------------------------------------------------------
real(r_kind),dimension(1:km_in,1:this%imL,1:this%jmL):: dBuf_SW
real(r_kind),dimension(1:km_in,1:this%imL,1:this%jmL):: dBuf_SE
real(r_kind),dimension(1:km_in,1:this%imL,1:this%jmL):: dBuf_NW
real(r_kind),dimension(1:km_in,1:this%imL,1:this%jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,ndata,i,j
integer(i_kind) nebpe

integer(i_kind):: mygen_dn,mygen_up
logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne
integer(i_kind):: itarg_up
integer:: g_ind
integer(i_kind), parameter :: DIR_SW=1, DIR_SE=2, DIR_NW=3, DIR_NE=4
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"

!-----------------------------------------------------------------------
   mygen_dn=1
   mygen_up=2
!
! Define generational flags
!
       g_ind=1

       lsendup_sw=Flsendup_sw(g_ind)
       lsendup_se=Flsendup_se(g_ind)
       lsendup_nw=Flsendup_nw(g_ind)
       lsendup_ne=Flsendup_ne(g_ind)

       itarg_up=Fitarg_up(g_ind)

!-----------------------------------------------------------------------

   if(my_hgen==mygen_up) then
      Warray(:,:,:) = 0.0d0
   endif

     ndata =km_in*imL*jmL

     sHandle(:) = MPI_REQUEST_NULL
     rHandle(:) = MPI_REQUEST_NULL

!
! --- Send data to SW portion of processors at higher generation
!
      if( lsendup_sw ) then
        nebpe = itarg_up
        if(nebpe == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_SW(:,i,j) = Harray(:,i,j)
             enddo
             enddo
        else
        call MPI_ISEND( Harray, ndata, dtype, nebpe, mype,  &
                       mpi_comm_comp, sHandle(DIR_SW), ierr)
        endif
      end if
!
! --- Receive SW portion of data at higher generation
!
      if( my_hgen==mygen_up .and. itargdn_sw >= 0 ) then
        nebpe = itargdn_sw
        if(nebpe /= mype) then
          call MPI_IRECV( dBuf_SW, ndata, dtype, nebpe, nebpe,          &
                          mpi_comm_comp, rHandle(DIR_SW), ierr)
        endif
      end if
!
! --- Send data to SE portion of processors at higher generation
!
      if( lsendup_se ) then
        nebpe = itarg_up
        if(nebpe == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_SE(:,i,j) = Harray(:,i,j)
             enddo
             enddo
        else
        call MPI_ISEND( Harray, ndata, dtype, nebpe, mype, &
                       mpi_comm_comp, sHandle(DIR_SE), ierr)
        endif
      end if
!
! --- Receive SE portion of data at higher generation
!
      if( my_hgen==mygen_up .and. itargdn_se >= 0 ) then
        nebpe = itargdn_se
        if(nebpe /= mype) then
          call MPI_IRECV( dBuf_SE, ndata, dtype, nebpe, nebpe,          &
                          mpi_comm_comp, rHandle(DIR_SE), ierr)
        endif
      end if
!
! --- Send data to NW portion of processors at higher generation
!
      if( lsendup_nw ) then
        nebpe = itarg_up
        if(nebpe == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_NW(:,i,j) = Harray(:,i,j)
             enddo
             enddo
        else
         call MPI_ISEND( Harray, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(DIR_NW), ierr)
        end if
      end if
!
! --- Receive NW portion of data at higher generation
!
      if( my_hgen==mygen_up .and. itargdn_nw >= 0 ) then
        nebpe = itargdn_nw
        if(nebpe /= mype) then
          call MPI_IRECV( dBuf_NW, ndata, dtype, nebpe, nebpe,          &
                          mpi_comm_comp, rHandle(DIR_NW), ierr)
        endif
      end if
!
! --- Send data to NE portion of processors at higher generation
!
      if( lsendup_ne ) then
        nebpe = itarg_up
        if(nebpe == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_NE(:,i,j) = Harray(:,i,j)
             enddo
             enddo
        else
        call MPI_ISEND( Harray, ndata, dtype, nebpe, mype, &
                      mpi_comm_comp, sHandle(DIR_NE), ierr)
        endif
      end if
!
! --- Receive NE portion of data at higher generation
!
      if( my_hgen==mygen_up .and. itargdn_ne >= 0 ) then
        nebpe = itargdn_ne
        if(nebpe /= mype) then
          call MPI_IRECV( dBuf_NE, ndata, dtype, nebpe, nebpe,          &
                         mpi_comm_comp, rHandle(DIR_NE), ierr)
        endif
      end if

      if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
      end if
      if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
      end if

      if( my_hgen==mygen_up .and. itargdn_sw >= 0 ) then
             do j=1,jmL
             do i=1,imL
                Warray(:,i,j)=dBuf_SW(:,i,j)
             enddo
             enddo
      end if

      if( my_hgen==mygen_up .and. itargdn_se >= 0 ) then
             do j=1,jmL
             do i=1,imL
               Warray(:,imL+i,j)=dBuf_SE(:,i,j)
             enddo
             enddo
      end if

      if( my_hgen==mygen_up .and. itargdn_nw >= 0 ) then
             do j=1,jmL
             do i=1,imL
               Warray(:,i,jmL+j)=dBuf_NW(:,i,j)
             enddo
             enddo
      end if

      if( my_hgen==mygen_up .and. itargdn_ne >= 0 ) then
             do j=1,jmL
             do i=1,imL
               Warray(:,imL+i,jmL+j)=dBuf_NE(:,i,j)
             enddo
             enddo
      end if

!-----------------------------------------------------------------------
endsubroutine upsend_all_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine upsend_all_gh &
!***********************************************************************
!                                                                      *
!         Upsend data from one grid generation to another              *
!         (Just for high grid generations)                             *
!                                                                      !
!                       - offset version -                             !
!                                                                      *
!***********************************************************************
(this,Harray,Warray,km_in,mygen_dn,mygen_up)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in
real(r_kind), dimension(km_in,1:this%imL,1:this%jmL),intent(in):: Harray
real(r_kind), dimension(km_in,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(out):: Warray
integer(i_kind),intent(in):: mygen_dn,mygen_up
!-----------------------------------------------------------------------
real(r_kind),dimension(1:km_in,1:this%imL,1:this%jmL):: dBuf_SW
real(r_kind),dimension(1:km_in,1:this%imL,1:this%jmL):: dBuf_SE
real(r_kind),dimension(1:km_in,1:this%imL,1:this%jmL):: dBuf_NW
real(r_kind),dimension(1:km_in,1:this%imL,1:this%jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,ndata,i,j
integer(i_kind) nebpe

logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne
integer(i_kind):: itarg_up
integer:: g_ind
integer(i_kind), parameter :: DIR_SW=1, DIR_SE=2, DIR_NW=3, DIR_NE=4
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"

!-----------------------------------------------------------------------
!
! Define generational flags
!
 
       g_ind=2 

       lsendup_sw=Flsendup_sw(g_ind).and.(my_hgen==mygen_dn)
       lsendup_se=Flsendup_se(g_ind).and.(my_hgen==mygen_dn)
       lsendup_nw=Flsendup_nw(g_ind).and.(my_hgen==mygen_dn)
       lsendup_ne=Flsendup_ne(g_ind).and.(my_hgen==mygen_dn)

       itarg_up=Fitarg_up(g_ind)                                          


!-----------------------------------------------------------------------

   if(my_hgen==mygen_up) then
      Warray(:,:,:) = 0.0d0
   endif

     ndata =km_in*imL*jmL

     sHandle(:) = MPI_REQUEST_NULL
     rHandle(:) = MPI_REQUEST_NULL

      if(  lsendup_sw ) then
        nebpe = itarg_up
        call MPI_ISEND( Harray, ndata, dtype, nebpe, mype,  &
                       mpi_comm_work, sHandle(DIR_SW), ierr)
      end if

      if( my_hgen==mygen_up .and. itargdn_sw >= 0 ) then
        nebpe = itargdn_sw
        call MPI_IRECV( dBuf_SW, ndata, dtype, nebpe, nebpe, &
                       mpi_comm_work, rHandle(DIR_SW), ierr)
      endif

      if( lsendup_se ) then
        nebpe = itarg_up
        call MPI_ISEND( Harray, ndata, dtype, nebpe, mype, &
                       mpi_comm_work, sHandle(DIR_SE), ierr)
      end if

      if( my_hgen==mygen_up .and. itargdn_se >= 0 ) then
        nebpe = itargdn_se
        call MPI_IRECV( dBuf_SE, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(DIR_SE), ierr)
      endif

      if( lsendup_nw ) then
        nebpe = itarg_up
         call MPI_ISEND( Harray, ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(DIR_NW), ierr)
      end if

      if( my_hgen==mygen_up .and. itargdn_nw >= 0 ) then
        nebpe = itargdn_nw
        call MPI_IRECV( dBuf_NW, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(DIR_NW), ierr)
      end if

      if( lsendup_ne ) then
        nebpe = itarg_up
        call MPI_ISEND( Harray, ndata, dtype, nebpe, mype, &
                       mpi_comm_work, sHandle(DIR_NE), ierr)
      end if

      if( my_hgen==mygen_up .and. itargdn_ne >= 0 ) then
        nebpe = itargdn_ne
        call MPI_IRECV( dBuf_NE, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(DIR_NE), ierr)
      endif

      if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
      end if
      if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
      end if

      if( my_hgen==mygen_up .and. itargdn_sw >= 0 ) then
             do j=1,jmL
             do i=1,imL
               Warray(:,i,j)=dBuf_SW(:,i,j)
             enddo
             enddo
      endif

      if( my_hgen==mygen_up .and. itargdn_se >= 0 ) then
             do j=1,jmL
             do i=1,imL
               Warray(:,imL+i,j)=dBuf_SE(:,i,j)
             enddo
             enddo
      endif

      if( my_hgen==mygen_up .and. itargdn_nw >= 0 ) then
             do j=1,jmL
             do i=1,imL
               Warray(:,i,jmL+j)=dBuf_NW(:,i,j)
             enddo
             enddo
      end if

      if( my_hgen==mygen_up .and. itargdn_ne >= 0 ) then
             do j=1,jmL
             do i=1,imL
               Warray(:,imL+i,jmL+j)=dBuf_NE(:,i,j)
             enddo
             enddo
      endif

!-----------------------------------------------------------------------
endsubroutine upsend_all_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine downsend_all_gh &
!***********************************************************************
!                                                                      *
!         Downsending data from low resolution pes    (mygen_up)       *
!         to the concurent high-resolution pes        (mygen_dn)       *
!         and add the existing and the recevied values                 *
!                                                                      !
!                       - offset version -                             !
!                                                                      *
!***********************************************************************
(this,Warray,Harray,km_in,mygen_up,mygen_dn)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in
real(r_kind), dimension(km_in,1:this%im,1:this%jm),intent(in):: Warray
real(r_kind), dimension(km_in,1:this%imL,1:this%jmL),intent(out):: Harray
integer, intent(in):: mygen_up,mygen_dn
!-----------------------------------------------------------------------
real(r_kind),dimension(1:km_in,1:this%imL,1:this%jmL):: dBuf_SW
real(r_kind),dimension(1:km_in,1:this%imL,1:this%jmL):: dBuf_SE
real(r_kind),dimension(1:km_in,1:this%imL,1:this%jmL):: dBuf_NW
real(r_kind),dimension(1:km_in,1:this%imL,1:this%jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,ndata,i,j
integer(i_kind) nebpe

logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne  
integer(i_kind):: itarg_up                                           
integer(i_kind):: g_ind
integer(i_kind), parameter :: DIR_SW=1, DIR_SE=2, DIR_NW=3, DIR_NE=4
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------

     Harray(:,:,:) = 0.0d0
!
! Define generational flags
!

         g_ind=2
       lsendup_sw=Flsendup_sw(g_ind).and.(my_hgen==mygen_dn)
       lsendup_se=Flsendup_se(g_ind).and.(my_hgen==mygen_dn)
       lsendup_nw=Flsendup_nw(g_ind).and.(my_hgen==mygen_dn)
       lsendup_ne=Flsendup_ne(g_ind).and.(my_hgen==mygen_dn)

       itarg_up=Fitarg_up(g_ind)

       ndata =km_in*imL*jmL

      dBuf_SW = 0.0d0
      dBuf_SE = 0.0d0
      dBuf_NW = 0.0d0
      dBuf_NE = 0.0d0

      sHandle(:) = MPI_REQUEST_NULL
      rHandle(:) = MPI_REQUEST_NULL

  if(my_hgen==mygen_up .and. itargdn_sw >= 0 ) then
        nebpe = itargdn_sw
        if(nebpe == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_SW(:,i,j) = Warray(:,i,j)
             enddo
             enddo
        else
        call MPI_ISEND( Warray(1,1,1), ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(DIR_SW), ierr )
        end if
  endif

  if(my_hgen==mygen_up .and.  itargdn_se >= 0 ) then
        nebpe = itargdn_se
        if(nebpe == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_SE(:,i,j) = Warray(:,imL+i,j)
             enddo
             enddo
        else
        call MPI_ISEND( Warray(1,imL+1,1), ndata, dtype, nebpe, mype,  &
                       mpi_comm_work, sHandle(DIR_SE), ierr )
        end if
  endif

  if(my_hgen==mygen_up .and. itargdn_nw >= 0 ) then
        nebpe = itargdn_nw
        if(nebpe == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_NW(:,i,j) = Warray(:,i,jmL+j)
             enddo
             enddo
        else
        call MPI_ISEND( Warray(1,1,jmL+1), ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(DIR_NW), ierr )
        end if
  endif

  if(my_hgen==mygen_up .and. itargdn_ne >= 0 ) then
        nebpe = itargdn_ne
        if(nebpe == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_NE(:,i,j) = Warray(:,imL+i,jmL+j)
             enddo
             enddo
        else
        call MPI_ISEND( Warray(1,imL+1,jmL+1), ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(DIR_NE), ierr )
        end if
  endif

  if( lsendup_sw ) then
        nebpe = itarg_up
        call MPI_IRECV( dBuf_SW, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_work, rHandle(DIR_SW), ierr)
  endif

  if( lsendup_se ) then
        nebpe = itarg_up
        call MPI_IRECV( dBuf_SE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_work, rHandle(DIR_SE), ierr)
  end if

  if( lsendup_nw ) then
        nebpe = itarg_up
        call MPI_IRECV( dBuf_NW, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_work, rHandle(DIR_NW), ierr)
  end if

  if( lsendup_ne ) then
        nebpe = itarg_up
        call MPI_IRECV( dBuf_NE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_work, rHandle(DIR_NE), ierr)
  end if

  if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
  end if
  if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
  end if

  if( lsendup_sw ) then
             do j=1,jmL
             do i=1,imL
               Harray(:,i,j)=dBuf_SW(:,i,j)  
             enddo
             enddo
  endif

  if( lsendup_se ) then
             do j=1,jmL
             do i=1,imL
               Harray(:,i,j)=dBuf_SE(:,i,j)
             enddo
             enddo
  endif

  if( lsendup_nw ) then
             do j=1,jmL
             do i=1,imL
               Harray(:,i,j)=dBuf_NW(:,i,j)
             enddo
             enddo
  endif

  if( lsendup_ne ) then
             do j=1,jmL
             do i=1,imL
               Harray(:,i,j)=dBuf_NE(:,i,j)
             enddo
             enddo
  endif

!-----------------------------------------------------------------------
endsubroutine downsend_all_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine downsend_all_g2 &
!***********************************************************************
!                                                                      *
!         Downsending data from low resolution pes    (mygen_up)       *
!         to the concurent high-resolution pes        (mygen_dn)       *
!         and add the existing and the recevied values                 *
!                                                                      *
!                       - offset version -                             *
!                                                                      *
!***********************************************************************
(this,Warray,Harray,km_in)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in
real(r_kind), dimension(km_in,1:this%im,1:this%jm),intent(in):: Warray
real(r_kind), dimension(km_in,1:this%imL,1:this%jmL),intent(out):: Harray
!-----------------------------------------------------------------------
real(r_kind),dimension(1:km_in,1:this%imL,1:this%jmL):: dBuf_SW
real(r_kind),dimension(1:km_in,1:this%imL,1:this%jmL):: dBuf_SE
real(r_kind),dimension(1:km_in,1:this%imL,1:this%jmL):: dBuf_NW
real(r_kind),dimension(1:km_in,1:this%imL,1:this%jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,ndata,i,j
integer(i_kind) nebpe

logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne  
integer:: mygen_up,mygen_dn
integer(i_kind):: itarg_up                                           
integer(i_kind):: g_ind
integer(i_kind), parameter :: DIR_SW=1, DIR_SE=2, DIR_NW=3, DIR_NE=4
!-----------------------------------------------------------------------
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!
! Define generational flags
!
    mygen_up=2
    mygen_dn=1

         g_ind=1
       lsendup_sw=Flsendup_sw(g_ind)
       lsendup_se=Flsendup_se(g_ind)
       lsendup_nw=Flsendup_nw(g_ind)
       lsendup_ne=Flsendup_ne(g_ind)

       itarg_up=Fitarg_up(g_ind)

      ndata =km_in*imL*jmL

      dBuf_SW = 0.0d0
      dBuf_SE = 0.0d0
      dBuf_NW = 0.0d0
      dBuf_NE = 0.0d0

      sHandle(:) = MPI_REQUEST_NULL
      rHandle(:) = MPI_REQUEST_NULL

!
! Send data down to generation 1
!
      if(my_hgen==mygen_up) then

        if(itargdn_sw >= 0) then
          nebpe = itargdn_sw
          if(nebpe == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_SW(:,i,j) = Warray(:,i,j)
             enddo
             enddo
          else
             call MPI_ISEND( Warray(1,1,1), ndata, dtype, nebpe, mype,  &
                             mpi_comm_comp, sHandle(DIR_SW), ierr )
          end if
        end if

        if(itargdn_se >= 0) then
          nebpe = itargdn_se
          if(nebpe == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_SE(:,i,j) = Warray(:,imL+i,j)
             enddo
             enddo
          else
             call MPI_ISEND( Warray(1,imL+1,1), ndata, dtype, nebpe, mype,  &
                             mpi_comm_comp, sHandle(DIR_SE), ierr )
          end if
        end if

        if(itargdn_nw >= 0) then
          nebpe = itargdn_nw
          if(nebpe == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_NW(:,i,j) = Warray(:,i,jmL+j)
             enddo
             enddo
          else
             call MPI_ISEND( Warray(1,1,jmL+1), ndata, dtype, nebpe, mype,  &
                             mpi_comm_comp, sHandle(DIR_NW), ierr )
          end if
        end if

        if(itargdn_ne >= 0) then
          nebpe = itargdn_ne
          if(nebpe == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_NE(:,i,j) = Warray(:,imL+i,jmL+j)
             enddo
             enddo
          else
             call MPI_ISEND( Warray(1,imL+1,jmL+1), ndata, dtype, nebpe, mype,  &
                             mpi_comm_comp, sHandle(DIR_NE), ierr )
          end if
        end if

      end if

!
! Receive data on generation one
!
      if( lsendup_sw .and. mype /= itarg_up ) then
        nebpe = itarg_up
        call MPI_IRECV( dBuf_SW, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(DIR_SW), ierr )
      end if

      if( lsendup_se .and. mype /= itarg_up ) then
        nebpe = itarg_up
        call MPI_IRECV( dBuf_SE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(DIR_SE), ierr )
      end if

      if( lsendup_nw .and. mype /= itarg_up ) then
        nebpe = itarg_up
        call MPI_IRECV( dBuf_NW, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(DIR_NW), ierr )
      end if

      if( lsendup_ne .and. mype /= itarg_up ) then
        nebpe = itarg_up
        call MPI_IRECV( dBuf_NE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(DIR_NE), ierr )
      end if

      if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
      end if
      if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
      end if

!
! Assign received and prescribed values
!
      if( lsendup_sw ) then
             do j=1,jmL
             do i=1,imL
               Harray(:,i,j)=dBuf_SW(:,i,j)
             enddo
             enddo
      end if

      if( lsendup_se ) then
             do j=1,jmL
             do i=1,imL
               Harray(:,i,j)=dBuf_SE(:,i,j)
             enddo
             enddo
      end if

      if( lsendup_nw ) then
             do j=1,jmL
             do i=1,imL
               Harray(:,i,j)=dBuf_NW(:,i,j)
             enddo
             enddo
      end if

      if( lsendup_ne ) then
             do j=1,jmL
             do i=1,imL
               Harray(:,i,j)=dBuf_NE(:,i,j)
             enddo
             enddo
      end if

!-----------------------------------------------------------------------
endsubroutine downsend_all_g2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine bocox_2d_g1 &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nbx lines of halos in x direction assuming mirror boundary  !
! conditions at the end of domain. Version for generation 1            !                                             
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!**********************************************************************!
(this,W,km_in,im_in,jm_in,nbx,nby)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby
real(r_kind),dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:):: sBuf_E,sBuf_W             &
                                             ,rBuf_E,rBuf_W           

integer(i_kind) itarg_w,itarg_e,imax,jmax
logical:: lwest,least

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax
integer(i_kind) g_ind,g
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          g_ind = 1

          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)

          imax = im_in       
          jmax = jm_in


!-----------------------------------------------------------------------
      ndatax = km_in*jmax*nbx

!----------------------------------------------------------------------
!
! SEND extended boundaries toward WEST and EAST 
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km_in,nbx,1:jmax), stat = iaerr )

                do j=1,jmax
                  do i=1,nbx
                    sBuf_W(:,i,j) = W(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_comp, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km_in,nbx,1:jmax), stat = iaerr )

                do j=1,jmax
                  do i=1,nbx
                    sBuf_E(:,i,j) = W(:,imax-nbx+i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_comp, sHandle(2), isend)

      end if

!
! RECEIVE boundaries from EAST and WEST 
!

! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km_in,nbx,1:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km_in,nbx,1:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if


!
! Assign received values from EAST and WEST
!

! From west

   if(lwest) then

     do j=1,jmax
     do i=1,nbx
       W(:,-nbx+i,j)= W(:,nbx+1-i,j)
     end do
     end do

   else 

     do j=1,jmax
     do i=1,nbx
       W(:,-nbx+i,j)= rBuf_W(:,i,j)
     enddo
     enddo


   endif

! From east

   if(least) then

     do j=1,jmax
     do i=1,nbx
       W(:,imax+i,j)=W(:,imax+1-i,j)
     end do
     end do

   else 

     do j=1,jmax
     do i=1,nbx
       W(:,imax+i,j)=rBuf_E(:,i,j)
     enddo
     enddo

   endif


!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_e >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      endif
      if( itarg_w >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      endif

!
!                           DEALLOCATE sBufferes
!

      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
endsubroutine bocox_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine bocox_2d_gh &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nbx lines of halos in x direction assuming mirror boundary  !
! conditions at the end of domain. Version for high generations        !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!**********************************************************************!
(this,W,km_in,im_in,jm_in,nbx,nby,Fimax_in,Fjmax_in,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby,mygen_min,mygen_max
real(r_kind),dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:):: sBuf_E,sBuf_W             &
                                             ,rBuf_E,rBuf_W           

integer(i_kind) itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax
integer(i_kind) g_ind,g
logical l_sidesend
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!
 
       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2 
         g = my_hgen
         l_sidesend=.true.
       else
         l_sidesend=.false.
       endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 


          if(least) then
            imax = Fimax_in(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax_in(g)
          else  
            jmax = jm
          endif


!-----------------------------------------------------------------------
      ndatax = km_in*jmax*nbx

!
!  SEND halos to WEST and EASTH
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km_in,nbx,1:jmax), stat = iaerr )

                do j=1,jmax
                  do i=1,nbx
                    sBuf_W(:,i,j) = W(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km_in,nbx,1:jmax), stat = iaerr )

                do j=1,jmax
                  do i=1,nbx
                    sBuf_E(:,i,j) = W(:,imax-nbx+i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(2), isend)

      end if

!
!     RECEIVE extended boundaries from EAST and WEST
!

! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km_in,nbx,1:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km_in,nbx,1:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if

!
! Assign received values from  WEST and EAST
!

! From west

   if(lwest) then

     do j=1,jmax
     do i=1,nbx
       W(:,-nbx+i,j)= W(:,nbx+1-i,j)
     end do
     end do

   else 

     do j=1,jmax
     do i=1,nbx
       W(:,-nbx+i,j)= rBuf_W(:,i,j)
     enddo
     enddo


   endif

! From east

   if(least) then

     do j=1,jmax
     do i=1,nbx
       W(:,imax+i,j)=W(:,imax-i,j)
     end do
     end do

   else 

     do j=1,jmax
     do i=1,nbx
       W(:,imax+i,j)=rBuf_E(:,i,j)
     enddo
     enddo

   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_e >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      endif
      if( itarg_w >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      endif

!
!                           DEALLOCATE sBufferes
!

      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
endsubroutine bocox_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine bocoy_2d_g1 &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nby lines of halos in y direction assuming mirror boundary  !
! conditions at the end of domain. Version for generation 1            !                                             
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!**********************************************************************!
(this,W,km_in,im_in,jm_in,nbx,nby)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby
real(r_kind),dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:):: sBuf_N,sBuf_S             &
                                             ,rBuf_N,rBuf_S

integer(i_kind) itarg_n,itarg_s,imax,jmax
logical:: lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatay
integer(i_kind) g_ind,g
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          g_ind = 1

          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)

          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 

          imax = im_in       
          jmax = jm_in


!-----------------------------------------------------------------------
      ndatay = km_in*imax*nby


!
!  SEND boundaries toward SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km_in,1:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=1,imax
                    sBuf_S(:,i,j) = W(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_comp, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km_in,1:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=1,imax
                    sBuf_N(:,i,j)=W(:,i,jmax-nby+j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_comp, sHandle(1), isend)

      end if

!
! RECEIVE boundaries from NORTH and SOUTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km_in,1:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km_in,1:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if
!
! Assign received values from NORTH and SOUTH
!

! --- from NORTH ---

   if( lnorth) then

     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j)=W(:,i,jmax+1-j)
     enddo
     enddo

   else

     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j)=rBuf_N(:,i,j)
     enddo
     enddo

   endif

! From SOUTH

   if(lsouth) then

     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j)=W(:,i,nby+1-j)
     end do
     end do

   else

     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j)=rBuf_S(:,i,j)
     enddo
     enddo

   endif



!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_s >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      endif
      if( itarg_n >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      endif

!
!                           DEALLOCATE sBufferes
!

      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
endsubroutine bocoy_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine bocoy_2d_gh &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nby lines of halos in y direction assuming mirror boundary  !
! conditions at the end of domain. Version for high generations        !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!**********************************************************************!
(this,W,km_in,im_in,jm_in,nbx,nby,Fimax_in,Fjmax_in,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby,mygen_min,mygen_max
real(r_kind),dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:):: sBuf_N,sBuf_S             &
                                             ,rBuf_N,rBuf_S

integer(i_kind) itarg_n,itarg_s,imax,jmax
logical:: lwest,least,lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatay
integer(i_kind) g_ind,g
logical l_sidesend
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!
 
       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2 
         g = my_hgen
         l_sidesend=.true.
       else
         l_sidesend=.false.
       endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 


          if(least) then
            imax = Fimax_in(g)
          else 
            imax = im_in       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax_in(g)
          else  
            jmax = jm_in
          endif


!-----------------------------------------------------------------------
      ndatay = km_in*imax*nby

!
!  SEND boundaries to SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km_in,1:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=1,imax
                    sBuf_S(:,i,j) = W(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km_in,1:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=1,imax
                    sBuf_N(:,i,j)=W(:,i,jmax-nby+j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_work, sHandle(1), isend)

      end if
!
!     RECEIVE boundaries from NORTH and SOUTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km_in,1:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km_in,1:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if

!
! Assign received values from NORTH and SOUTH
!

! --- from NORTH ---

   if( lnorth) then

     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j)=W(:,i,jmax+1-j)
     enddo
     enddo

   else

     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j)=rBuf_N(:,i,j)
     enddo
     enddo

   endif

! From south

   if(lsouth) then

     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j)=W(:,i,nby+1-j)
     end do
     end do

   else

     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j)=rBuf_S(:,i,j)
     enddo
     enddo

   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_s >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      endif
      if( itarg_n >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      endif

!
!                           DEALLOCATE sBufferes
!

      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
endsubroutine bocoy_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine bocoTx_2d_g1 &
!***********************************************************************
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nbx lines close to edges of the subdomins from neighboring  !
! halos in x direction assuming mirror boundary conditions             !
! Version for generation 1                                             !                                             
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(this,W,km_in,im_in,jm_in,nbx,nby)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby
real(r_kind), dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:):: sBuf_E,sBuf_W             &
                                             ,rBuf_E,rBuf_W   

integer(i_kind) itarg_w,itarg_e,imax,jmax
logical lwest,least

integer(i_kind) sHandle(2),rHandle(2),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!
! Limit comminications to selected number of generations
!


         g_ind=1
!
! from mg_domain
!
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)

          imax = im_in    
          jmax = jm_in


!----------------------------------------------------------------------
      ndatax =km_in*jmax*nbx

!
! SEND halos toward WEST and EAST
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km_in,1:nbx,1:jmax), stat = iaerr )

              do j=1,jmax
              do i=1-nbx,0
                sBuf_W(:,i+nbx,j) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_comp, sHandle(1), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km_in,1:nbx,1:jmax), stat = iaerr )

              do j=1,jmax
              do i=1,nbx
                sBuf_E(:,i,j) = W(:,imax+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_comp, sHandle(2), isend)

      end if

!
! RECEIVE halos from EAST and WEST
!

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e


          allocate( rBuf_E(1:km_in,1:nbx,1:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w


         allocate( rBuf_W(1:km_in,1:nbx,1:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )


      end if

!
! Assign received halos from WEST and EAST to interrior of domains
!

! From west

   if(lwest) then
     do j=1,jmax
     do i=1,nbx
       W(:,i,j)= W(:,i,j)+W(:,1-i,j)
     end do
     end do
   else
     do j=1,jmax
     do i=1,nbx
      W(:,i,j)= W(:,i,j)+rBuf_W(:,i,j)
     end do
     end do
   endif

! From east

   if(least) then
     do j=1,jmax
     do i=1,nbx
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+W(:,imax+nbx+1-i,j)
     end do
     end do
   else 
     do j=1,jmax
     do i=1,nbx  
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+rBuf_E(:,i,j)
     end do
     end do
   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_w  >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      endif
      if( itarg_e  >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      endif

!
!                           DEALLOCATE sBufferes
!

      if( itarg_w  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if
      if( itarg_e  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if

!-----------------------------------------------------------------------
endsubroutine bocoTx_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine bocoTx_2d_gh &
!***********************************************************************
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nbx lines close to edges of the subdomins from neighboring  !
! halos in x direction assuming mirror boundary conditions             !
! Version for high generations                                         ! 
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(this,W,km_in,im_in,jm_in,nbx,nby,Fimax_in,Fjmax_in,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby,mygen_min,mygen_max
real(r_kind), dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:):: sBuf_E,sBuf_W             &
                                             ,rBuf_E,rBuf_W   
integer(i_kind) itarg_w,itarg_e,imax,jmax
logical lwest,least,lnorth

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!
! Limit comminications to selected number of generations
!

       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2
         g = my_hgen
         l_sidesend=.true.
       else 
         l_sidesend=.false.
       endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!
! from mg_domain
!
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)

          lnorth  = Flnorth(g_ind)


          if(least) then
            imax = Fimax_in(g)
          else 
            imax = im_in       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax_in(g)
          else  
            jmax = jm_in
          endif


!----------------------------------------------------------------------
      ndatax =km_in*jmax*nbx
!
! SEND halos toward WEST and EAST
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km_in,1:nbx,1:jmax), stat = iaerr )

              do j=1,jmax
              do i=1-nbx,0
                sBuf_W(:,i+nbx,j) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km_in,1:nbx,1:jmax), stat = iaerr )

              do j=1,jmax
              do i=1,nbx
                sBuf_E(:,i,j) = W(:,imax+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(2), isend)

      end if

!
! RECEIVE halos from EAST and WEST
!

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km_in,1:nbx,1:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km_in,1:nbx,1:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if
!
! Assign received values from WEST and EAST
!

! From west

   if(lwest) then
     do j=1,jmax
     do i=1,nbx
       W(:,i,j)= W(:,i,j)+W(:,1-i,j)
     end do
     end do
   else
     do j=1,jmax
     do i=1,nbx
      W(:,i,j)= W(:,i,j)+rBuf_W(:,i,j)
     end do
     end do
   endif

! From east

   if(least) then
     do j=1,jmax
     do i=1,nbx
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+W(:,imax+nbx+1-i,j)
     end do
     end do
   else 
     do j=1,jmax
     do i=1,nbx  
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+rBuf_E(:,i,j)
     end do
     end do
   endif

!-----------------------------------------------------------------------

!                           DEALLOCATE rBufferes

      if( itarg_w  >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      end if
      if( itarg_e  >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      end if

!                           DEALLOCATE sBufferes

      if( itarg_w  >= 0 ) then
         call MPI_WAIT( sHandle(4), istat, ierr )
      end if
      if( itarg_e  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
endsubroutine bocoTx_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine bocoTy_2d_g1 &
!***********************************************************************
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nby lines close to edges of the subdomins from neighboring  !
! halos in y direction assuming mirror boundary conditions             !
! Version for generation 1                                             !                                             
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(this,W,km_in,im_in,jm_in,nbx,nby)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby
real(r_kind), dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:):: sBuf_N,sBuf_S             &
                                             ,rBuf_N,rBuf_S

integer(i_kind) itarg_n,itarg_s,imax,jmax
logical lsouth,lnorth                                       

integer(i_kind) sHandle(2),rHandle(2),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!
! Limit comminications to selected number of generations
!

         g_ind=1
!
! from mg_domain
!
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)

          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)

          imax = im_in    
          jmax = jm_in


!----------------------------------------------------------------------
      ndatay =km_in*imax*nby

!
! SEND SOUTH and NORTH halos
!
! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km_in,1:imax,1:nby), stat = iaerr )

              do j=1-nby,0
              do i=1,imax
                sBuf_S(:,i,j+nby) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_comp, sHandle(1), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km_in,1:imax,1:nby), stat = iaerr )

              do j=1,nby
              do i=1,imax
                sBuf_N(:,i,j)=W(:,i,jmax+j)
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_comp, sHandle(2), isend)

      end if

!
! RECEIVE halos from NORTH and SOUTH
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km_in,1:imax,1:nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km_in,1:imax,1:nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )


      end if

!
! ASSIGN received values from SOUTH and NORTH
!

! From south

   if(lsouth) then
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+W(:,i,1-j)
     end do
     end do
   else
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+rBuf_S(:,i,j)
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+W(:,i,jmax+nby+1-j)
     enddo
     enddo
   else
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+rBuf_N(:,i,j)
     enddo
     enddo
   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_s  >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      end if
      if( itarg_n  >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      end if

!
!                           DEALLOCATE sBufferes
!

      if( itarg_s  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if
      if( itarg_n  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if


!-----------------------------------------------------------------------
endsubroutine bocoTy_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine bocoTy_2d_gh &
!***********************************************************************
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nby lines close to edges of the subdomins from neighboring  !
! halos in y direction assuming mirror boundary conditions             !
! Version for high generations                                         !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(this,W,km_in,im_in,jm_in,nbx,nby,Fimax_in,Fjmax_in,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby,mygen_min,mygen_max
real(r_kind), dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:)::  sBuf_N,sBuf_S            &
                                              ,rBuf_N,rBuf_S
integer(i_kind) itarg_n,itarg_s,itarg_e,imax,jmax
logical least,lsouth,lnorth                                       

integer(i_kind) sHandle(2),rHandle(2),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!


       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2
         g = my_hgen
         l_sidesend=.true.
       else 
         l_sidesend=.false.
       endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!
! from mg_domain
!
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)

          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)


          if(least) then
            imax = Fimax_in(g)
          else 
            imax = im_in       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax_in(g)
          else  
            jmax = jm_in
          endif


!----------------------------------------------------------------------

      ndatay =km_in*imax*nby
!
! SEND halos toward SOUTH and NORTH
!
! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km_in,1:imax,1:nby), stat = iaerr )

              do j=1-nby,0
              do i=1,imax
                sBuf_S(:,i,j+nby) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(1), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km_in,1:imax,1:nby), stat = iaerr )

              do j=1,nby
              do i=1,imax
                sBuf_N(:,i,j)=W(:,i,jmax+j)
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_work, sHandle(2), isend)

      end if

!
! RECEIVE halos from NORTH and SOUTH
!
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km_in,1:imax,1:nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km_in,1:imax,1:nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )


      end if

!
! Assign received values from SOUTH and NORTH
!

! From south

   if(lsouth) then
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+W(:,i,1-j)
     end do
     end do
   else
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+rBuf_S(:,i,j)
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+W(:,i,jmax+nby+1-j)
     enddo
     enddo
   else
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+rBuf_N(:,i,j)
     enddo
     enddo
   endif

!-----------------------------------------------------------------------

!                           DEALLOCATE rBufferes

      if( itarg_s  >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      end if
      if( itarg_n  >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      end if

!                           DEALLOCATE sBufferes

      if( itarg_s  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if
      if( itarg_n  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
endsubroutine bocoTy_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine boco_2d_loc &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies (nbx,nby) lines of halos in (x,y) directions assuming       !
! mirror boundary conditions. Version for localiztion                  !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!**********************************************************************!
(this,W,km_in,im_in,jm_in,nbx,nby,Fimax_in,Fjmax_in,g)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby,g
real(r_kind),dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:)::                           &
                                  sBuf_N,sBuf_E,sBuf_S,sBuf_W           &
                                 ,rBuf_N,rBuf_E,rBuf_S,rBuf_W

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
integer(i_kind) g_ind
logical l_sidesend
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!

         l_sidesend=.true.

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain
!
          itarg_n = Fitarg_n_loc(g)
          itarg_s = Fitarg_s_loc(g)
          itarg_w = Fitarg_w_loc(g)
          itarg_e = Fitarg_e_loc(g)

          lwest   = Flwest_loc(g)
          least   = Fleast_loc(g)
          lsouth  = Flsouth_loc(g)
          lnorth  = Flnorth_loc(g)


!
!  Keep this for now but use only Mod(nxm,8)=Mod(nym,8)=0
!

          if(least) then
            imax = Fimax_in(g)
          else
            imax = im_in       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax_in(g)
          else
            jmax = jm_in
          endif


!-----------------------------------------------------------------------
      ndatay = km_in*imax*nby
      ndatax = km_in*(jmax+2*nby)*nbx


!
!  SEND boundaries to SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km_in,1:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=1,imax
                    sBuf_S(:,i,j) = W(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km_in,1:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=1,imax
                    sBuf_N(:,i,j)=W(:,i,jmax-nby+j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_work, sHandle(1), isend)

      end if
!
!     RECEIVE boundaries from NORTH and SOUTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km_in,1:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km_in,1:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if

!
! Assign received values from NORTH and SOUTH
!


! From south

   if(lsouth) then

     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j)=W(:,i,nby+1-j)
     end do
     end do

   else

     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j)=rBuf_S(:,i,j)
     enddo
     enddo

   endif


! --- from NORTH ---

   if( lnorth) then

     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j)=W(:,i,jmax+1-j)
     enddo
     enddo

   else

     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j)=rBuf_N(:,i,j)
     enddo
     enddo

   endif

!
!  SEND extended boundaries to WEST and EASTH
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km_in,nbx,1-nby:jmax+nby), stat = iaerr )

                do j=1-nby,jmax+nby
                  do i=1,nbx
                    sBuf_W(:,i,j) = W(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km_in,nbx,1-nby:jmax+nby), stat = iaerr )

                do j=1-nby,jmax+nby
                  do i=1,nbx
                    sBuf_E(:,i,j) = W(:,imax-nbx+i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(2), isend)

      end if

!
!     RECEIVE extended boundaries from EAST and WEST
!

! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km_in,nbx,1-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km_in,nbx,1-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if

!
! Assign received values from  WEST and EAST
!

! From west

   if(lwest) then

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,-nbx+i,j)= W(:,nbx+1-i,j)
     end do
     end do

   else

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,-nbx+i,j)= rBuf_W(:,i,j)
     enddo
     enddo


   endif

! From east

   if(least) then

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax+i,j)=W(:,imax+1-i,j)
     end do
     end do

   else

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax+i,j)=rBuf_E(:,i,j)
     enddo
     enddo

   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!
      if( itarg_w >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      end if
      if( itarg_e >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      end if
      if( itarg_s >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      end if
      if( itarg_n >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      end if

!
!                           DEALLOCATE sBufferes
!

      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
endsubroutine boco_2d_loc

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine bocoT_2d_loc &
!***********************************************************************
!                                                                      !
!  Supply n-lines inside of domains, including edges, with halos from  !
!  the surrounding domains.  Assume mirror boundary conditions at the  !
!  boundaries of the domain. Vesrion for localization.                 !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(this,W,km_in,im_in,jm_in,nbx,nby,Fimax_in,Fjmax_in,g)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby,g
real(r_kind), dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:)::                           &
                                        sBuf_N,sBuf_E,sBuf_S,sBuf_W     &
                                       ,rBuf_N,rBuf_E,rBuf_S,rBuf_W   
integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical lwest,least,lsouth,lnorth                                       

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
logical l_sidesend
integer(i_kind) g_ind,k
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!


         g_ind=g
         l_sidesend=.true.


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!
! from mg_domain
!
          itarg_n = Fitarg_n_loc(g_ind)
          itarg_s = Fitarg_s_loc(g_ind)
          itarg_w = Fitarg_w_loc(g_ind)
          itarg_e = Fitarg_e_loc(g_ind)

          lwest   = Flwest_loc(g_ind)
          least   = Fleast_loc(g_ind)
          lsouth  = Flsouth_loc(g_ind)
          lnorth  = Flnorth_loc(g_ind)


          if(least) then
            imax = Fimax_in(g)
          else 
            imax = im_in       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax_in(g)
          else  
            jmax = jm_in
          endif


!----------------------------------------------------------------------
      ndatax =km_in*(jmax+2*nby)*nbx
      ndatay =km_in*imax*nby

!
! SEND extended halos toward WEST and EAST
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km_in,1:nbx,1-nby:jmax+nby), stat = iaerr )

              do j=1-nby,jmax+nby
              do i=1,nbx
                sBuf_W(:,i,j) = W(:,-nbx+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km_in,1:nbx,1-nby:jmax+nby), stat = iaerr )

              do j=1-nby,jmax+nby
              do i=1,nbx
                sBuf_E(:,i,j) = W(:,imax+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(2), isend)

      end if

!
! RECEIVE extended halos from EAST and WEST
!

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km_in,1:nbx,1-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km_in,1:nbx,1-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if
!
! Assign received values from WEST and EAST
!

! From west

   if(lwest) then
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,i,j)= W(:,i,j)+W(:,1-i,j)
     end do
     end do
   else
     do j=1-nby,jmax+nby
     do i=1,nbx
      W(:,i,j)= W(:,i,j)+rBuf_W(:,i,j)
     end do
     end do
   endif

! From east

   if(least) then
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+W(:,imax+1+nbx-i,j)
     end do
     end do
   else 
     do j=1-nby,jmax+nby
     do i=1,nbx  
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+rBuf_E(:,i,j)
     end do
     end do
   endif

!
! SEND halos toward SOUTH and NORTH
!
! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km_in,1:imax,1:nby), stat = iaerr )

              do j=1,nby  
              do i=1,imax
                sBuf_S(:,i,j) = W(:,i,-nby+j)
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km_in,1:imax,1:nby), stat = iaerr )

              do j=1,nby
              do i=1,imax
                sBuf_N(:,i,j)=W(:,i,jmax+j)
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_work, sHandle(1), isend)

      end if

!
! RECEIVE halos from NORTH and SOUTH
!
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km_in,1:imax,1:nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km_in,1:imax,1:nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )


      end if

!
! Assign received values from SOUTH and NORTH
!

! From south

   if(lsouth) then
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+W(:,i,1-j)
     end do
     end do
   else
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+rBuf_S(:,i,j)
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+W(:,i,jmax+1+nby-j)
     enddo
     enddo
   else
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+rBuf_N(:,i,j)
     enddo
     enddo
   endif

!-----------------------------------------------------------------------

!                           DEALLOCATE rBufferes

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)
        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!                           DEALLOCATE sBufferes

      if( itarg_w  >= 0 ) then
         call MPI_WAIT( sHandle(4), istat, ierr )
      end if
      if( itarg_e  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if
      if( itarg_s  >= 0 ) then
         call MPI_WAIT( sHandle(3), istat, ierr )
      end if
      if( itarg_n  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
endsubroutine bocoT_2d_loc

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine upsend_loc_g12 &
!***********************************************************************
!                                                                      !
!         Upsend data from generation one to generation two            !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(this,V_in,H,km_4_in,flag)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_4_in,flag
real(r_kind), dimension(km_4_in,1:this%imL,1:this%jmL),intent(in):: V_in
real(r_kind), dimension(km_4_in,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(out):: H
!-----------------------------------------------------------------------
real(r_kind),dimension(1:km_4_in,1:this%imL,1:this%jmL):: dBuf_SW
real(r_kind),dimension(1:km_4_in,1:this%imL,1:this%jmL):: dBuf_SE
real(r_kind),dimension(1:km_4_in,1:this%imL,1:this%jmL):: dBuf_NW
real(r_kind),dimension(1:km_4_in,1:this%imL,1:this%jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,ndata,i,j
integer(i_kind):: mygen_dn,mygen_up
integer(i_kind):: itarg_up
logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne
integer(i_kind), parameter :: DIR_SW=1, DIR_SE=2, DIR_NW=3, DIR_NE=4
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------
   mygen_dn=1
   mygen_up=2
!
! Define generational flags
!

       itarg_up=Fitargup_loc12(flag)

       lsendup_sw = lsendup_sw_loc
       lsendup_se = lsendup_se_loc
       lsendup_nw = lsendup_nw_loc
       lsendup_ne = lsendup_ne_loc
!-----------------------------------------------------------------------

      H(:,:,:) = 0.0d0

     ndata =km_4_in*imL*jmL

     sHandle(:) = MPI_REQUEST_NULL
     rHandle(:) = MPI_REQUEST_NULL

      if(  lsendup_sw ) then
        if(itarg_up == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_SW(:,i,j) = V_in(:,i,j)
             enddo
             enddo
        else
        call MPI_ISEND( V_in, ndata, dtype, itarg_up, mype,  &
                       mpi_comm_comp, sHandle(DIR_SW), ierr)
        endif
      endif
!
      if( itargdn_sw_loc21 >= 0 ) then
        if(itargdn_sw_loc21 /= mype) then
          call MPI_IRECV( dBuf_SW, ndata, dtype, itargdn_sw_loc21, itargdn_sw_loc21,          &
                          mpi_comm_comp, rHandle(DIR_SW), ierr)
        endif
      endif

      if( lsendup_se ) then
        if(itarg_up == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_SE(:,i,j) = V_in(:,i,j)
             enddo
             enddo
        else
        call MPI_ISEND( V_in, ndata, dtype, itarg_up, mype, &
                       mpi_comm_comp, sHandle(DIR_SE), ierr)
        endif
      end if

      if( itargdn_se_loc21 >= 0 ) then
        if(itargdn_se_loc21 /= mype) then
          call MPI_IRECV( dBuf_SE, ndata, dtype, itargdn_se_loc21, itargdn_se_loc21,          &
                          mpi_comm_comp, rHandle(DIR_SE), ierr)
        endif
      end if

      if( lsendup_nw ) then
        if(itarg_up == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_NW(:,i,j) = V_in(:,i,j)
             enddo
             enddo
        else
         call MPI_ISEND( V_in, ndata, dtype, itarg_up, mype,  &
                        mpi_comm_comp, sHandle(DIR_NW), ierr)
        end if
      end if

      if( itargdn_nw_loc21 >= 0 ) then
        if(itargdn_nw_loc21 /= mype) then
          call MPI_IRECV( dBuf_NW, ndata, dtype, itargdn_nw_loc21, itargdn_nw_loc21,          &
                          mpi_comm_comp, rHandle(DIR_NW), ierr)
        endif
      end if

      if( lsendup_ne ) then
        if(itarg_up == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_NE(:,i,j) = V_in(:,i,j)
             enddo
             enddo
        else
        call MPI_ISEND( V_in, ndata, dtype, itarg_up, mype, &
                      mpi_comm_comp, sHandle(DIR_NE), ierr)
        endif
      end if

      if( itargdn_ne_loc21 >= 0 ) then
        if(itargdn_ne_loc21 /= mype) then
          call MPI_IRECV( dBuf_NE, ndata, dtype, itargdn_ne_loc21, itargdn_ne_loc21,          &
                         mpi_comm_comp, rHandle(DIR_NE), ierr)
        endif
      end if

      if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
      end if
      if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
      end if

      if( itargdn_sw_loc21 >= 0 ) then
             do j=1,jmL
             do i=1,imL
                H(:,i,j)=dBuf_SW(:,i,j)
             enddo
             enddo
      endif

      if( itargdn_se_loc21 >= 0 ) then
             do j=1,jmL
             do i=1,imL
               H(:,imL+i,j)=dBuf_SE(:,i,j)
             enddo
             enddo
      endif

      if( itargdn_nw_loc21 >= 0 ) then
             do j=1,jmL
             do i=1,imL
               H(:,i,jmL+j)=dBuf_NW(:,i,j)
             enddo
             enddo
      endif

      if( itargdn_ne_loc21 >= 0 ) then
             do j=1,jmL
             do i=1,imL
               H(:,imL+i,jmL+j)=dBuf_NE(:,i,j)
             enddo
             enddo
      endif

!-----------------------------------------------------------------------
endsubroutine upsend_loc_g12

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine upsend_loc_g23 &
!***********************************************************************
!                                                                      !
!         Upsend data from generation three to generation four         !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(this,V_in,H,km_16_in,flag)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_16_in,flag
real(r_kind), dimension(km_16_in,1:this%imL,1:this%jmL),intent(in):: V_in
real(r_kind), dimension(km_16_in,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(out):: H
!-----------------------------------------------------------------------
real(r_kind),dimension(1:km_16_in,1:this%imL,1:this%jmL):: dBuf_SW
real(r_kind),dimension(1:km_16_in,1:this%imL,1:this%jmL):: dBuf_SE
real(r_kind),dimension(1:km_16_in,1:this%imL,1:this%jmL):: dBuf_NW
real(r_kind),dimension(1:km_16_in,1:this%imL,1:this%jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,ndata,i,j
integer(i_kind):: mygen_dn,mygen_up
integer(i_kind):: itarg_up
logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne
integer(i_kind), parameter :: DIR_SW=1, DIR_SE=2, DIR_NW=3, DIR_NE=4
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------
   mygen_dn=2
   mygen_up=3
!
! Define generational flags
!

       itarg_up=Fitargup_loc23(flag)

       lsendup_sw = lsendup_sw_loc
       lsendup_se = lsendup_se_loc
       lsendup_nw = lsendup_nw_loc
       lsendup_ne = lsendup_ne_loc
!-----------------------------------------------------------------------

      H(:,:,:) = 0.0d0

     ndata =km_16_in*imL*jmL

     sHandle(:) = MPI_REQUEST_NULL
     rHandle(:) = MPI_REQUEST_NULL

      if(  lsendup_sw ) then
        if(itarg_up == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_SW(:,i,j) = V_in(:,i,j)
             enddo
             enddo
        else
        call MPI_ISEND( V_in, ndata, dtype, itarg_up, mype,  &
                       mpi_comm_comp, sHandle(DIR_SW), ierr)
        endif
      endif

      if( itargdn_sw_loc32 >= 0 ) then
        if(itargdn_sw_loc32 /= mype) then
          call MPI_IRECV( dBuf_SW, ndata, dtype, itargdn_sw_loc32, itargdn_sw_loc32,          &
                          mpi_comm_comp, rHandle(DIR_SW), ierr)
        endif
      endif

      if( lsendup_se ) then
        if(itarg_up == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_SE(:,i,j) = V_in(:,i,j)
             enddo
             enddo
        else
        call MPI_ISEND( V_in, ndata, dtype, itarg_up, mype, &
                       mpi_comm_comp, sHandle(DIR_SE), ierr)
        endif
      end if

      if( itargdn_se_loc32 >= 0 ) then
        if(itargdn_se_loc32 /= mype) then
          call MPI_IRECV( dBuf_SE, ndata, dtype, itargdn_se_loc32, itargdn_se_loc32,          &
                          mpi_comm_comp, rHandle(DIR_SE), ierr)
        endif
      end if

      if( lsendup_nw ) then
        if(itarg_up == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_NW(:,i,j) = V_in(:,i,j)
             enddo
             enddo
        else
         call MPI_ISEND( V_in, ndata, dtype, itarg_up, mype,  &
                        mpi_comm_comp, sHandle(DIR_NW), ierr)
        end if
      end if

      if( itargdn_nw_loc32 >= 0 ) then
        if(itargdn_nw_loc32 /= mype) then
          call MPI_IRECV( dBuf_NW, ndata, dtype, itargdn_nw_loc32, itargdn_nw_loc32,          &
                          mpi_comm_comp, rHandle(DIR_NW), ierr)
        endif
      end if

      if( lsendup_ne ) then
        if(itarg_up == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_NE(:,i,j) = V_in(:,i,j)
             enddo
             enddo
        else
        call MPI_ISEND( V_in, ndata, dtype, itarg_up, mype, &
                      mpi_comm_comp, sHandle(DIR_NE), ierr)
        endif
      end if

      if( itargdn_ne_loc32 >= 0 ) then
        if(itargdn_ne_loc32 /= mype) then
          call MPI_IRECV( dBuf_NE, ndata, dtype, itargdn_ne_loc32, itargdn_ne_loc32,          &
                         mpi_comm_comp, rHandle(DIR_NE), ierr)
        endif
      end if

      if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
      end if
      if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
      end if

      if( itargdn_sw_loc32 >= 0 ) then
             do j=1,jmL
             do i=1,imL
                H(:,i,j)=dBuf_SW(:,i,j)
             enddo
             enddo
      endif

      if( itargdn_se_loc32 >= 0 ) then
             do j=1,jmL
             do i=1,imL
               H(:,imL+i,j)=dBuf_SE(:,i,j)
             enddo
             enddo
      endif

      if( itargdn_nw_loc32 >= 0 ) then
             do j=1,jmL
             do i=1,imL
               H(:,i,jmL+j)=dBuf_NW(:,i,j)
             enddo
             enddo
      endif

      if( itargdn_ne_loc32 >= 0 ) then
             do j=1,jmL
             do i=1,imL
               H(:,imL+i,jmL+j)=dBuf_NE(:,i,j)
             enddo
             enddo
      endif

!-----------------------------------------------------------------------
endsubroutine upsend_loc_g23

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine upsend_loc_g34 &
!***********************************************************************
!                                                                      !
!         Upsend data from generation three to generation four         !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(this,V_in,H,km_64_in,flag)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_64_in,flag
real(r_kind), dimension(km_64_in,1:this%imL,1:this%jmL),intent(in):: V_in
real(r_kind), dimension(km_64_in,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(out):: H
!-----------------------------------------------------------------------
real(r_kind),dimension(1:km_64_in,1:this%imL,1:this%jmL):: dBuf_SW
real(r_kind),dimension(1:km_64_in,1:this%imL,1:this%jmL):: dBuf_SE
real(r_kind),dimension(1:km_64_in,1:this%imL,1:this%jmL):: dBuf_NW
real(r_kind),dimension(1:km_64_in,1:this%imL,1:this%jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,ndata,i,j
integer(i_kind):: mygen_dn,mygen_up
integer(i_kind):: itarg_up
logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne
integer(i_kind), parameter :: DIR_SW=1, DIR_SE=2, DIR_NW=3, DIR_NE=4
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------
   mygen_dn=3
   mygen_up=4
!
! Define generational flags
!

       itarg_up=Fitargup_loc34(flag)

       lsendup_sw = lsendup_sw_loc
       lsendup_se = lsendup_se_loc
       lsendup_nw = lsendup_nw_loc
       lsendup_ne = lsendup_ne_loc
!-----------------------------------------------------------------------

      H(:,:,:) = 0.0d0

     ndata =km_64_in*imL*jmL

     sHandle(:) = MPI_REQUEST_NULL
     rHandle(:) = MPI_REQUEST_NULL

      if(  lsendup_sw ) then
        if(itarg_up == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_SW(:,i,j) = V_in(:,i,j)
             enddo
             enddo
        else
        call MPI_ISEND( V_in, ndata, dtype, itarg_up, mype,  &
                       mpi_comm_comp, sHandle(DIR_SW), ierr)
        endif
      endif

      if( itargdn_sw_loc43 >= 0 ) then
        if(itargdn_sw_loc43 /= mype) then
          call MPI_IRECV( dBuf_SW, ndata, dtype, itargdn_sw_loc43, itargdn_sw_loc43,          &
                          mpi_comm_comp, rHandle(DIR_SW), ierr)
        endif
      endif

      if( lsendup_se ) then
        if(itarg_up == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_SE(:,i,j) = V_in(:,i,j)
             enddo
             enddo
        else
        call MPI_ISEND( V_in, ndata, dtype, itarg_up, mype, &
                       mpi_comm_comp, sHandle(DIR_SE), ierr)
        endif
      end if

      if( itargdn_se_loc43 >= 0 ) then
        if(itargdn_se_loc43 /= mype) then
          call MPI_IRECV( dBuf_SE, ndata, dtype, itargdn_se_loc43, itargdn_se_loc43,          &
                          mpi_comm_comp, rHandle(DIR_SE), ierr)
        endif
      end if

      if( lsendup_nw ) then
        if(itarg_up == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_NW(:,i,j) = V_in(:,i,j)
             enddo
             enddo
        else
         call MPI_ISEND( V_in, ndata, dtype, itarg_up, mype,  &
                        mpi_comm_comp, sHandle(DIR_NW), ierr)
        end if
      end if

      if( itargdn_nw_loc43 >= 0 ) then
        if(itargdn_nw_loc43 /= mype) then
          call MPI_IRECV( dBuf_NW, ndata, dtype, itargdn_nw_loc43, itargdn_nw_loc43,          &
                          mpi_comm_comp, rHandle(DIR_NW), ierr)
        endif
      end if

      if( lsendup_ne ) then
        if(itarg_up == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_NE(:,i,j) = V_in(:,i,j)
             enddo
             enddo
        else
        call MPI_ISEND( V_in, ndata, dtype, itarg_up, mype, &
                      mpi_comm_comp, sHandle(DIR_NE), ierr)
        endif
      end if

      if( itargdn_ne_loc43 >= 0 ) then
        if(itargdn_ne_loc43 /= mype) then
          call MPI_IRECV( dBuf_NE, ndata, dtype, itargdn_ne_loc43, itargdn_ne_loc43,          &
                         mpi_comm_comp, rHandle(DIR_NE), ierr)
        endif
      end if

      if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
      end if
      if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
      end if

      if( itargdn_sw_loc43 >= 0 ) then
             do j=1,jmL
             do i=1,imL
                H(:,i,j)=dBuf_SW(:,i,j)
             enddo
             enddo
      endif

      if( itargdn_se_loc43 >= 0 ) then
             do j=1,jmL
             do i=1,imL
               H(:,imL+i,j)=dBuf_SE(:,i,j)
             enddo
             enddo
      endif

      if( itargdn_nw_loc43 >= 0 ) then
             do j=1,jmL
             do i=1,imL
               H(:,i,jmL+j)=dBuf_NW(:,i,j)
             enddo
             enddo
      endif

      if( itargdn_ne_loc43 >= 0 ) then
             do j=1,jmL
             do i=1,imL
               H(:,imL+i,jmL+j)=dBuf_NE(:,i,j)
             enddo
             enddo
      endif

!-----------------------------------------------------------------------
endsubroutine upsend_loc_g34

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine downsend_loc_g43 &
!***********************************************************************
!                                                                      *
!         Downsending data from low resolution pes    (mygen_up)       *
!         to the concurent high-resolution pes        (mygen_dn)       *
!         and add the existing and the recevied values                 *
!                                                                      !
!                       - offset version -                             !
!                                                                      *
!***********************************************************************
(this,W,Z,km_64_in,flag)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_64_in,flag
real(r_kind), dimension(km_64_in,1:this%im,1:this%jm),intent(in):: W
real(r_kind), dimension(km_64_in,1:this%imL,1:this%jmL),intent(out):: Z
!-----------------------------------------------------------------------
real(r_kind),dimension(1:km_64_in,1:this%imL,1:this%jmL):: dBuf_SW
real(r_kind),dimension(1:km_64_in,1:this%imL,1:this%jmL):: dBuf_SE
real(r_kind),dimension(1:km_64_in,1:this%imL,1:this%jmL):: dBuf_NW
real(r_kind),dimension(1:km_64_in,1:this%imL,1:this%jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,ndata,i,j
integer(i_kind):: itarg_up                                           
logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne
integer(i_kind), parameter :: DIR_SW=1, DIR_SE=2, DIR_NW=3, DIR_NE=4
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------

     Z(:,:,:) = 0.0d0
!
! Define generational flags
!

       itarg_up=Fitargup_loc34(flag)

       lsendup_sw = lsendup_sw_loc
       lsendup_se = lsendup_se_loc
       lsendup_nw = lsendup_nw_loc
       lsendup_ne = lsendup_ne_loc

       ndata =km_64_in*imL*jmL

      dBuf_SW = 0.0d0
      dBuf_SE = 0.0d0
      dBuf_NW = 0.0d0
      dBuf_NE = 0.0d0

      sHandle(:) = MPI_REQUEST_NULL
      rHandle(:) = MPI_REQUEST_NULL

     if(itargdn_sw_loc43 >= 0) then
        if(itargdn_sw_loc43 == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_SW(:,i,j) = W(:,i,j)
             enddo
             enddo
        else
        call MPI_ISEND( W(1,1,1), ndata, dtype, itargdn_sw_loc43, mype,  &
                        mpi_comm_work, sHandle(DIR_SW), ierr )
        end if
     endif

      if( lsendup_sw .and. itarg_up /= mype ) then
        call MPI_IRECV( dBuf_SW, ndata, dtype, itarg_up, itarg_up, &
                        mpi_comm_work, rHandle(DIR_SW), ierr)
      endif

     if(itargdn_se_loc43 >= 0) then
        if(itargdn_se_loc43 == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_SE(:,i,j) = W(:,imL+i,j)
             enddo
             enddo
        else
        call MPI_ISEND( W(1,imL+1,1), ndata, dtype, itargdn_se_loc43, mype,  &
                       mpi_comm_work, sHandle(DIR_SE), ierr )
        end if
     endif

      if( lsendup_se .and. itarg_up /= mype ) then
        call MPI_IRECV( dBuf_SE, ndata, dtype, itarg_up, itarg_up, &
                        mpi_comm_work, rHandle(DIR_SE), ierr)
      end if

     if(itargdn_nw_loc43 >= 0) then
        if(itargdn_nw_loc43 == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_NW(:,i,j) = W(:,i,jmL+j)
             enddo
             enddo
        else
        call MPI_ISEND( W(1,1,jmL+1), ndata, dtype, itargdn_nw_loc43, mype,  &
                        mpi_comm_work, sHandle(DIR_NW), ierr )
        end if
     endif

      if( lsendup_nw .and. itarg_up /= mype ) then
        call MPI_IRECV( dBuf_NW, ndata, dtype, itarg_up, itarg_up, &
                       mpi_comm_work, rHandle(DIR_NW), ierr)
      end if

     if(itargdn_ne_loc43 >= 0) then
        if(itargdn_ne_loc43 == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_NE(:,i,j) = W(:,imL+i,jmL+j)
             enddo
             enddo
        else
        call MPI_ISEND( W(1,imL+1,jmL+1), ndata, dtype, itargdn_ne_loc43, mype,  &
                        mpi_comm_work, sHandle(DIR_NE), ierr )
        end if
     endif

      if( lsendup_ne .and. itarg_up /= mype ) then
        call MPI_IRECV( dBuf_NE, ndata, dtype, itarg_up, itarg_up, &
                        mpi_comm_work, rHandle(DIR_NE), ierr)
      end if

      if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
      end if
      if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
      end if

      if( lsendup_sw ) then
             do j=1,jmL
             do i=1,imL
               Z(:,i,j)=dBuf_SW(:,i,j)  
             enddo
             enddo
      endif

      if( lsendup_se ) then
             do j=1,jmL
             do i=1,imL
               Z(:,i,j)=dBuf_SE(:,i,j)
             enddo
             enddo
      endif

      if( lsendup_nw ) then
             do j=1,jmL
             do i=1,imL
               Z(:,i,j)=dBuf_NW(:,i,j)
             enddo
             enddo
      endif

      if( lsendup_ne ) then
             do j=1,jmL
             do i=1,imL
               Z(:,i,j)=dBuf_NE(:,i,j)
             enddo
             enddo
      endif

!-----------------------------------------------------------------------
endsubroutine downsend_loc_g43

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine downsend_loc_g32 &
!***********************************************************************
!                                                                      *
!         Downsending data from low resolution pes    (mygen_up)       *
!         to the concurent high-resolution pes        (mygen_dn)       *
!         and add the existing and the recevied values                 *
!                                                                      !
!                       - offset version -                             !
!                                                                      *
!***********************************************************************
(this,Z,H,km_16_in,flag)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_16_in,flag
real(r_kind), dimension(km_16_in,1:this%im,1:this%jm),intent(in):: Z
real(r_kind), dimension(km_16_in,1:this%imL,1:this%jmL),intent(out):: H
!-----------------------------------------------------------------------
real(r_kind),dimension(1:km_16_in,1:this%imL,1:this%jmL):: dBuf_SW
real(r_kind),dimension(1:km_16_in,1:this%imL,1:this%jmL):: dBuf_SE
real(r_kind),dimension(1:km_16_in,1:this%imL,1:this%jmL):: dBuf_NW
real(r_kind),dimension(1:km_16_in,1:this%imL,1:this%jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,ndata,i,j
integer(i_kind):: itarg_up                                           
logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne
integer(i_kind), parameter :: DIR_SW=1, DIR_SE=2, DIR_NW=3, DIR_NE=4
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------

     H(:,:,:) = 0.0d0
!
! Define generational flags
!

       itarg_up=Fitargup_loc23(flag)

       lsendup_sw = lsendup_sw_loc
       lsendup_se = lsendup_se_loc
       lsendup_nw = lsendup_nw_loc
       lsendup_ne = lsendup_ne_loc

       ndata =km_16_in*imL*jmL

      dBuf_SW = 0.0d0
      dBuf_SE = 0.0d0
      dBuf_NW = 0.0d0
      dBuf_NE = 0.0d0

      sHandle(:) = MPI_REQUEST_NULL
      rHandle(:) = MPI_REQUEST_NULL

  if( itargdn_sw_loc32 >= 0 ) then
        if(itargdn_sw_loc32 == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_SW(:,i,j) = Z(:,i,j)
             enddo
             enddo
        else
        call MPI_ISEND( Z(1,1,1), ndata, dtype, itargdn_sw_loc32, mype,  &
                        mpi_comm_work, sHandle(DIR_SW), ierr )
        end if
  endif

  if( lsendup_sw .and. itarg_up /= mype ) then
        call MPI_IRECV( dBuf_SW, ndata, dtype, itarg_up, itarg_up, &
                        mpi_comm_work, rHandle(DIR_SW), ierr)
  endif

  if( itargdn_se_loc32 >= 0 ) then
        if(itargdn_se_loc32 == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_SE(:,i,j) = Z(:,imL+i,j)
             enddo
             enddo
        else
        call MPI_ISEND( Z(1,imL+1,1), ndata, dtype, itargdn_se_loc32, mype,  &
                       mpi_comm_work, sHandle(DIR_SE), ierr )
        end if
  endif

  if( lsendup_se .and. itarg_up /= mype ) then
        call MPI_IRECV( dBuf_SE, ndata, dtype, itarg_up, itarg_up, &
                        mpi_comm_work, rHandle(DIR_SE), ierr)
  end if

  if( itargdn_nw_loc32 >= 0 ) then
        if(itargdn_nw_loc32 == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_NW(:,i,j) = Z(:,i,jmL+j)
             enddo
             enddo
        else
        call MPI_ISEND( Z(1,1,jmL+1), ndata, dtype, itargdn_nw_loc32, mype,  &
                        mpi_comm_work, sHandle(DIR_NW), ierr )
        end if
  endif

  if( lsendup_nw .and. itarg_up /= mype ) then
        call MPI_IRECV( dBuf_NW, ndata, dtype, itarg_up, itarg_up, &
                        mpi_comm_work, rHandle(DIR_NW), ierr)
  end if

  if( itargdn_ne_loc32 >= 0 ) then
        if(itargdn_ne_loc32 == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_NE(:,i,j) = Z(:,imL+i,jmL+j)
             enddo
             enddo
        else
        call MPI_ISEND( Z(1,imL+1,jmL+1), ndata, dtype, itargdn_ne_loc32, mype,  &
                        mpi_comm_work, sHandle(DIR_NE), ierr )
        end if
  endif

  if( lsendup_ne .and. itarg_up /= mype ) then
        call MPI_IRECV( dBuf_NE, ndata, dtype, itarg_up, itarg_up, &
                        mpi_comm_work, rHandle(DIR_NE), ierr)
  end if

  if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
  end if
  if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
  end if

  if( lsendup_sw ) then
             do j=1,jmL
             do i=1,imL
               H(:,i,j)=dBuf_SW(:,i,j)  
             enddo
             enddo
  endif

  if( lsendup_se ) then
             do j=1,jmL
             do i=1,imL
               H(:,i,j)=dBuf_SE(:,i,j)
             enddo
             enddo
  endif

  if( lsendup_nw ) then
             do j=1,jmL
             do i=1,imL
               H(:,i,j)=dBuf_NW(:,i,j)
             enddo
             enddo
  endif

  if( lsendup_ne ) then
             do j=1,jmL
             do i=1,imL
               H(:,i,j)=dBuf_NE(:,i,j)
             enddo
             enddo
  endif

!-----------------------------------------------------------------------
endsubroutine downsend_loc_g32

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine downsend_loc_g21 &
!***********************************************************************
!                                                                      *
!         Downsending data from low resolution pes    (mygen_up)       *
!         to the concurent high-resolution pes        (mygen_dn)       *
!         and add the existing and the recevied values                 *
!                                                                      !
!                       - offset version -                             !
!                                                                      *
!***********************************************************************
(this,H,V_out,km_4_in,flag)
!-----------------------------------------------------------------------
use mpi
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_4_in,flag
real(r_kind), dimension(km_4_in,1:this%im,1:this%jm),intent(in):: H
real(r_kind), dimension(km_4_in,1:this%imL,1:this%jmL),intent(out):: V_out
!-----------------------------------------------------------------------
real(r_kind),dimension(1:km_4_in,1:this%imL,1:this%jmL):: dBuf_SW
real(r_kind),dimension(1:km_4_in,1:this%imL,1:this%jmL):: dBuf_SE
real(r_kind),dimension(1:km_4_in,1:this%imL,1:this%jmL):: dBuf_NW
real(r_kind),dimension(1:km_4_in,1:this%imL,1:this%jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4)
integer(i_kind) ierr,ndata,i,j
integer(i_kind):: itarg_up                                           
logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne
integer(i_kind), parameter :: DIR_SW=1, DIR_SE=2, DIR_NW=3, DIR_NE=4
include "type_parameter_locpointer.inc"
include "type_intstat_locpointer.inc"
include "type_parameter_point2this.inc"
include "type_intstat_point2this.inc"
!-----------------------------------------------------------------------

     V_out(:,:,:) = 0.0d0
!
! Define generational flags
!

       itarg_up=Fitargup_loc12(flag)

       lsendup_sw = lsendup_sw_loc
       lsendup_se = lsendup_se_loc
       lsendup_nw = lsendup_nw_loc
       lsendup_ne = lsendup_ne_loc

       ndata =km_4_in*imL*jmL

      dBuf_SW = 0.0d0
      dBuf_SE = 0.0d0
      dBuf_NW = 0.0d0
      dBuf_NE = 0.0d0

      sHandle(:) = MPI_REQUEST_NULL
      rHandle(:) = MPI_REQUEST_NULL

  if( itargdn_sw_loc21 >= 0 ) then
        if(itargdn_sw_loc21 == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_SW(:,i,j) = H(:,i,j)
             enddo
             enddo
        else
        call MPI_ISEND( H(1,1,1), ndata, dtype, itargdn_sw_loc21, mype,  &
                        mpi_comm_work, sHandle(DIR_SW), ierr )
        end if
  endif

      if( lsendup_sw .and. itarg_up /= mype ) then
        call MPI_IRECV( dBuf_SW, ndata, dtype, itarg_up, itarg_up, &
                        mpi_comm_work, rHandle(DIR_SW), ierr)
      endif

  if( itargdn_se_loc21 >= 0 ) then
        if(itargdn_se_loc21 == mype) then
             do j=1,jmL
             do i=1,imL
               dBuf_SE(:,i,j) = H(:,imL+i,j)
             enddo
             enddo
        else
        call MPI_ISEND( H(1,imL+1,1), ndata, dtype, itargdn_se_loc21, mype,  &
                       mpi_comm_work, sHandle(DIR_SE), ierr )
        end if
  endif

      if( lsendup_se .and. itarg_up /= mype ) then
        call MPI_IRECV( dBuf_SE, ndata, dtype, itarg_up, itarg_up, &
                        mpi_comm_work, rHandle(DIR_SE), ierr)
      end if

  if( itargdn_nw_loc21 >= 0 ) then
        if(itargdn_nw_loc21 == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_NW(:,i,j) = H(:,i,jmL+j)
             enddo
             enddo
        else
        call MPI_ISEND( H(1,1,jmL+1), ndata, dtype, itargdn_nw_loc21, mype,  &
                        mpi_comm_work, sHandle(DIR_NW), ierr )
        end if
  endif

      if( lsendup_nw .and. itarg_up /= mype ) then
        call MPI_IRECV( dBuf_NW, ndata, dtype, itarg_up, itarg_up, &
                       mpi_comm_work, rHandle(DIR_NW), ierr)
      end if

  if( itargdn_ne_loc21 >= 0 ) then
        if(itargdn_ne_loc21 == mype) then
             do j=1,jmL
             do i=1,imL
                dBuf_NE(:,i,j) = H(:,imL+i,jmL+j)
             enddo
             enddo
        else
        call MPI_ISEND( H(1,imL+1,jmL+1), ndata, dtype, itargdn_ne_loc21, mype,  &
                        mpi_comm_work, sHandle(DIR_NE), ierr )
        end if
  endif

      if( lsendup_ne .and. itarg_up /= mype ) then
        call MPI_IRECV( dBuf_NE, ndata, dtype, itarg_up, itarg_up, &
                        mpi_comm_work, rHandle(DIR_NE), ierr)
      end if

      if( any(rHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, rHandle, MPI_STATUSES_IGNORE, ierr)
      end if
      if( any(sHandle /= MPI_REQUEST_NULL) ) then
        call MPI_WAITALL(4, sHandle, MPI_STATUSES_IGNORE, ierr)
      end if

      if( lsendup_sw ) then
             do j=1,jmL
             do i=1,imL
               V_out(:,i,j)=dBuf_SW(:,i,j)  
             enddo
             enddo
      endif

      if( lsendup_se ) then
             do j=1,jmL
             do i=1,imL
               V_out(:,i,j)=dBuf_SE(:,i,j)
             enddo
             enddo
      endif

      if( lsendup_nw ) then
             do j=1,jmL
             do i=1,imL
               V_out(:,i,j)=dBuf_NW(:,i,j)
             enddo
             enddo
      endif

      if( lsendup_ne ) then
             do j=1,jmL
             do i=1,imL
               V_out(:,i,j)=dBuf_NE(:,i,j)
             enddo
             enddo
      endif

!-----------------------------------------------------------------------
endsubroutine downsend_loc_g21

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end submodule mg_bocos
































