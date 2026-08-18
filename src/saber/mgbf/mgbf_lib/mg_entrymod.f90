submodule(mg_intstate) mg_entrymod
!$$$  submodule documentation block
!                .      .    .                                       .
! module:   mg_entrymod
!   prgmmr: rancic           org: NCEP/EMC            date: 2020
!
! abstract:  Initialize and finialize multigrid Beta filter
!            for modeling of background error covariance
!
! module history log:
!   2023-04-19  lei     - object-oriented coding
!   2024-01-11  rancic  - optimization for ensemble localization
!   2024-02-20  yokota  - refactoring to apply for GSI
!
! Subroutines Included:
!   mg_initialize -
!   mg_finalize -
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

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine mg_initialize(this,n_owned_anl,anl_lonlat1d,inputfilename,obj_parameter,mpi_comm)
implicit none
!**********************************************************************!
!                                                                      !
!   Initialization subroutine                                          !
!                                                     M. Rancic (2020) !
!***********************************************************************
class (mg_intstate_type), intent(inout) :: this
integer(i_kind),optional,intent(in):: n_owned_anl
real(r_kind),optional,intent(in):: anl_lonlat1d(:,:)
character(len=*),optional,intent(in) :: inputfilename

class(mg_parameter_type),optional,intent(in):: obj_parameter
integer(i_kind),intent(in) :: mpi_comm
integer(i_kind) :: owned_comm
integer(i_kind) :: ierr

!---------------------------------------------------------------------------
!
!               Firs set of subroutines is called only once and serves to
!               initialte the MGBF run
!
!---------------------------------------------------------------------------

!****
!**** Initialize run multigrid Beta filter parameters
!****
call MPI_Comm_dup(mpi_comm, owned_comm, ierr)
this%mpi_comm_comp=owned_comm
if (present(inputfilename)) then
   call this%init_mg_parameter(inputfilename)
else if (present(obj_parameter)) then
   this%mg_parameter_type=obj_parameter
   ! The derived-type assignment also copies the source communicator handle.
   ! Restore the private communicator owned by this MGBF instance.
   this%mpi_comm_comp=owned_comm
end if

 if (present(anl_lonlat1d)) then
    if (size(anl_lonlat1d,2) /= 2 .or. size(anl_lonlat1d,1) <  n_owned_anl) then
      write(6,*)"thinkdeb size(anl_lonlat1d,2) ",size(anl_lonlat1d,2)
      write(6,*)"thinkdeb size(anl_lonlat1d,1) ",size(anl_lonlat1d,1)
      write(6,*)"thinkdeb n_owned_anl ) ", n_owned_anl
      call flush(6)
      error stop "anl_lonlat1d has wrong shape"
    end if

 end if

!****
!**** Initialize MPI
!****
if(this%nxm*this%nym>1) call this%init_mg_MPI

!***
!*** Initialize integration domain
!***
      write(6,*)"thinkdeb in mg_entry,  ", 3
      call flush(6)
call this%init_mg_domain
if(this%l_loc) then
   call this%init_domain_loc
end if

      write(6,*)"thinkdeb in mg_entry,  ", 4
      call flush(6)
!---------------------------------------------------------------------------
!
!               All others are function of km2,km3,km,nm,mm,im,jm
!               and needs to be called separately for each application
!
!---------------------------------------------------------------------------
!***
!*** Define km and WORKA array based on input from mg_parameters and
!*** depending on specific application
!***

!***
!*** Allocate variables, define weights, prepare mapping
!*** between analysis and filter grid
!***

call this%allocate_mg_intstate
      write(6,*)"thinkdeb in mg_entry,  ", 5
      call flush(6)

call this%def_offset_coef
      write(6,*)"thinkdeb in mg_entry,  ", 6
      call flush(6)
if(present(n_owned_anl).and.present(anl_lonlat1d)) then
call this%def_mg_weights(n_owned_anl=n_owned_anl,lonlat1d_anl=anl_lonlat1d)
else
call this%def_mg_weights
end if
      write(6,*)"thinkdeb in mg_entry,  ", 7
      call flush(6)

if(this%mgbf_line) then
   write(6,*)"thinkdeb init_mg_line is called"
   call this%init_mg_line
end if
      write(6,*)"thinkdeb in mg_entry,  ", 8
      call flush(6)

call this%lsqr_mg_coef
      write(6,*)"thinkdeb in mg_entry,  ", 9
      call flush(6)

call this%lwq_vertical_coef(this%lm_a,this%lm,this%cvf1,this%cvf2,this%cvf3,this%cvf4,this%lref)

      write(6,*)"thinkdeb in mg_entry,  ", 10
      call flush(6)
!***
!*** Just for testing of standalone version. In GSI WORKA will be given
!*** through a separate subroutine
!***



!-----------------------------------------------------------------------
end subroutine mg_initialize

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module subroutine mg_finalize(this)
!**********************************************************************!
!                                                                      !
!   Finalize multigrid Beta Function                                   !
!                                                     M. Rancic (2020) !
!***********************************************************************
implicit none
class (mg_intstate_type), intent(inout) :: this

real(r_kind), allocatable, dimension(:,:):: PA, VA
integer(i_kind):: n,m,L
integer(i_kind):: ierr
integer:: nm,mm,lm
!-----------------------------------------------------------------------

if(this%ldelta) then
   !
   ! Horizontal cross-section
   !
   nm=this%nm
   mm=this%mm
   lm=this%lm
end if

if(this%nxm*this%nym>1) call this%barrierMPI

call this%deallocate_mg_intstate
if (this%mpi_comm_comp /= MPI_COMM_NULL) then
   call MPI_Comm_free(this%mpi_comm_comp, ierr)
end if

!-----------------------------------------------------------------------
end subroutine mg_finalize
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end submodule mg_entrymod
