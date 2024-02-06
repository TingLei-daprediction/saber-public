! (C) Copyright 2022 United States Government as represented by the Administrator of the National
!     Aeronautics and Space Administration
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module mgbf_covariance_mod

! atlas
use atlas_module,                   only: atlas_fieldset, atlas_field
use atlas_module,    only: atlas_functionspace

! fckit
use fckit_mpi_module,               only: fckit_mpi_comm
use fckit_configuration_module,     only: fckit_configuration

! oops
use kinds,                          only: r_kind
use random_mod

! saber
!clt use mgbf_grid_mod,                   only: mgbf_grid
use type_mgbf_mod,                   only: mgbf_type

implicit none
private
public mgbf_covariance


! Fortran class header
type :: mgbf_covariance
  type(mgbf_type) :: mgbf_driver 
  logical :: noMGBF
  logical :: bypassMGBFbe
  logical :: cv   ! cv=.true.; sv=.false.
  integer :: mp_comm_world
  integer :: rank
!clt  integer :: lat2,lon2 ! these belog to mgbf_grid
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: randomize
    procedure, public :: multiply
    procedure, public :: multiply_ad
end type mgbf_covariance

character(len=*), parameter :: myname='mgbf_covariance_mod'

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, comm, config, background, firstguess)

! Arguments
class(mgbf_covariance),     intent(inout) :: self
type(fckit_mpi_comm),      intent(in)    :: comm
type(fckit_configuration), intent(in)    :: config
type(atlas_fieldset),      intent(in)    :: background
type(atlas_fieldset),      intent(in)    :: firstguess

! Locals
character(len=*), parameter :: myname_=myname//'*create'
character(len=:), allocatable :: nml,bef
logical :: central
integer :: layout(2)

type(atlas_field) :: afield
real(kind=r_kind), pointer :: t(:,:)

! Hold communicator
! -----------------
!self%mp_comm_world=comm%communicator()

! Create the grid
! ---------------
!clt call self%grid%create(config, comm)
self%rank = comm%rank()

call config%get_or_die("debuggingxx bypass mgbf", self%noMGBF)
#if 0
if (.not. self%noMGBF) then
  call config%get_or_die("saber central block", central)
  if (.not. central) then
     call abor1_ftn(myname_//": not ready to handle sqrt(B) case")
  endif
  call config%get_or_die("debugging deep bypass mgbf B error", self%bypassMGBFbe)

! Get required name of resources for MGBF B error
! ----------------------------------------------
  call config%get_or_die("mgbf berror namelist file",  nml)
  call config%get_or_die("mgbf error covariance file", bef)

! Initialize MGBF-Berror components
! --------------------------------
  layout=self%grid%layout
! layout=-1
  call mgbfbclim_init(self%cv)
endif
#endif 
call  self%mgbf_driver%mgbf_init()
! Get background (temporary test of the functionality)
!cltafield = background%field('air_temperature')
!clt call afield%data(t)

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

! Arguments
class(mgbf_covariance) :: self

! Locals

#if 0
if (.not. self%noMGBF) then
   call mgbfbclim_final(.false.)
endif
#endif 

! Delete the grid
! ---------------
!clt call self%grid%delete()

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine randomize(self, fields)

! Arguments
class(mgbf_covariance), intent(inout) :: self
type(atlas_fieldset),  intent(inout) :: fields

! Locals
type(atlas_field) :: afield
real(kind=r_kind), pointer :: psi(:,:), chi(:,:), t(:,:), q(:,:), qi(:,:), ql(:,:), o3(:,:)
real(kind=r_kind), pointer :: ps(:)

integer, parameter :: rseed = 3

! Get Atlas field
afield = fields%field('stream_function')
call afield%data(psi)

afield = fields%field('velocity_potential')
call afield%data(chi)

afield = fields%field('air_temperature')
call afield%data(t)

afield = fields%field('surface_pressure')
call afield%data(ps)

afield = fields%field('specific_humidity')
call afield%data(q)

afield = fields%field('cloud_liquid_ice')
call afield%data(qi)

afield = fields%field('cloud_liquid_water')
call afield%data(ql)

afield = fields%field('ozone_mass_mixing_ratio')
call afield%data(o3)


! Set fields to random numbers
call normal_distribution(psi, 0.0_r_kind, 1.0_r_kind, rseed)


end subroutine randomize

! --------------------------------------------------------------------------------------------------

subroutine multiply(self, fields)
! Arguments
class(mgbf_covariance), intent(inout) :: self
type(atlas_fieldset),  intent(inout) :: fields
type(atlas_functionspace) :: afunctionspace

! Locals
type(atlas_field) :: afield
real(kind=r_kind), pointer :: t(:,:)

!clt now noly consider t
!  afield = fields%field('air_temperature')
!  call afield%data(t)
  call self%mgbf_driver%mgbf_apply(fields)
! Halo exchange
!afunctionspace = afield%functionspace()
!call afunctionspace%halo_exchange(afield)
  

end subroutine multiply

! --------------------------------------------------------------------------------------------------

subroutine multiply_ad(self, fields)

! Arguments
class(mgbf_covariance), intent(inout) :: self
type(atlas_fieldset),  intent(inout) :: fields

! This routine only needed when B = G^T G (sqrt-factored)

! To do list for this method
! 1. Convert fields (Atlas fieldsets) to MGBF bundle
! 2. Call MGBF covariance operator adjoint (sqrt version)
!        afield = fields%field('stream_function')
!        call afield%data(var3d)
!        var3d=0.0_r_kind

end subroutine multiply_ad

! --------------------------------------------------------------------------------------------------

end module mgbf_covariance_mod
