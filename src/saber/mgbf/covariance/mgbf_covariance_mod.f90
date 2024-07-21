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
use mgbf_kinds,                          only: r_kind,i_kind
use random_mod

! saber
!clt use mgbf_grid_mod,                   only: mgbf_grid
use mg_intstate , only:            mg_intstate_type

implicit none
private
public mgbf_covariance


! Fortran class header
type :: mgbf_covariance
  type(mg_intstate_type) :: intstate 
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
character(len=:), allocatable :: mgbf_nml,centralblockname
logical :: central
integer :: layout(2)

type(atlas_field) :: afield

! Hold communicator
! -----------------
!self%mp_comm_world=comm%communicator()

! Create the grid
! ---------------
!clt call self%grid%create(config, comm)
self%rank = comm%rank()

!clt call config%get_or_die("debuggingxx bypass mgbf", self%noMGBF)
call config%get_or_die("mgbf namelist file ",  mgbf_nml)
!if (.not. self%noMGBF) then
  call config%get_or_die("saber block name", centralblockname)
!  if (.not. central) then
!     call abor1_ftn(myname_//": not ready to handle sqrt(B) case")
!  endif
!  call config%get_or_die("debugging deep bypass mgbf B error", self%bypassMGBFbe)

! Get required name of resources for MGBF B error
! ----------------------------------------------
!  call config%get_or_die("mgbf berror namelist file",  mgbf_nml)
!//  call config%get_or_die("mgbf error covariance file", bef)

! Initialize MGBF-Berror components
! --------------------------------
! layout=-1
!clt endif
call  self%intstate%mg_initialize(mgbf_nml)  !mgbf_nml like mgbeta.nml
! Get background (temporary test of the functionality)
!cltafield = background%field('air_temperature')
!clt call afield%data(t)

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

! Arguments
class(mgbf_covariance) :: self

! Locals

!clt //if (.not. self%noMGBF) then
   call self%intstate%mg_finalize()
!clt endif

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
real(kind=r_kind), pointer :: ptr_2d(:,:)
real(kind=r_kind), pointer :: ptr_3d(:,:,:)
integer(kind=i_kind):: nz,ilev,isize
real(kind=r_kind), allocatable :: work_mgbf(:,:,:)
real(kind=r_kind), allocatable :: work_mgbf2(:,:,:)
real(kind=r_kind), allocatable :: work2d_mgbf(:,:)
integer(kind=i_kind) :: dim2d(2),dim3d(3)
integer(kind=i_kind):: myrank,nxloc,nyloc,nzloc
integer(kind=i_kind):: i,j,k,ij
integer(kind=i_kind):: n2d
logical  ::  l3d_encountered  


!clt now noly consider t
!  afield = fields%field('air_temperature')
!  call afield%data(t)
!*** From the analysis to first generation of filter grid
          n2d=0
          l3d_encountered=.false.
          allocate(work_mgbf(self%intstate%km_a_all,self%intstate%nm,self%intstate%mm))
          allocate(work_mgbf2(self%intstate%km_a_all,self%intstate%nm,self%intstate%mm))
          allocate(work2d_mgbf(self%intstate%km_a_all,self%intstate%nm*self%intstate%mm))
          dim2d=shape(work2d_mgbf)
          dim3d=shape(work_mgbf)
          nxloc=dim3d(2)
          nyloc=dim3d(3)
          nzloc=dim3d(1)
          work_mgbf2=0.0

             ilev=1
          do isize=1,fields%size()
             
             afield= fields%field(isize)  !clttodo
             if(afield%rank() == 2)  then
               nz=afield%levels()
             write(6,*)'thinkdeb55 isize/name/nz is ',isize,' ',afield%name(),' ',nz
               call afield%data(ptr_2d)
               work2d_mgbf(ilev:ilev+nz-1,:)=ptr_2d 
               ilev=ilev+nz
               if(nz >  1) l3d_encountered=.true.
               if(nz == 1) then 
                  if(l3d_encountered )  stop  !  is required 2d fields are saved consecutively 
                 n2d=n2d+1
               endif
       do k=1,64 ! #nzloc
         do i=1,nxloc*nyloc
             if(ptr_2d(k,i) .gt.0.001) then 
                write(6,*)'thinkdeb666ptr, non zeror k,ij work2d_mgbf ',i,k+64*(isize-1),' ',ptr_2d(k,i)
             endif
         enddo 
       enddo 
             elseif (afield%rank() == 3) then  
               write(6,*)'this case needs more work, stop' ! a better exption handling to be added
               call flush(6)
               stop 
               call afield%data(ptr_3d)
               nz=afield%levels()
               work_mgbf(ilev:ilev+nz-1,:,:)=ptr_3d 
               ilev=ilev+nz
             else
               write(6,*)'wrong in mgbf_covariance_mod.f90 ' !todo  
               stop
             endif 
          enddo
          if(self%intstate%km2.ne.n2d) then 
             write(6,*)'The numbers of 2d variables is different from  mgbf-expected ,stop'
             stop   ! a better exception handling is to be added
          endif
       do k=1,nzloc
          work_mgbf(k,:,:) =reshape(work2d_mgbf(k,:),[dim3d(2),dim3d(3)])
       enddo

       write(6,*)"thinkdeb666-1"
       do k=1,nzloc
         do i=1,nxloc*nyloc
             if(work2d_mgbf(k,i) .gt.0.001) then 
                write(6,*)'thinkdeb666, non zeror k,ij work2d_mgbf ',i,k,' ',work2d_mgbf(k,i)
             endif
         enddo 
       enddo 
       do k=1,nzloc
         do j=1,nxloc
           do i=1,nxloc
             if(work_mgbf(k,i,j) .gt.0.001) then 
                write(6,*)'thinkdeb666, non zeror k,i,j work_mgbf ',i,j,k,' ',work_mgbf(k,i,j)
             endif
           enddo
         enddo 
       enddo 
          call self%intstate%anal_to_filt_allmap(work_mgbf)
         if(1.gt.0) then
          call self%intstate%filtering_procedure(self%intstate%mgbf_proc,1)
         endif
         
!cltorg          call self%intstate%filt_to_anal_allmap(work_mgbf)
          call self%intstate%filt_to_anal_allmap(work_mgbf2)
          work_mgbf=work_mgbf2
        do k=1,nzloc
          work2d_mgbf(k,:)=reshape(work_mgbf(k,:,:),[dim2d(2)])
        enddo
             ilev=1
          do isize=1,fields%size()
  
             afield=fields%field(isize)  !clttodo
             if(afield%rank() == 2) then 
               call afield%data(ptr_2d)
               nz=afield%levels()
             write(6,*)'thinkdeb552 isize/name/nz is ',isize,' ',afield%name(),' ',nz
               ptr_2d(1:nz,:)=work2d_mgbf(ilev:ilev+nz-1,:) 
               ilev=ilev+nz
             elseif (afield%rank() == 3) then  
               call afield%data(ptr_3d)
               nz=afield%levels()
               write(6,*)'wrong in mgbf_covariance_mod.f90 todo ' !todo  
               stop
                 

!clt               ptr_3d=work2d_mgbf(ilev:ilev+nz-1,:) 
               ilev=ilev+nz
             else
               write(6,*)'wrong in mgbf_covariance_mod.f90 ' !todo  
               stop
             endif 
           enddo





          deallocate(work_mgbf)
          deallocate(work2d_mgbf)

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
