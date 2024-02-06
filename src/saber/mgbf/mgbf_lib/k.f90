!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module type_mgbf_mod
!***********************************************************************
!                                                                      !
!   Multigrid Beta filter for modeling background error covariance     !
!                                                                      !
!                                                     M. Rancic (2020) !
!***********************************************************************
use mpi
use kinds, only: r_kind,i_kind
use mg_entrymod, only:mg_entrymod_type, mg_initialize,mg_finalize
use mg_mppstuff, only: finishMPI,mype
use mg_filtering, only: mg_filtering_procedure
use mg_transfer, only: mg_transfer_type,anal_to_filt_all,filt_to_anal_all 
use mg_parameter, only: mgbf_proc
use type_fieldset, only: fieldset_type
implicit none
type mgbf_type 
   type(mg_entrymod_type):: mg_entrymod
   type(mg_transfer_type):: mg_transfer
   contains
    procedure,pass:: mgbf_init
    procedure,pass:: mgbf_apply
    procedure,nopass:: mgbf_finalize
end type mgbf_type
!-----------------------------------------------------------------------
contains

subroutine mgbf_init(this)
    class (mgbf_type),intent(in)::this 
!***
!*** Initialzie multigrid Beta filter                                   
!***
      write(6,*)'thinkdeb in type_mgbf before mg_initialize in mgbf_init'
      call flush(6)
          call this%mg_entrymod%mg_initialize
      write(6,*)'thinkdeb in type_mgbf after mg_initialize in mgbf_init'
      call flush(6)
 
end subroutine mgbf_init

!***
!*** From the analysis to first generation of filter grid
!***
 subroutine mgbf_apply(this,fieldset)
 use mg_intstate,only: worka 
 use atlas_module,    only: atlas_fieldset,atlas_field,atlas_functionspace
 use mg_parameter, only: km,n0,nm,m0,mm
 type(atlas_functionspace) :: afunctionspace
    class (mgbf_type),intent(in):: this
    type(atlas_field) :: afield
    type(atlas_fieldset),intent(inout) :: fieldset !< Fieldset
    real(kind=r_kind), pointer :: t(:,:)
    real(kind=r_kind), allocatable :: trev(:,:,:)
    integer(i_kind)::i,j,k,ij,ii,jj,nx,ny
    integer(i_kind)::jedi_nx,jedi_ny,istart,jstart
  if(mype == 0) then
     jedi_nx=51
     jedi_ny=27  ! they include halo points
  else
     jedi_nx=50
     jedi_ny=26
  endif
    nx=nm-n0  !clt i#+1
    ny=mm-m0  !clt#+1
  istart=(jedi_nx-nx)/2
  jstart=(jedi_ny-ny)/2
  istart=0;jstart=0
   allocate(trev(jedi_nx,jedi_ny,km))
      write(6,*)'thinkdeb in type_mgbf.f90 apply begin '
      write(6,*)'thinkdebtype_mgbf.f90applybe nx..',nx,ny,jedi_nx,jedi_ny,istart,jstart
      call flush(6)
     afield = fieldset%field('air_temperature')
     call afield%data(t)
      write(6,*)'thinkdeb in type_mgbf.f90 apply begin2km ny,nx ',km,ny,nx,size(t,dim=2)
      write(6,*)'thinkdeb in type_mgbf.f90 apply begin2km ny,nx2 ',km,ny,nx,size(t,1)
      call flush(6)   
     do k=1,km 
     ij=1
     do j=1,jedi_ny
      do i=1,jedi_nx
!       i=ii+n0-1
!       j=jj+m0-1
!clt       worka(k,i,j)=t(k,ij)
       trev(i,j,k)=t(k,ij)
       ij=ij+1
      enddo
     enddo
     enddo
     worka=0
      do k=1,km
       do j=1,ny
         do i=1,nx
           worka(k,i,j)=trev(istart+i,jstart+j,k)  ! no flipping of i,j)
         enddo
       enddo
      enddo
!     worka(:,n0,:)=worka(:,n0+1,:)
!     worka(:,:,m0)=worka(:,:,m0)

      write(6,*)'thinkdeb in type_mgbf.f90 apply begin 3 '
      call flush(6)
    
    call this%mg_transfer%anal_to_filt_all !cltthink (fieldset)
      write(6,*)'thinkdeb in type_mgbf.f90 apply begin 4 '
      call flush(6)


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!***
!*** Adjoint test if needed
!***

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!***
!*** Filtering
!***
!======================================================================

!clt       call mgbf_obj_in%mg_transfer%mg_filtering_procedure(mgbf_proc,fieldset)
      write(6,*)'thinkdeb in type_mgbf.f90 apply begin 5 '
      call flush(6)
       call mg_filtering_procedure(mgbf_proc)  !cltthink ,fieldset)

!======================================================================

!***
!*** From first generation of filter grid to analysis grid (x-directoin)
!***
      write(6,*)'thinkdeb in type_mgbf.f90 apply begin 6'
      call flush(6)

          call this%mg_transfer%filt_to_anal_all  !cltthink (fieldset)


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!***
!*** Adjoint test if needed
!***
      
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! Halo exchange
      write(6,*)'thinkdeb in type_mgbf.f90apply begin 7',jedi_ny,ny,worka(10,25,ny)
      write(6,*)'thinkdeb in type_mgbf.f90apply begin 7',jedi_ny,ny,worka(10,25,ny-1)
      call flush(6)
      do k=1,km
      trev(:,:,k)=worka(k,1,1)
      trev(1:nx,jedi_ny,k)=worka(k,1:nx,ny)
      enddo 
      write(6,*)'thinkdeb intype_mgbf.f90apply begin 7.1',jedi_ny,ny,trev(25,jedi_ny,10)
      write(6,*)'thinkdeb intype_mgbf.f90apply begin 7.1',jedi_ny,ny,trev(25,jedi_ny-1,10)
      write(6,*)'thinkdeb intype_mgbf.f90apply begin 7.1',jedi_ny,ny,trev(25,jedi_ny-2,10)
      write(6,*)'thinkdeb intype_mgbf.f90apply begin 7.1',jedi_ny,ny,trev(25,jedi_ny-3,10)
      write(6,*)'thinkdeb intype_mgbf.f90apply begin 7.1',jedi_ny,ny,trev(25,jedi_ny-4,10)
      write(6,*)'thinkdeb intype_mgbf.f90apply begin 7.1',jedi_ny,ny,trev(25,jedi_ny-5,10)
      do k=1,km
       do j=1,ny
         do i=1,nx
!cltbfore           worka(k,i,j)=trev(j,i,k)  ! flipping of i,j)
          trev(istart+i,jstart+j,k)= worka(k,i,j)  !=trev(j,i,k)  ! flipping of i,j)
         enddo
       enddo
      enddo
!clt      trev(1:jedi_nx,jedi_ny,10)=0.093
      write(6,*)'thinkdeb intype_mgbf.f90apply begin 7.2',jedi_ny,ny,trev(25,jedi_ny,10)
      write(6,*)'thinkdeb intype_mgbf.f90apply begin 7.2',jedi_ny,ny,trev(25,jedi_ny-1,10)
      write(6,*)'thinkdeb intype_mgbf.f90apply begin 7.2',jedi_ny,ny,trev(25,jedi_ny-2,10)
      write(6,*)'thinkdeb intype_mgbf.f90apply begin 7.2',jedi_ny,ny,trev(25,jedi_ny-3,10)
      write(6,*)'thinkdeb intype_mgbf.f90apply begin 7.2',jedi_ny,ny,trev(25,jedi_ny-4,10)
      write(6,*)'thinkdeb intype_mgbf.f90apply begin 7.2',jedi_ny,ny,trev(25,jedi_ny-5,10)
     do k=1,km 
     ij=1
     do j=1,jedi_ny
      do i=1,jedi_nx
!       i=ii+n0-1
!       j=jj+m0-1
!cltbefore       trev(i,j,k)=t(k,ij)
       t(k,ij)=trev(i,j,k)
       ij=ij+1
      enddo
     enddo
     enddo

deallocate(trev)
afunctionspace = afield%functionspace()
call afunctionspace%halo_exchange(afield)

      write(6,*)'thinkdeb in type_mgbf.f90 apply begin 8 '
      call flush(6)

   
      
!==================== Forward (Smoothing step) ========================
!***
!*** DONE! Deallocate variables
end subroutine mgbf_apply
!***
subroutine mgbf_finalize(this)
    class (mgbf_type),intent(in)::this 
       call this%mg_entrymod%mg_finalize
end subroutine mgbf_finalize


!-----------------------------------------------------------------------
end module type_mgbf_mod
