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
          call this%mg_entrymod%mg_initialize
 
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
    integer(i_kind)::i,j,k,ij,ii,jj,nx,ny
    nx=nm-n0+1
    ny=mm-m0+1
     afield = fieldset%field('air_temperature')
     call afield%data(t)
     do k=1,km 
     ij=1
     do jj=1,ny
      do ii=1,nx
       i=ii-n0+1
       j=jj-m0+1
       worka(k,i,j)=t(ij,k)
      enddo
     enddo
     enddo

    
    call this%mg_transfer%anal_to_filt_all !cltthink (fieldset)


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
       call mg_filtering_procedure(mgbf_proc)  !cltthink ,fieldset)

!======================================================================

!***
!*** From first generation of filter grid to analysis grid (x-directoin)
!***

          call this%mg_transfer%filt_to_anal_all  !cltthink (fieldset)


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!***
!*** Adjoint test if needed
!***
      
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! Halo exchange

afunctionspace = afield%functionspace()
call afunctionspace%halo_exchange(afield)

   
      
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
