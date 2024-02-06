!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        submodule(mg_intstate) mg_entrymod
!***********************************************************************
!                                                                      !
!   Initialize and finialize multigrid Beta filter for modeling of     !
!   background error covariance                                        !
!                                                     M. Rancic (2020) !
!***********************************************************************
use mpi
use kinds, only: r_kind,i_kind

!cltmoved to mg_parameter integer(i_kind):: km,km2,km3
                        contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module subroutine mg_initialize(this,inputfilename,obj_parameter)
implicit none
!**********************************************************************!
!                                                                      !
!   Initialization subroutine                                          !
!                                                     M. Rancic (2020) !
!***********************************************************************
!type (type_mg_entrymod):: this
class (mg_intstate_type):: this
character*(*),optional,intent(in) :: inputfilename
class(mg_parameter_type),optional,intent(in)::obj_parameter

!cltreal(r_kind), allocatable, dimension(:,:):: PA
!#include  "type_parameter_locpointer.inc"
!#include  "type_parameter_point2this.inc"


!---------------------------------------------------------------------------
!
!               Firs set of subroutines is called only once and serves to 
!               initialte the MGBF run                                 
! 
!---------------------------------------------------------------------------

!****
!**** Initialize run multigrid Beta filter parameters
!****
      if (present(inputfilename)) then  
      call this%init_mg_parameter(inputfilename)
      elseif (present(obj_parameter)) then
      this%mg_parameter_type=obj_parameter
      endif

!****
!**** Initialize MPI
!****

      call this%init_mg_MPI

!***
!*** Initialize integration domain
!***

      call this%init_mg_domain


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

    if(this%l_filt) then
      this%km2 = this%km2_f
      this%km3 = this%km3_f
    else 
      this%km2 = this%km2_e
      this%km3 = this%km3_e
    endif
       write(6,*)'thinkdeb33 ',this%km2,this%km3,this%lm
!cltdebug      this%km2=0;this%km3=0  !cltthinktodo this is not defined in the test case
                              !using
                              !/scratch1/NCEPDEV/da/Miodrag.Rancic/Mars_Jul05_2022/RUN/mgbf.nml_offset  
      this%km = this%km2+this%lm*this%km3

!***
!*** Allocate variables, define weights, prepare mapping 
!*** between analysis and filter grid
!***

      call this%allocate_mg_intstate  !(this%km)  !cltthink

      call this%def_offset_coef

      call this%def_mg_weights

      if( this%mgbf_line) then
         call this%init_mg_line
      endif

      call this%lsqr_mg_coef 

!for now        call lwq_vertical_coef(lm ,lmf,cvf1,cvf2,cvf3,cvf4,lref)
!for now        call lwq_vertical_coef(lmf,lmh,cvh1,cvh2,cvh3,cvh4,lref_h)

!***
!*** Just for testing of standalone version. In GSI WORKA will be given
!*** through a separate subroutine 
!***

!    call input_3d(WORKA(     1:  lm,:,:),1,1,     1,mm,nm,  lm,mm0,4,3)
!    call input_3d(WORKA(  lm+1:2*lm,:,:),1,1,  lm+1,mm,nm,2*lm,mm0,6,5)
!    call input_3d(WORKA(2*lm+1:3*lm,:,:),1,1,2*lm+1,mm,nm,3*lm,mm0,2,1)
!    call input_3d(WORKA(3*lm+1:4*lm,:,:),1,1,3*lm+1,mm,nm,4*lm,mm0,3,2)
!    call input_3d(WORKA(4*lm+1:5*lm,:,:),1,1,4*lm+1,mm,nm,5*lm,mm0,7,3)
!    call input_3d(WORKA(5*lm+1:6*lm,:,:),1,1,5*lm+1,mm,nm,6*lm,mm0,4,5)

!    call input_3d(WORKA(6*lm+1:6*lm+1,:,:),1,1,6*lm+1,mm,nm,6*lm+1,mm0,2,1)
!    call input_3d(WORKA(6*lm+2:6*lm+2,:,:),1,1,6*lm+2,mm,nm,6*lm+2,mm0,4,1)
!    call input_3d(WORKA(6*lm+3:6*lm+3,:,:),1,1,6*lm+3,mm,nm,6*lm+3,mm0,5,1)
!    call input_3d(WORKA(6*lm+4:6*lm+4,:,:),1,1,6*lm+4,mm,nm,6*lm+4,mm0,7,1)

!clt     WORKA(:,:,:)=0.
!TEST
!             call finishMPI
!TEST

!-----------------------------------------------------------------------
                        endsubroutine mg_initialize

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                 module       subroutine mg_finalize(this)
!**********************************************************************!
!                                                                      !
!   Finalize multigrid Beta Function                                   !
!                                                     M. Rancic (2020) !
!***********************************************************************
!clt #use mg_parameter, only: nm,mm
implicit none
class (mg_intstate_type)::this

real(r_kind), allocatable, dimension(:,:):: PA, VA
integer(i_kind):: n,m,L
integer:: nm,mm,lm
!-----------------------------------------------------------------------

if(this%ldelta) then

!
! Horizontal cross-section
!
nm=this%nm
mm=this%mm
lm=this%lm
endif

     call this%barrierMPI


          call this%deallocate_mg_intstate          

!-----------------------------------------------------------------------
                        endsubroutine mg_finalize
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        end submodule mg_entrymod
