!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        program RBETA_TEST 
!***********************************************************************
!                                                                      !
!   Multigrid Beta filter for modeling background error covariance     !
!                                                                      !
!                                                     M. Rancic (2020) !
!***********************************************************************
use mpi
use kinds, only: r_kind,i_kind
!clt use mg_entrymod, only: mg_initialize,mg_finalize
!clt use mg_mppstuff, only: finishMPI,mype
!clt use mg_filtering, only: mg_filtering_procedure
!clt use mg_transfer, only: anal_to_filt_all,filt_to_anal_all 
!clt use mg_parameter, only: mgbf_proc
use  mg_intstate
use mg_timers
use mg_input

implicit none
type (mg_intstate_type):: obj_mgbf
type (mg_intstate_type):: obj2_mgbf
real(r_kind), allocatable, dimension(:,:):: PA
real(r_kind), allocatable,dimension(:,:,:):: WORKA
        integer :: mype,unitnum
        character*4 :: file_str
  integer(i_kind):: ierr

!-----------------------------------------------------------------------

                                                   call btim(   total_tim)
                                                   call btim(    init_tim)
  
      call MPI_INIT(ierr)

!***
!*** Initialzie multigrid Beta filter                                   
if(1.gt.0) then
!***
          call obj_mgbf%mg_initialize("mgbeta.nml")

                                                     call etim(    init_tim)
!clt          write(6,*)"worka dim ",obj_mgbf%km,obj_mgbf%n0,obj_mgbf%nm,obj_mgbf%m0,obj_mgbf%mm   
          allocate(WORKA(obj_mgbf%km,obj_mgbf%n0:obj_mgbf%nm,obj_mgbf%m0:obj_mgbf%mm))                               ; WORKA=0.
if(obj_mgbf%ldelta) then

 allocate(PA(1:obj_mgbf%nm,1:obj_mgbf%mm))

    PA = 0.
    call input_spec1_2d(obj_mgbf, PA,obj_mgbf%nxm/2,obj_mgbf%mym/2,'md')

!    WORKA(3*lm+1:4*lm,:,:)=0.
    WORKA(3*obj_mgbf%lm+obj_mgbf%lm/2,:,:)=PA(:,:)


deallocate(PA)

endif
!***
!*** From the analysis to first generation of filter grid
!***
                                                   call btim(    an2filt_tim)

          call obj_mgbf%anal_to_filt_all(WORKA)
                                                   call etim(    an2filt_tim)


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!***
!*** Adjoint test if needed
!***

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!***
!*** Filtering
!***
!======================================================================

       call obj_mgbf%mg_filtering_procedure(obj_mgbf%mgbf_proc)  !clt to be changed

!======================================================================

!***
!*** From first generation of filter grid to analysis grid (x-directoin)
!***

                                                   call btim(   filt2an_tim)
          call obj_mgbf%filt_to_anal_all(WORKA)

                                                   call etim(   filt2an_tim)
        mype=obj_mgbf%mype
        unitnum=25+mype
        write(6,*)WORKA(1,1,1)
        write(file_str,"(I4.4)") mype
   
        open(unit=unitnum,file='mpi'//file_str//'version-worka.bin',access="sequential",form='unformatted',STATUS='replace')
        write(unitnum)WORKA
!clt        if(any(WORKA .gt.0.01)) then
        if(mype==35) then
       write(6,*)'thinkdebworka ',size(WORKA)
       write(6,*)WORKA
       call flush(6)
      endif

        close (unitnum)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!***
!*** Adjoint test if needed
!***
      
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
   
      
!==================== Forward (Smoothing step) ========================
!***
!*** DONE! Deallocate variables
!***
                                                   call btim(   output_tim)
       call obj_mgbf%mg_finalize

                                                   call etim(   output_tim)
                                                   call etim(   total_tim)


        deallocate(WORKA)
!***
!*** Print wall clock and cpu timing
!***
!clt for another obj2_mgbf
endif !1 gt 2



!clt      call obj_mgbf%finishMPI
          write(6,*)'thinkdeb to run for obj_2'
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          call obj2_mgbf%mg_initialize("mgbeta.nml")

          write(6,*)"worka dim2 ",obj2_mgbf%km,obj2_mgbf%n0,obj2_mgbf%nm,obj2_mgbf%m0,obj2_mgbf%mm   
          allocate(WORKA(obj2_mgbf%km,obj2_mgbf%n0:obj2_mgbf%nm,obj2_mgbf%m0:obj2_mgbf%mm))                               ; WORKA=0.
if(obj2_mgbf%ldelta) then

 allocate(PA(1:obj2_mgbf%nm,1:obj2_mgbf%mm))

    PA = 0.
    call input_spec1_2d(obj2_mgbf, PA,obj2_mgbf%nxm/2,obj2_mgbf%mym/2,'md')

!    WORKA(3*lm+1:4*lm,:,:)=0.
    WORKA(3*obj2_mgbf%lm+obj2_mgbf%lm/2,:,:)=PA(:,:)


deallocate(PA)

endif
!***
!*** From the analysis to first generation of filter grid
!***
!                                                   call btim(    an2filt_tim)

          call obj2_mgbf%anal_to_filt_all(WORKA)
                                                   !call etim(    an2filt_tim)


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!***
!*** Adjoint test if needed
!***

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!***
!*** Filtering
!***
!======================================================================

       call obj2_mgbf%mg_filtering_procedure(obj2_mgbf%mgbf_proc)  !clt to be changed

!======================================================================

!***
!*** From first generation of filter grid to analysis grid (x-directoin)
!***

!                                                   call btim(   filt2an_tim)
          call obj2_mgbf%filt_to_anal_all(WORKA)

!                                                   call etim(   filt2an_tim)
        mype=obj2_mgbf%mype
        unitnum=25+mype
        write(6,*)WORKA(1,1,1)
        write(file_str,"(I4.4)") mype
   
        open(unit=unitnum,file='mpi'//file_str//'version-worka.bin',access="sequential",form='unformatted',STATUS='replace')
        write(unitnum)WORKA
!clt        if(any(WORKA .gt.0.01)) then
        if(mype==35) then
       write(6,*)'thinkdebworka2 ',size(WORKA)
       write(6,*)WORKA
       call flush(6)
      endif

        close (unitnum)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!***
!*** Adjoint test if needed
!***
      
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
   
      
!==================== Forward (Smoothing step) ========================
!***
!*** DONE! Deallocate variables
!***
!                                                   call btim(   output_tim)
       call obj2_mgbf%mg_finalize

!                                                   call etim(   output_tim)
!                                                   call etim(   total_tim)


!***
!*** Print wall clock and cpu timing
!***
      call print_mg_timers("version0-timing_cpu.csv", print_cpu, obj2_mgbf%mype)
        
      call MPI_FINALIZE(ierr)


!-----------------------------------------------------------------------
                        endprogram RBETA_TEST 
