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
use mg_entrymod, only: mg_initialize,mg_finalize
use mg_mppstuff, only: finishMPI,mype
use mg_filtering, only: mg_filtering_procedure
use mg_transfer, only: anal_to_filt_all,filt_to_anal_all 
use mg_parameter, only: mgbf_proc
use mg_timers

implicit none

!-----------------------------------------------------------------------

                                                   call btim(   total_tim)
                                                   call btim(    init_tim)
!***
!*** Initialzie multigrid Beta filter                                   
!***
          call mg_initialize

                                                   call etim(    init_tim)
!***
!*** From the analysis to first generation of filter grid
!***
                                                   call btim(    an2filt_tim)

          call anal_to_filt_all
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

       call mg_filtering_procedure(mgbf_proc)

!======================================================================

!***
!*** From first generation of filter grid to analysis grid (x-directoin)
!***

                                                   call btim(   filt2an_tim)
          call filt_to_anal_all

                                                   call etim(   filt2an_tim)

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
       call mg_finalize

                                                   call etim(   output_tim)
                                                   call etim(   total_tim)


!***
!*** Print wall clock and cpu timing
!***
      call print_mg_timers("timing_cpu.csv", print_cpu)



      call finishMPI


!-----------------------------------------------------------------------
                        endprogram RBETA_TEST 
