message( STATUS "- GSI GFS" )

# dirac_gsi_gfs_global
saber_add_test( TARGET saber_dirac_mgbf-1
                MPI 4 
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_mgbf-1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )

# Generate the common ten-member ensemble used by the serial and parallel
# ensemble-localization tests.
saber_add_test( TARGET saber_randomization_mgbf_ensemble
                MPI 4
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/randomization_mgbf_ensemble.yaml
                     testdata/randomization_mgbf_ensemble.log
                DEPENDS saber_quench_error_covariance_toolbox.x )

# Eight-rank reference: ten members localized by one four-rank MGBF component.
saber_add_test( TARGET saber_dirac_mgbf-2
                MPI 8
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_mgbf-2.yaml testdata/dirac_mgbf-2.log
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_mgbf_ensemble )

# The same ten members split across two four-rank MGBF communicators.
saber_add_test( TARGET saber_dirac_mgbf-3
                MPI 8
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_mgbf-3.yaml testdata/dirac_mgbf-3.log
                DEPENDS saber_quench_error_covariance_toolbox.x
                TEST_DEPENDS saber_randomization_mgbf_ensemble )

# Runtime equivalence check, avoiding two independently maintained references.
saber_add_test( TARGET saber_compare_diagnostics_mgbf_ensemble_parallel
                TYPE SCRIPT
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_compare_dirac_diagnostics.py
                ARGS testinput/compare_diagnostics_mgbf_ensemble_parallel.yaml
                TEST_DEPENDS saber_dirac_mgbf-2 saber_dirac_mgbf-3 )
