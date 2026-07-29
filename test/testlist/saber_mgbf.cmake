message( STATUS "- GSI GFS" )

# dirac_gsi_gfs_global
saber_add_test( TARGET saber_dirac_mgbf-1
                MPI 4 
                OMP 1
                COMMAND ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x
                ARGS testinput/dirac_mgbf-1.yaml
                DEPENDS saber_quench_error_covariance_toolbox.x )
