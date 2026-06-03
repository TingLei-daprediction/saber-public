# (C) Copyright 2024 DOC/NOAA
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
set (build_saber_mgbf 1)
if( build_saber_mgbf )
    add_subdirectory(mgbf_lib)

    list(APPEND mgbf_src_files_list
    # Covariance block
    covariance/MGBF_Covariance.h
    covariance/MGBF_Covariance.cc
    covariance/MGBF_Covariance.interface.F90
    covariance/MGBF_Covariance.interface.h
    covariance/mgbf_covariance_mod.f90

    # Interpolation block

    # Unstructured interpolation code ported from oops (until new interp code can be used)

    # Utilities

    )
endif()

set( mgbf_src_files

${mgbf_src_files_list}

PARENT_SCOPE
)
