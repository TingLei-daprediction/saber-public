# (C) Copyright 2022 United States Government as represented by the Administrator of the National
#     Aeronautics and Space Administration
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
file(GLOB jbfiles mgbf_lib/*.f90)
set (jbfilenames "")
foreach ( _fname ${jbfiles} ) 
   get_filename_component( basefilename ${_fname} NAME )
   list ( APPEND jbfilenames  mgbf_lib/${basefilename} )
endforeach ()
set (build_saber_mgbf 1)
if( build_saber_mgbf )
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
${jbfilenames}

PARENT_SCOPE
)
