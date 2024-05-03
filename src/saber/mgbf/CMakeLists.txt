# (C) Copyright 2022 United States Government as represented by the Administrator of the National
#     Aeronautics and Space Administration
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


add_subdirectory(mgbf_lib)  # Include the mgbf_lib directory

#set (jbfilenames "mgbf_lib/jp_pbfil.f90" )
set (build_saber_mgbf 1)
if( build_saber_mgbf )
    list(APPEND mgbf_src_files_list

    # Covariance block
    covariance/MGBF_Covariance.h
    covariance/MGBF_Covariance.cc
    covariance/MGBF_Covariance.interface.F90
    covariance/MGBF_Covariance.interface.h
    covariance/mgbf_covariance_mod.f90

    # Grid
#clt    grid/MGBF_Grid.h
#clt    grid/MGBF_Grid.interface.F90
#clt    grid/MGBF_Grid.interface.h
#clt    grid/mgbf_grid_mod.f90

    # Interpolation block
    #clt don't need interpolation block  for being now
#    interpolation/MGBF_Interpolation.h

    # Unstructured interpolation code ported from oops (until new interp code can be used)
#    interpolation/unstructured_interp/saber_unstructured_interpolation_mod.F90
#    interpolation/unstructured_interp/UnstructuredInterpolation.cc
#    interpolation/unstructured_interp/UnstructuredInterpolation.h
#    interpolation/unstructured_interp/UnstructuredInterpolation.interface.F90
#    interpolation/unstructured_interp/UnstructuredInterpolation.interface.h

    # Utilities
#    utils/mgbf_utils_mod.f90

    )
endif()
#clt find_package(mgbf_lib  REQUIRED )
message (STATUS "thinkdeb1 " ${mgbf_src_files_list} )

set( mgbf_src_files

${mgbf_src_files_list}

)
 message (STATUS "thinkdeb2.4"   ${mgbf_src_files} )
add_library(mgbf_interface STATIC ${mgbf_src_files} )
target_link_libraries(mgbf_interface PRIVATE mgbf_lib)
target_include_directories(mgbf_interface
  PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/mgbf_lib/modules)
#target_include_directories(mgbf_interface public ../blocks)
