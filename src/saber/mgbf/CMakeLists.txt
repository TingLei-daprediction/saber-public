# (C) Copyright 2024 DOC/NOAA
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
set (build_saber_mgbf 1)
if( build_saber_mgbf )
    list(APPEND mgbf_src_files_list

    # MGBF library modules. Keep parent modules before their submodules so
    # Intel Fortran can generate the .mod/.smod files before they are used.
    mgbf_lib/mgbf_kind.f90
    mgbf_lib/mgbf_util.f90
    mgbf_lib/jp_pkind.f90
    mgbf_lib/jp_pkind2.f90
    mgbf_lib/jp_pietc.f90
    mgbf_lib/jp_pietc_s.f90
    mgbf_lib/jp_pmat.f90
    mgbf_lib/jp_pmat4.f90
    mgbf_lib/jp_pbfil2.f90
    mgbf_lib/jp_pbfil3.f90
    mgbf_lib/phint.f90
    mgbf_lib/phint1.f90
    mgbf_lib/mg_timers.f90
    mgbf_lib/mg_tools.f90

    # Parent modules.
    mgbf_lib/mg_parameter.f90
    mgbf_lib/mg_intstate.f90

    # Submodules of mg_parameter.
    mgbf_lib/jp_pbfil.f90
    mgbf_lib/mg_domain.f90
    mgbf_lib/mg_domain_loc.f90
    mgbf_lib/mg_mppstuff.f90

    # Submodules of mg_intstate.
    mgbf_lib/mg_entrymod.f90
    mgbf_lib/mg_bocos.f90
    mgbf_lib/mg_filtering.f90
    mgbf_lib/mg_generations.f90
    mgbf_lib/mg_interpolate.f90
    mgbf_lib/mg_transfer.f90

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
