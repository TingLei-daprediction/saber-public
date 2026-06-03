# EMC/NCEP/NOAA 2024
# Explicit ordering required: parent modules must compile before their submodules
# so Intel Fortran can find the correct .smod files.
#
include(GNUInstallDirs)

# Dependency tree:
#   mg_parameter.f90  ->  jp_pbfil.f90, mg_domain.f90, mg_domain_loc.f90,
#                         mg_mppstuff.f90
#   mg_intstate.f90   ->  mg_bocos.f90, mg_entrymod.f90, mg_filtering.f90,
#                         mg_generations.f90, mg_interpolate.f90, mg_transfer.f90
set( mgbf_lib_src_files
  # standalone utility modules (no submodule dependencies)
  mgbf_kind.f90
  mgbf_util.f90
  jp_pkind.f90
  jp_pkind2.f90
  jp_pietc.f90
  jp_pietc_s.f90
  jp_pmat.f90
  jp_pmat4.f90
  jp_pbfil2.f90
  jp_pbfil3.f90
  phint.f90
  phint1.f90
  mg_timers.f90
  mg_tools.f90
  # parent modules (must come before their submodules)
  mg_parameter.f90
  mg_intstate.f90
  # submodules of mg_parameter
  jp_pbfil.f90
  mg_domain.f90
  mg_domain_loc.f90
  mg_mppstuff.f90
  # submodules of mg_intstate
  mg_entrymod.f90
  mg_bocos.f90
  mg_filtering.f90
  mg_generations.f90
  mg_interpolate.f90
  mg_transfer.f90
)
add_library(mgbflib STATIC ${mgbf_lib_src_files})
set_target_properties(mgbflib PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
target_include_directories(mgbflib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>)
target_link_libraries(mgbflib PUBLIC MPI::MPI_Fortran)
install(TARGETS mgbflib
        EXPORT saber-targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
