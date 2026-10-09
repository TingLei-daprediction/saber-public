message( STATUS "- MGBF" )

# Basic MGBF Dirac test.
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
# These parallel-ensemble tests also require the companion OOPS change that
# preserves Atlas field metadata during redistribution to a subcommunicator.
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

# -----------------------------------------------------------------------------
# l_loc_filter_once: filter each variable group once in localization.
#   V2 (equivalence): every switch-on run must reproduce the switch-off run of
#       the same build, value by value over the whole B*Dirac output file
#       (tools/saber_compare_mgbf_filter_once.py). The absolute part of the
#       bound is scaled by the summed single-Dirac switch-off outputs. No
#       reference files are needed.
#   V3 (guards): invalid setups must abort with the expected message.
# Namelists are derived at configure time from the dirac_mgbf.nml test data by
# appending items just before the line that closes &parameters_mgbeta (a later
# namelist item overrides an earlier one). The V2 runs force
# l_for_localization=.true., km2=0 and an explicit km3 so off and on runs see
# the same setup. MGBF-internal scales (nscale > 1) are not exercised in V2.
# -----------------------------------------------------------------------------
set( MGBF_FO_PYTHON "python3" CACHE STRING
     "Python 3 interpreter provided by the user for the MGBF filter-once compare tests (standard library only)" )

set( _fo_src     ${jedi_model_data_tier_1_saber}/dirac_mgbf.nml )
set( _fo_data    ${CMAKE_CURRENT_BINARY_DIR}/testdata )
set( _fo_toolbox ${CMAKE_BINARY_DIR}/bin/saber_quench_error_covariance_toolbox.x )
set( _fo_on      "l_loc_filter_once = .true.," )
set( _fo_loc     "l_for_localization = .true., km2 = 0," )

# Write testdata/<out_name>: dirac_mgbf.nml with <items> appended to &parameters_mgbeta.
# Sets _FO_NML_OK to FALSE in the caller when the closing line cannot be found.
function( mgbf_fo_group_nml out_name items )
  file( READ ${_fo_src} _text )
  string( TOLOWER "${_text}" _lower )
  string( FIND "${_lower}" "&parameters_mgbeta" _head )
  if( _head EQUAL -1 )
    set( _FO_NML_OK FALSE PARENT_SCOPE )
    return()
  endif()
  string( SUBSTRING "${_lower}" ${_head} -1 _tail )
  string( REGEX MATCH "\n[ \t]*/" _close "${_tail}" )
  if( NOT _close )
    set( _FO_NML_OK FALSE PARENT_SCOPE )
    return()
  endif()
  string( FIND "${_tail}" "${_close}" _pos )
  math( EXPR _cut "${_head} + ${_pos}" )
  string( SUBSTRING "${_text}" 0 ${_cut} _before )
  string( SUBSTRING "${_text}" ${_cut} -1 _after )
  file( WRITE ${_fo_data}/${out_name} "${_before}\n  ${items}${_after}" )
endfunction()

# Write testdata/<out_name>: an SDL/VDL init namelist (&parameters_mgbf_init).
function( mgbf_fo_init_nml out_name nscale nvargrp groups cor ivargroup )
  file( WRITE ${_fo_data}/${out_name}
"&parameters_mgbf_init
  nscale = ${nscale}, nvargrp = ${nvargrp},
  readin_mgbf_nml_group = ${groups},
  readin_multigrp_cor = ${cor},
  readin_ivargroup = ${ivargroup},
/
" )
endfunction()

# V2 Dirac run from testinput/mgbf_filter_once_dirac.yaml.in.
#   input: all (Diracs in all three variables) or sf / vp / ps (one Dirac in
#          streamfunction, velocity potential or surface pressure)
function( mgbf_fo_dirac_test name nml_key nml_file input omp adjoint_test )
  set( FO_NAME ${name} )
  set( FO_NML_KEY ${nml_key} )
  set( FO_NML_FILE ${nml_file} )
  set( FO_ADJOINT_TEST ${adjoint_test} )
  if( input STREQUAL "all" )
    set( FO_DIRAC_LON "[135.0, 120.0, 150.0]" )
    set( FO_DIRAC_LAT "[40.0, 30.0, 45.0]" )
    set( FO_DIRAC_LEVEL "[30, 10, 1]" )
    set( FO_DIRAC_VARIABLE "[air_horizontal_streamfunction, air_horizontal_velocity_potential, air_pressure_at_surface]" )
  elseif( input STREQUAL "sf" )
    set( FO_DIRAC_LON "[135.0]" )
    set( FO_DIRAC_LAT "[40.0]" )
    set( FO_DIRAC_LEVEL "[30]" )
    set( FO_DIRAC_VARIABLE "[air_horizontal_streamfunction]" )
  elseif( input STREQUAL "vp" )
    set( FO_DIRAC_LON "[120.0]" )
    set( FO_DIRAC_LAT "[30.0]" )
    set( FO_DIRAC_LEVEL "[10]" )
    set( FO_DIRAC_VARIABLE "[air_horizontal_velocity_potential]" )
  elseif( input STREQUAL "ps" )
    set( FO_DIRAC_LON "[150.0]" )
    set( FO_DIRAC_LAT "[45.0]" )
    set( FO_DIRAC_LEVEL "[1]" )
    set( FO_DIRAC_VARIABLE "[air_pressure_at_surface]" )
  else()
    message( FATAL_ERROR "mgbf_fo_dirac_test: unknown input '${input}'" )
  endif()
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/testinput/mgbf_filter_once_dirac.yaml.in
                  ${CMAKE_CURRENT_BINARY_DIR}/testinput/${name}.yaml @ONLY )
  file( MAKE_DIRECTORY ${_fo_data}/${name} )
  ecbuild_add_test( TARGET saber_${name}
                    MPI 4
                    OMP ${omp}
                    COMMAND ${_fo_toolbox}
                    ARGS testinput/${name}.yaml testdata/${name}.log
                    DEPENDS saber_quench_error_covariance_toolbox.x
                    TEST_DEPENDS saber_randomization_mgbf_ensemble )
endfunction()

# V2 compare: <result> run against <reference> run; the remaining arguments are
# the single-Dirac switch-off runs that scale the absolute part of the bound.
function( mgbf_fo_compare_test name reference result )
  set( _args --reference testdata/${reference} --result testdata/${result} )
  set( _deps saber_${reference} saber_${result} )
  foreach( contribution ${ARGN} )
    list( APPEND _args --contribution testdata/${contribution} )
    list( APPEND _deps saber_${contribution} )
  endforeach()
  ecbuild_add_test( TARGET saber_${name}
                    TYPE SCRIPT
                    COMMAND ${MGBF_FO_PYTHON}
                    ARGS ${CMAKE_BINARY_DIR}/bin/saber_compare_mgbf_filter_once.py
                         ${_args} --rtol 1.0e-12 --atol 1.0e-12
                    TEST_DEPENDS ${_deps} )
endfunction()

# One V2 configuration <cfg>: switch-off and switch-on runs for every input,
# each switch-on run compared with its switch-off run. With omp2 TRUE, an
# OMP 2 switch-on run of the combined input is also compared with the OMP 1
# switch-off run (the MGBF loops have no reductions, so threads do not change
# the arithmetic).
function( mgbf_fo_v2_config cfg nml_key nml_off nml_on adjoint_test omp2 )
  foreach( input all sf vp ps )
    mgbf_fo_dirac_test( mgbf_fo_${cfg}_${input}_off "${nml_key}" ${nml_off} ${input} 1 ${adjoint_test} )
    mgbf_fo_dirac_test( mgbf_fo_${cfg}_${input}_on  "${nml_key}" ${nml_on}  ${input} 1 ${adjoint_test} )
  endforeach()
  set( _contributions mgbf_fo_${cfg}_sf_off mgbf_fo_${cfg}_vp_off mgbf_fo_${cfg}_ps_off )
  mgbf_fo_compare_test( compare_mgbf_fo_${cfg}_all
                        mgbf_fo_${cfg}_all_off mgbf_fo_${cfg}_all_on ${_contributions} )
  foreach( input sf vp ps )
    mgbf_fo_compare_test( compare_mgbf_fo_${cfg}_${input}
                          mgbf_fo_${cfg}_${input}_off mgbf_fo_${cfg}_${input}_on
                          mgbf_fo_${cfg}_${input}_off )
  endforeach()
  if( omp2 )
    mgbf_fo_dirac_test( mgbf_fo_${cfg}_all_on_omp2 "${nml_key}" ${nml_on} all 2 ${adjoint_test} )
    mgbf_fo_compare_test( compare_mgbf_fo_${cfg}_all_omp2
                          mgbf_fo_${cfg}_all_off mgbf_fo_${cfg}_all_on_omp2 ${_contributions} )
  endif()
endfunction()

# V3 guard run from testinput/mgbf_filter_once_guard.yaml.in; passes when the
# expected abort message appears.
function( mgbf_fo_guard_test name nml_key nml_file expected )
  set( FO_NAME ${name} )
  set( FO_NML_KEY ${nml_key} )
  set( FO_NML_FILE ${nml_file} )
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/testinput/mgbf_filter_once_guard.yaml.in
                  ${CMAKE_CURRENT_BINARY_DIR}/testinput/${name}.yaml @ONLY )
  file( MAKE_DIRECTORY ${_fo_data}/${name} )
  ecbuild_add_test( TARGET saber_${name}
                    MPI 4
                    OMP 1
                    COMMAND ${_fo_toolbox}
                    ARGS testinput/${name}.yaml testdata/${name}.log
                    DEPENDS saber_quench_error_covariance_toolbox.x )
  if( TEST saber_${name} )
    set_tests_properties( saber_${name} PROPERTIES PASS_REGULAR_EXPRESSION "${expected}" )
  endif()
endfunction()

if( NOT EXISTS ${_fo_src} )
  message( WARNING "MGBF filter-once tests skipped: ${_fo_src} not found" )
else()
  set( _FO_NML_OK TRUE )
  set( _nml   "mgbf namelist file" )
  set( _init  "mgbf sdl and vdl init namelist file" )

  # One group: all three variables (km3 = 3)
  mgbf_fo_group_nml( mgbf_fo_1grp_off.nml "${_fo_loc} km3 = 3," )
  mgbf_fo_group_nml( mgbf_fo_1grp_on.nml  "${_fo_loc} km3 = 3, ${_fo_on}" )

  # Two groups: (streamfunction, velocity potential) and (surface pressure),
  # non-zero cross-group weights
  mgbf_fo_group_nml( mgbf_fo_grp1_off.nml "${_fo_loc} km3 = 2," )
  mgbf_fo_group_nml( mgbf_fo_grp2_off.nml "${_fo_loc} km3 = 1," )
  mgbf_fo_group_nml( mgbf_fo_grp1_on.nml  "${_fo_loc} km3 = 2, ${_fo_on}" )
  mgbf_fo_group_nml( mgbf_fo_grp2_on.nml  "${_fo_loc} km3 = 1, ${_fo_on}" )
  mgbf_fo_init_nml( mgbf_fo_2grp_off_init.nml 1 2
                    "'testdata/mgbf_fo_grp1_off.nml', 'testdata/mgbf_fo_grp2_off.nml'"
                    "1.0, 0.5, 0.5, 1.0" "2, 3" )
  mgbf_fo_init_nml( mgbf_fo_2grp_on_init.nml 1 2
                    "'testdata/mgbf_fo_grp1_on.nml', 'testdata/mgbf_fo_grp2_on.nml'"
                    "1.0, 0.5, 0.5, 1.0" "2, 3" )

  # V3 namelists
  mgbf_fo_group_nml( mgbf_fo_guard_static_b.nml "l_for_localization = .false., km2 = 0, km3 = 3, ${_fo_on}" )
  mgbf_fo_group_nml( mgbf_fo_guard_km2.nml      "l_for_localization = .true., km2 = 1, km3 = 2, ${_fo_on}" )
  mgbf_fo_group_nml( mgbf_fo_guard_km3.nml      "${_fo_loc} km3 = 0, ${_fo_on}" )
  mgbf_fo_group_nml( mgbf_fo_guard_n_ens.nml    "${_fo_loc} km3 = 3, n_ens = 2, ${_fo_on}" )
  mgbf_fo_group_nml( mgbf_fo_guard_scale2.nml   "${_fo_loc} km3 = 2, ${_fo_on}" )
  # group 1 on, group 2 off: also proves that an omitted switch is not inherited
  mgbf_fo_init_nml( mgbf_fo_guard_mixed_init.nml 1 2
                    "'testdata/mgbf_fo_grp1_on.nml', 'testdata/mgbf_fo_grp2_off.nml'"
                    "1.0, 0.0, 0.0, 1.0" "2, 3" )
  mgbf_fo_init_nml( mgbf_fo_guard_ivargroup_init.nml 1 2
                    "'testdata/mgbf_fo_grp1_on.nml', 'testdata/mgbf_fo_grp2_on.nml'"
                    "1.0, 0.0, 0.0, 1.0" "1, 3" )
  mgbf_fo_init_nml( mgbf_fo_guard_scales_init.nml 2 1
                    "'testdata/mgbf_fo_1grp_on.nml', 'testdata/mgbf_fo_guard_scale2.nml'"
                    "1.0" "3" )

  if( NOT _FO_NML_OK )
    message( WARNING "MGBF filter-once tests skipped: could not locate the closing line of "
                     "&parameters_mgbeta in ${_fo_src}" )
  else()
    # V2
    # One group (single-group branch)
    mgbf_fo_v2_config( 1grp "${_nml}" testdata/mgbf_fo_1grp_off.nml testdata/mgbf_fo_1grp_on.nml
                       true TRUE )
    # Two groups with cross-group weights 0.5 (multi-group branch). These weights
    # make L non-symmetric, so no adjoint test here.
    mgbf_fo_v2_config( 2grp "${_init}" testdata/mgbf_fo_2grp_off_init.nml testdata/mgbf_fo_2grp_on_init.nml
                       false TRUE )
    # Not testable: one group with a non-unit multigrp_cor(1,1). That weight is
    # only set through the init namelist, and with nscale = nvargrp = 1 create()
    # replaces the group namelist by the init file itself (pre-existing), so the
    # case only arises with MGBF-internal scales, which V2 does not exercise.

    # V3
    mgbf_fo_guard_test( mgbf_fo_guard_static_b "${_nml}" testdata/mgbf_fo_guard_static_b.nml
                        "l_loc_filter_once: requires l_for_localization=.true.; a static B" )
    mgbf_fo_guard_test( mgbf_fo_guard_km2 "${_nml}" testdata/mgbf_fo_guard_km2.nml
                        "l_loc_filter_once: km2 must be 0" )
    mgbf_fo_guard_test( mgbf_fo_guard_km3 "${_nml}" testdata/mgbf_fo_guard_km3.nml
                        "l_loc_filter_once: km3 .variables in this group. must be >= 1" )
    mgbf_fo_guard_test( mgbf_fo_guard_n_ens "${_nml}" testdata/mgbf_fo_guard_n_ens.nml
                        "l_loc_filter_once: n_ens must be 1" )
    mgbf_fo_guard_test( mgbf_fo_guard_mixed "${_init}" testdata/mgbf_fo_guard_mixed_init.nml
                        "l_loc_filter_once must be the same in every group namelist" )
    mgbf_fo_guard_test( mgbf_fo_guard_ivargroup "${_init}" testdata/mgbf_fo_guard_ivargroup_init.nml
                        "readin_ivargroup does not match the km3 of the group namelists" )
    mgbf_fo_guard_test( mgbf_fo_guard_scales "${_init}" testdata/mgbf_fo_guard_scales_init.nml
                        "variables in group differ from scale 1" )
  endif()
endif()
