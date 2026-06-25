# (C) Copyright 2026 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


# Set compiler flags for basic build types,
# for compilers where this is not provided by ecbuild.
include(build_type_compiler_flags)

# Set JEDI's common compiler flags
include(jedi_common_compiler_flags)

# Set SABER-specific compiler flags
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  list(APPEND _gnu_fortran_debug_flags
    "-fcheck=all"
    "-fimplicit-none"
    "-ffpe-trap=invalid,zero,overflow,denormal"
    "-finit-derived"
    "-finit-integer=-999"
    "-finit-real=snan")
  string(JOIN " " _gnu_fortran_debug_flags ${_gnu_fortran_debug_flags})
  ecbuild_add_fortran_flags("${_gnu_fortran_debug_flags}" BUILD DEBUG)
endif()
if(HAVE_WARNING)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    ecbuild_add_fortran_flags("-Wall -Wextra")
  endif()
endif()

# --------------------------------------------------------------------------
# SABER debug instrumentation knob (Intel). Edit the flag string below.
# Appended last, so it overrides earlier flags (e.g. -heap-arrays 0 beats
# the inherited -heap-arrays 32). Verify with:
#   grep Fortran_FLAGS <build>/saber/src/saber/CMakeFiles/saber.dir/flags.make
# --------------------------------------------------------------------------
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  ecbuild_add_fortran_flags("-heap-arrays 0" BUILD DEBUG)
  # Optional extra diagnostics (uncomment one at a time):
  #   ecbuild_add_fortran_flags("-gen-interfaces -warn interfaces") # catch dummy/actual mismatches at compile
  #   ecbuild_add_fortran_flags("-fsanitize=address -fno-omit-frame-pointer") # also needs -DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address
endif()
