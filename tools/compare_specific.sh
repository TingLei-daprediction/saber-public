#!/bin/sh
#----------------------------------------------------------------------
# Shell script: compare_dirac
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
test1=$1
test2=$2
suffix=$3

# Initialize exit status
status=0

# BUMP tests
if test "${test1%%_*}" = "bump" ; then
   # NCCMP Parameters
   nthreads=4
   tolerance=1.e-5

   # Build file names
   file1=testref/${test1}/test_1-1_${suffix}.nc
   file2=testref/${test2}/test_1-1_${suffix}.nc

   # Compare files with NCCMP
   if type "nccmp" > /dev/null ; then
      echo "Command: nccmp -dfFmqS --threads=${nthreads} -T ${tolerance} ${file1} ${file2}"
      nccmp -dfFmqS --threads=${nthreads} -T ${tolerance} ${file1} ${file2}
      exit_code=$?
      if test "${exit_code}" != "0" ; then
         echo "\e[31mTest failed checking: "${file#testdata/}"\e[0m"
         status=1
      fi
   else
      echo "Cannot find command: nccmp"
   fi
fi

# Exit
exit ${status}
