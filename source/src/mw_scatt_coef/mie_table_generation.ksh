#!/usr/bin/ksh
#
#  To be run from the Mie table source directory (src/mw_scatt_coef/) either interactively:
#
#    ./mie_table_generation.ksh
#
#  or on a supercomputer, e.g.
#
#    llsubmit ./mie_table_generation.ksh
#
# The following options are for the ECMWF IBM and will be ignored elsewhere
#############################################################################
# @ shell          = /usr/bin/ksh
# @ class          = debug
# @ job_type       = serial
# @ job_name       = mie
# @ output         = mie.out
# @ error          = mie.out
# @ resources      = ConsumableCpus(64) ConsumableMemory(6gb)
# @ node           = 1 
# @ tasks_per_node = 1 
# @ cpu_limit      = 30:05:00 
# @ wall_clock_limit = 30:05:00 
# @ core_limit     = 4096
# @ notification   = never
# @ queue
#############################################################################

# Set number of threads if using RTTOV compiled with OpenMP
export OMP_NUM_THREADS=4

# Executable
EXEC=../../bin/rttov_scatt_make_coef.exe


# Remaining options are ECMWF-only

# DEBUG=1 option will only work properly at ECMWF
DEBUG=0

# OMP=1 option to run multithreaded at ECMWF
OMP=0
ECMWF=0

if [[  $DEBUG == 1 ]] ; then
  set -xv 
fi

if [[ $OMP == 1 ]] ; then

  # Multithreading
  if [[ $ECMWF == 1 ]] ; then

    # ECMWF IBM 
    export OMP_NUM_THREADS=64
    export XLSMPOPTS="parthds=8:stack=500000000 : spins=500000 : yields=50000 "  

  fi
fi

if [[ $DEBUG == 0 ]] ; then

  $EXEC

else

  # Debugging on the ECMWF IBM using totalview
  PATH=/usr/local/apps/toolworks/totalview.8.7.0-3/bin:$PATH
  PATH=/usr/local/lpp/poe.5104.petva2.pmr52259:$PATH
  export DISPLAY=hombre:0.0
  totalview $EXEC

fi
