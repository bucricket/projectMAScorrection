#---------------------------------------------------------
# Makefile for C and C++ examples of calling RTTOV
#---------------------------------------------------------

# First compile RTTOV. Then:
# - set FC to the Fortran compiler used to build RTTOV and
#     specify appropriate compiler flags.
# - modify the paths below for the RTTOV libraries and HDF5
#     or NetCDF libraries as necessary.


#--- Compilers -------------------------------------------

# Select Fortran compiler used to compile RTTOV:
FC=gfortran -g -fPIC
# FC=ifort -g -fPIC -nofor_main
# FC=pgf90 -g -fPIC -Mnomain
# FC=nagfor -g -PIC


# C and C++ compilers
CC=cc -w -g -fPIC -std=c99
CPP=c++ -w -g -pg -fPIC -std=c++0x

#---------------------------------------------------------


# --- RTTOV compiled with HDF5 ---------------------------

# RTTOV libraries including HDF5
RTTOV_LIBS=-L../lib -lrttov11_wrapper -lrttov11_other \
                    -lrttov11_emis_atlas -lrttov11_brdf_atlas \
                    -lrttov11_parallel -lrttov11_coef_io \
                    -lrttov11_hdf -lrttov11_main

# HDF5 library required if RTTOV was compiled with HDF5
HDF5_PREFIX  = /path/to/HDF5/install

FFLAGS_HDF5  = -I$(HDF5_PREFIX)/include
LDFLAGS_HDF5 = -L$(HDF5_PREFIX)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz
#---------------------------------------------------------


# --- RTTOV compiled without HDF5 ------------------------

# RTTOV libraries without HDF5
# RTTOV_LIBS=-L../lib -lrttov11_wrapper -lrttov11_other \
#                     -lrttov11_emis_atlas -lrttov11_brdf_atlas \
#                     -lrttov11_parallel -lrttov11_coef_io \
#                     -lrttov11_main

# NetCDF library only required if RTTOV was not compiled with HDF5
# NETCDF_PREFIX  = /path/to/NetCDF/install

# FFLAGS_NETCDF  = -I$(NETCDF_PREFIX)/include
# LDFLAGS_NETCDF = -L$(NETCDF_PREFIX)/lib -lnetcdff -lnetcdf
#---------------------------------------------------------


all : interface_example_c interface_example_cpp Rttov_example RttovSafe_example

clean :
	\rm -f interface_example_c interface_example_cpp Rttov_example RttovSafe_example *.o

interface_example_c:  interface_example_c.c
	$(CC) -I. -c interface_example_c.c
	$(FC) -I. -o interface_example_c interface_example_c.o   $(RTTOV_LIBS) $(LDFLAGS_HDF5) $(LDFLAGS_NETCDF) -lm

interface_example_cpp:  interface_example_cpp.cpp
	$(CPP) -I. -c interface_example_cpp.cpp
	$(FC) -I. -o interface_example_cpp interface_example_cpp.o   $(RTTOV_LIBS) $(LDFLAGS_HDF5) $(LDFLAGS_NETCDF) -lm -lstdc++

Options.o : Options.h Options.cpp
	$(CPP) -I. -c Options.cpp

Profile.o : Profile.h Profile.cpp
	$(CPP) -I. -c Profile.cpp

Profiles.o : Profiles.h Profiles.cpp
	$(CPP) -I. -c Profiles.cpp

Rttov.o : Rttov.cpp Rttov.h Options.h Profiles.h
	$(CPP) -I. -c Rttov.cpp

RttovSafe.o : RttovSafe.cpp RttovSafe.h Options.h Profiles.h Rttov.h
	$(CPP) -I. -c RttovSafe.cpp

Rttov_example:  Rttov_example.cpp Rttov.o Profiles.o Options.o
	$(CPP) -I. -c Rttov_example.cpp
	$(FC) -I. -o Rttov_example Rttov_example.o Rttov.o Options.o Profiles.o ${RTTOV_LIBS} ${LDFLAGS_HDF5} $(LDFLAGS_NETCDF) -lm -lstdc++

RttovSafe_example:  RttovSafe_example.cpp RttovSafe.o Rttov.o Profile.o Options.o
	$(CPP) -I. -c RttovSafe_example.cpp
	$(FC) -I. -o RttovSafe_example RttovSafe_example.o RttovSafe.o Rttov.o Options.o Profiles.o Profile.o ${RTTOV_LIBS} ${LDFLAGS_HDF5} $(LDFLAGS_NETCDF) -lm -lstdc++
