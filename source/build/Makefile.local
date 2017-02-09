#!/bin/bash

echo
echo '========================================================================='
echo ' RTTOV compilation script'
echo '========================================================================='
echo
echo "This script compiles RTTOV, optionally with external libraries (HDF5 or netCDF)."
echo "It should be run from the src/ directory."
echo
echo "Most RTTOV functionality is available without external libraries: only those"
echo "features listed below require an external library to compile."
echo
echo "If compiling with HDF5 (recommended) or netCDF you must first edit the"
echo "build/Makefile.local file with the location of the required library."
echo
echo "RTTOV features which have external dependencies:"
echo "Reading HDF5 coefficient files : requires the HDF5 library."
echo "Emissivity/BRDF atlases        : require either the HDF5 or netCDF library."
echo "Python/C/C++ interfaces        : require the emissivity/BRDF atlases."
echo "RTTOV GUI                      : requires the HDF5 library."
echo "Python interface and RTTOV GUI : f2py must be installed."
echo
echo '========================================================================='
echo

# Functions for interpreting input

yesno () {
    result="n"
    if [ $# -eq 1 ]; then
        if [ $1 = "y" ] || [ $1 = "Y" ]; then
            result="y"
        fi
    fi
    echo $result
}

onezero () {
    if [ $1 = "y" ]; then
        echo 1
    else
        echo 0
    fi
}


# Check we are in src/ directory; if not then try to find it

if [[ $(pwd) != *"src" || ! -d main/ ]]; then
    if [[ -d src/main/ ]]; then
        echo "Not in RTTOV src/ directory, changing to src/"
        cd src
    elif [[ -d ../src/main/ ]]; then
        echo "Not in RTTOV src/ directory, changing to src/"
        cd ../src
    elif [[ -d ../../src/main/ ]]; then
        echo "Not in RTTOV src/ directory, changing to src/"
        cd ../../src
    else
        echo "This script should be run from the RTTOV src/ directory."
        exit 1
    fi
    echo
fi


# Specify build flags and installation directory

echo "Compiler flag files available in build/arch/:"
echo
ls ../build/arch/
echo
echo "Specify required compiler flag file (leave blank for default: gfortran)"
read -p "> " myarch
if [[ $myarch = "" ]]; then
    myarch="gfortran"
elif [[ ! -f ../build/arch/$myarch ]]; then
    echo "Error: $myarch not found in build/arch/ directory."
    echo "Either specify an existing file or add your own using an existing file as a template."
    exit 1
fi

# changed by MS so that theres no path choice
installdir="./"

# Sort out external libraries: set options for regenerating Makefiles and specify the build target

target=
hdf5="n"
f2py="n"
gui="n"
echo "Have you updated the file build/Makefile.local with the location of your HDF5 or netCDF installation? (y/n)"
read -p "> " extlib
extlib=$(yesno $extlib)

if [[ $extlib = "y" ]]; then
    target=all

    echo "Are you compiling with the HDF5 library (y) or the netCDF library (n)? (y/n)"
    read -p "> " hdf5
    hdf5=$(yesno $hdf5)

    echo "Testing for f2py..."
    f2py -h 1> /dev/null 2> /dev/null
    if [[ $? -eq 0 ]]; then
        if [[ $hdf5 = "y" ]]; then
            echo "...f2py detected: do you want to compile the Python wrapper and RTTOV GUI? (y/n)"
        else
            echo "...f2py detected: do you want to compile the Python wrapper? (y/n)"
        fi
        read -p "> " f2py
        f2py=$(yesno $f2py)
        if [[ $hdf5 = "y" ]]; then
            gui=$f2py
        fi
    else
        echo "...f2py does not appear to be installed (f2py -h failed), Python interface and RTTOV GUI cannot be compiled"
    fi
fi


# Check for a previous build

clean="n"
if [[ -d ../$installdir/bin ]]; then
    echo
    echo "Previous build detected in $installdir: perform a clean compilation? (y/n)"
    echo "Choose y if you have changed the compilation options since the previous build and re-compilation fails."
    read -p "> " clean
    clean=$(yesno $clean)
    echo
fi


# Additional flags for make

makeflags=
echo
echo "Specify any additional flags to pass to make (e.g. -j); leave blank if unsure"
read -p "> " makeflags
echo


# Summarise the user inputs

echo
echo "==========================="
echo " RTTOV COMPILATION SUMMARY "
echo "==========================="
echo
echo "Compiling with flags    : $myarch"
echo "Compiling in directory  : $installdir"
echo
echo "RTTOV features available:"
echo "HDF5 coefficient I/O    : $hdf5"
echo "Emissivity/BRDF atlases : $extlib"
echo "C/C++ wrapper           : $extlib"
echo "Python wrapper          : $f2py"
echo "RTTOV GUI               : $gui"
echo
if [[ $extlib = "y" ]]; then
    if [[ $hdf5 = "y" ]]; then
        echo "You must use the **HDF5** format IR emissivity/BRDF atlas files."
    else
        echo "You must use the **netCDF** format IR emissivity/BRDF atlas files."
    fi
    echo "These are available from http://nwpsaf.eu/deliverables/rtm/rtm_rttov11.html"
    echo
fi

cmd_mkfile="../build/Makefile.PL RTTOV_HDF=$(onezero $hdf5) RTTOV_F2PY=$(onezero $f2py)"
cmd_clean="make ARCH=$myarch INSTALLDIR=$installdir clean $makeflags"
cmd_build="make ARCH=$myarch INSTALLDIR=$installdir $target $makeflags"

echo "Regenerating Makefiles using:"
echo "$ $cmd_mkfile"
echo
echo "Compiling RTTOV using:"
if [[ $clean = "y" ]]; then
    echo "$ $cmd_clean"
fi
echo "$ $cmd_build"

echo
echo "OK to continue and compile RTTOV? (y/n)"
read -p "> " yn
if [[ $(yesno $yn) = "n" ]]; then
    exit
fi


# Compile RTTOV

echo
echo "Regenerating Makefiles..."
$cmd_mkfile

echo
echo "Compiling RTTOV..."
if [[ $clean = "y" ]]; then
    $cmd_clean
fi
$cmd_build

if [[ $? -eq 0 ]]; then
    echo
    echo "RTTOV compiled successfully"
    echo
fi
