#!/bin/bash

# =================================================================

# Script to download RTTOV coefficient files from the website:
# run the script from within your rtcoef_rttov11/ directory.

# The function "download_extract" downloads all files with a given
# file extension from the specified directory. It then unzips and
# untars the files as necessary. Intermediate (e.g. .tar) files
# are deleted so only coefficient files are left.

# It uses the local directory structure of rtcoef_rttov11/ to
# determine the rttov*pred*L/ directories on the server.

# Filetypes are distinguished as follows:
# - HDF5 files end in .H5
# - tarballs of smaller coef files end in .tar.bz2
# - other large inidividually zipped files end in .dat.bz2

# =================================================================


echo '========================================================================='
echo ' RTTOV coefficient download script'
echo '========================================================================='
echo
echo '* Individual coefficient files may be downloaded from the RTTOV website:'
echo '     http://nwpsaf.eu/deliverables/rtm/rttov11_coefficients.html'
echo '  See the website for information about the different types of coefficient files.'
echo
echo '* This script can be used to download all coefficient files or all' \
     'files of a particular type.'
echo '  You must call this from within the rtcoef_rttov11 directory:'
echo '  the coefficient files are unzipped and extracted in the standard' \
     'locations within this directory.'
echo
echo '* Be aware: EXISTING COEFFICIENT FILES MAY BE OVERWRITTEN WITHOUT WARNING.'
echo
echo '* You only need to download the coefficient files for the instruments you'
echo '  are simulating and the kinds of simulations you are carrying out.'
echo
echo '* For larger files HDF5 format is recommended over ASCII format.'
echo
echo '* The VIS/IR/MW optical depth coefficient files are included with the RTTOV distribution.'
echo '  However the tarballs on the website may contain updated coefficient files.'
echo
echo 'File types are:'
echo '  - Optical depth (rtcoef) coefficient files for MW, IR and VIS/IR'
echo '  - IR aerosol (scaer) and cloud (sccld) coefficient files'
echo '  - RTTOV-SCATT MW mietable files'
echo '  - HDF5/ASCII optical depth (rtcoef) files for hi-res IR sounders'
echo '  - HDF5/ASCII hi-res IR sounder aerosol (scaer) and cloud (sccld) coefficient files'
echo '  - HDF5/ASCII PC-RTTOV coefficient files'
echo
echo '========================================================================='
echo


protocol='http://'
server='nwpsaf.eu/'
path='downloads/rtcoef_rttov11/'

yesno () {
    result=0
    if [ $# -eq 1 ]; then
        if [ $1 = "y" ] || [ $1 = "Y" ]; then
            result=1
        fi
    fi
    echo $result
}

download_extract () {
    # Downloads all files from a directory on server with given file extension
    # Zipped files are unzipped, tarballs are extracted and .tar files deleted
    # Tarball extraction overwrites existing files
    # Other files will be renamed rather than being overwritten

    # First argument is file extension e.g. ".tar.bz2" or ".dat.bz2" or ".H5"
    # Second argument is directory e.g. "cldaer" or "rttov7pred54L"

    d=$2
    echo
    echo "Downloading files with extension $1 from ${protocol}${server}${path}${d}, this may take some time..."
    rm -fr $server
    wget -q -r -np -l1 -A$1 ${protocol}${server}${path}${d}
    if [[ -d ${server}${path}${d} ]]; then
        if [[ $(ls ${server}${path}${d}) ]]; then
            for f in $(ls ${server}${path}${d}); do
                echo "File found: $(basename $f)"
                if [[ -f ${d}/$(basename $f) ]]; then
                    mv ${d}/$(basename $f) ${d}/$(basename $f)_old_$(date +%Y%m%d_%H%M%S)
                    echo "Renamed existing file ${d}/$(basename $f)"
                fi
            done
            mv ${server}${path}${d}/* $d
        else
            echo "No files found with names: ${protocol}${server}${path}${d}/*$1"
        fi
    else
        echo "Directory not found on server: ${protocol}${server}${path}${d}"
    fi
    rm -fr $server

    cd $d
    for f in $(ls); do
        if [[ $f = *".bz2" ]]; then
            echo "Unzipping $f..."
            bunzip2 $f 2> /dev/null
            local result=$?
            if [ $result -eq 1 ]; then
                echo "Unzipped file is already present, not unzipping $f"
            elif [ $result -ne 0 ]; then
                echo "Error unzipping, not unzipping $f"
            fi
        fi
    done

    for f in $(ls); do
        if [[ $f = *".tar" ]]; then
            echo "Extracting tarball $f..."
            tar xf $f 2> /dev/null
            local result=$?
            if [ $result -ne 0 ]; then
                echo "Tar extraction failed, not extracting $f"
            fi
            rm -f $f
        fi
    done
    cd ../
}


# Check we are in rtcoef_rttov11/ directory; if not then try to find it

if [[ $(pwd) != *"rtcoef_rttov11" ]]; then
    if [[ -d ../rtcoef_rttov11 ]]; then
        echo "Not in RTTOV rtcoef_rttov11/ directory, changing to rtcoef_rttov11/"
        cd ../rtcoef_rttov11
    elif [[ -d ../../rtcoef_rttov11 ]]; then
        echo "Not in RTTOV rtcoef_rttov11/ directory, changing to rtcoef_rttov11/"
        cd ../../rtcoef_rttov11
    else
        echo "This script must be run from the RTTOV rtcoef_rttov11/ directory."
        exit 1
    fi
    echo
fi

echo "For the larger coefficient files do you want HDF5 (y) or ASCII (n) format? "
read -p "> " yn
hdf=0
ascii=1
if [ $(yesno $yn) -eq 1 ]; then
    hdf=1
    ascii=0
fi

echo "Download all files (y) or specify files to download (n)? "
read -p "> " yn
if [ $(yesno $yn) -eq 1 ]; then
    get_visirmw=1
    get_ircldaer=1
    get_mietables=1
    get_hiresascii=$ascii
    get_hireshdf5=$hdf
    get_hirescldaerascii=$ascii
    get_hirescldaerhdf5=$hdf
    get_pcascii=$ascii
    get_pchdf5=$hdf
else
    echo "Download VIS/IR/MW rtcoef files? (y/n) "
    read -p "> " yn
    get_visirmw=$(yesno $yn)

    echo "Download IR cld/aer coef files? (y/n) "
    read -p "> " yn
    get_ircldaer=$(yesno $yn)

    echo "Download MW mietable files? (y/n) "
    read -p "> " yn
    get_mietables=$(yesno $yn)

    echo "Download hi-res IR rtcoef files? (y/n) "
    read -p "> " yn
    get_hiresascii=0
    get_hireshdf5=0
    if [[ $(yesno $yn) -eq 1 ]]; then
        get_hiresascii=$ascii
        get_hireshdf5=$hdf
    fi

    echo "Download hi-res IR cld/aer coef files? (y/n) "
    read -p "> " yn
    get_hirescldaerascii=0
    get_hirescldaerhdf5=0
    if [[ $(yesno $yn) -eq 1 ]]; then
        get_hirescldaerascii=$ascii
        get_hirescldaerhdf5=$hdf
    fi

    echo "Download PC coef files? (y/n) "
    read -p "> " yn
    get_pcascii=0
    get_pchdf5=0
    if [[ $(yesno $yn) -eq 1 ]]; then
        get_pcascii=$ascii
        get_pchdf5=$hdf
    fi
fi
echo


# Download VIS/IR/MW coef files:
if [ $get_visirmw -eq 1 ]; then
    for d in $(ls); do
        if [[ -d $d && $d = "rttov"*"pred"*"L" ]]; then
            download_extract ".tar.bz2" $d
        fi
    done
fi

# Download non-hires IR cld/aer files:
if [ $get_ircldaer -eq 1 ]; then
    download_extract ".tar.bz2" "cldaer"
fi

# Download Mietable files:
if [ $get_mietables -eq 1 ]; then
    download_extract ".dat.bz2" "mietable"
fi

# Download ASCII hi-res IR sounder files:
if [ $get_hiresascii -eq 1 ]; then
    for d in $(ls); do
        if [[ -d $d && $d = "rttov"*"pred"*"L" ]]; then
            download_extract ".dat.bz2" $d
        fi
    done
fi

# Download HDF5 hi-res IR sounder files:
if [ $get_hireshdf5 -eq 1 ]; then
    for d in $(ls); do
        if [[ -d $d && $d = "rttov"*"pred"*"L" ]]; then
            download_extract ".H5" $d
        fi
    done
fi

# Download ASCII hires IR cld/aer files:
if [ $get_hirescldaerascii -eq 1 ]; then
    download_extract ".dat.bz2" "cldaer"
fi

# Download HDF5 hires IR cld/aer files:
if [ $get_hirescldaerhdf5 -eq 1 ]; then
    download_extract ".H5" "cldaer"
fi

# Download ASCII PC coef files:
if [ $get_pcascii -eq 1 ]; then
    download_extract ".dat.bz2" "pc"
fi

# Download HDF5 PC coef files:
if [ $get_pchdf5 -eq 1 ]; then
    download_extract ".H5" "pc"
fi

echo
echo 'Please visit http://nwpsaf.eu/deliverables/rtm/rttov11_coefficients.html for information about the coefficient files'
echo