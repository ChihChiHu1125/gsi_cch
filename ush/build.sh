#!/bin/bash

set -eu

while getopts "r" option;
do
 case $option in
  r)
   echo "Recieved -r flag, will recompile without clean"
   export BUILD_CLEAN="NO"
   ;;
  :)
   echo "option -$OPTARG needs an argument"
   ;;
  *)
   echo "invalid option -$OPTARG, exiting..."
   exit
   ;;
 esac
done

# Get the root of the cloned GSI directory
readonly DIR_ROOT=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )/.." && pwd -P)

# User Options
BUILD_TYPE=${BUILD_TYPE:-"Release"}
CMAKE_OPTS=${CMAKE_OPTS:-}
COMPILER=${COMPILER:-"intel"}
BUILD_DIR=${BUILD_DIR:-"${DIR_ROOT}/build"}
INSTALL_PREFIX=${INSTALL_PREFIX:-"${DIR_ROOT}/install"}
GSI_MODE=${GSI_MODE:-"Regional"}  # By default build Regional GSI (for regression testing)
ENKF_MODE=${ENKF_MODE:-"GFS"}     # By default build Global EnKF  (for regression testing)
REGRESSION_TESTS=${REGRESSION_TESTS:-"YES"} # Build regression test suite

#==============================================================================#

# Detect machine (sets MACHINE_ID)
source $DIR_ROOT/ush/detect_machine.sh

# Load modules
#set -x
source $DIR_ROOT/ush/module-setup.sh
module use $DIR_ROOT/modulefiles
module load gsi_$MACHINE_ID
module list
module unload ncdiag
module unload netcdf
module unload hdf5
module load netcdf/4.7.4
module load hdf5/1.10.6
module load ncdiag
module list

# Set CONTROLPATH variables for Regression testing on supported MACHINE_ID
if [[ $MACHINE_ID = wcoss ]] ; then
    CONTROLPATH="/da/save/Michael.Lueken/svn1/build"
elif [[ $MACHINE_ID = wcoss_dell_p3 ]] ; then
    CONTROLPATH="/gpfs/dell2/emc/modeling/noscrub/Michael.Lueken/svn1/install/bin"
elif [[ $MACHINE_ID = hera.intel ]] ; then
    CONTROLPATH="/scratch1/NCEPDEV/da/Michael.Lueken/svn1/install/bin"
    #export CRTM_LIB=/scratch2/GFDL/gfdlscr/Mingjing.Tong/CRTM/REL-2.4.0/crtm_v2.4.0/lib/libcrtm.a
    #export CRTM_INC=/scratch2/GFDL/gfdlscr/Mingjing.Tong/CRTM/REL-2.4.0/crtm_v2.4.0/include
    #export CRTM_LIB=/scratch2/GFDL/gfdlscr/Mingjing.Tong/CRTM/crtm_im/crtm_im/lib64/libcrtm_static.a
    #export CRTM_INC=/scratch2/GFDL/gfdlscr/Mingjing.Tong/CRTM/crtm_im/crtm_im/module/crtm/Intel/19.1.2.20200623
fi

# Collect BUILD Options
CMAKE_OPTS+=" -DCMAKE_BUILD_TYPE=$BUILD_TYPE"

# Install destination for built executables, libraries, CMake Package config
CMAKE_OPTS+=" -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX"

# Configure for GSI and EnKF
CMAKE_OPTS+=" -DGSI_MODE=$GSI_MODE -DENKF_MODE=${ENKF_MODE}"

# Build regression test suite (on supported MACHINE_ID where CONTROLPATH exists)
[[ ${REGRESSION_TESTS} =~ [yYtT] ]] && CMAKE_OPTS+=" -DBUILD_REG_TESTING=ON -DCONTROLPATH=${CONTROLPATH:-}"

# Re-use or create a new BUILD_DIR (Default: create new BUILD_DIR)
if [[ ${BUILD_CLEAN:-"YES"} =~ [yYtT] ]] ; then
   rm -rf $BUILD_DIR
   mkdir -p $BUILD_DIR && cd $BUILD_DIR
fi

# Configure, build, install
set -x
if [[ ${BUILD_CLEAN:-"YES"} =~ [yYtT] ]] ; then
   cmake $CMAKE_OPTS $DIR_ROOT
else
   cd $BUILD_DIR
fi
make -j ${BUILD_JOBS:-8} VERBOSE=${BUILD_VERBOSE:-}
make install
set +x

exit
