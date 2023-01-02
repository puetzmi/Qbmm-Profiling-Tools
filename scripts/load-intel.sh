# This script sets the environment to compile and run applications using the Intel(R) OneAPI compilers. If Intel(R) OneAPI is not installed in the default location '/opt/intel/oneapi' the environment variable 'INTEL_ONEAPI_ROOT' must be set to the right directory.
#
#!/usr/bin/env bash

# Source configuration file
SCRIPTPATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
source $(dirname $SCRIPTPATH)/intel-env-config.sh

# Find Intel(R) OneAPI directory
if [[ -z "${INTEL_ONEAPI_ROOT}" ]]; then
    if [[ -d $INTEL_ONEAPI_DEFAULT ]]; then
        export INTEL_ONEAPI_ROOT=$INTEL_ONEAPI_DEFAULT
    else
        >&2 echo "Intel OneAPI not found. Set the environment variable 'INTEL_ONEAPI_ROOT' to the correct location."
        return
    fi
fi

# Set environment
source $INTEL_ONEAPI_ROOT/mkl/$INTEL_MKL_VERSION/env/vars.sh
source $INTEL_ONEAPI_ROOT/compiler/$INTEL_COMPILER_VERSION/env/vars.sh
source $INTEL_ONEAPI_ROOT/mpi/$INTEL_MPI_VERSION/env/vars.sh
export CC=icx
export CXX=icpx
