# This script sets the environment to compile and run applications using the Intel(R) OneAPI compilers. If Intel(R) OneAPI is not installed in the default location '/opt/intel/oneapi' the environment variable 'INTEL_ONEAPI_ROOT' must be set to the right directory.
#
#!/usr/bin/env bash

INTEL_ONEAPI_DEFAULT=/opt/intel/oneapi

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
source $INTEL_ONEAPI_ROOT/mkl/latest/env/vars.sh
source $INTEL_ONEAPI_ROOT/compiler/latest/env/vars.sh
source $INTEL_ONEAPI_ROOT/mpi/latest/env/vars.sh
export CC=icx
export CXX=icpx
