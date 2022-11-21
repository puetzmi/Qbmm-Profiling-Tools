# This script sets the environment to compile and run applications using the GNU compilers. If the environment variable 'MKLROOT' is not set to the MKL directory, either the variable 'INTEL_ONEAPI_ROOT' must be set to the Intel(R) OneAPI directory (where 'mkl' is located) or Intel(R) OneAPI must be installed in the default location '/opt/intel/oneapi'.
#
#!/usr/bin/env bash

INTEL_ONEAPI_DEFAULT=/opt/intel/oneapi

# Find Intel(R) OneAPI directory if MKLROOT is not set
if [[ -z "${MKLROOT}" ]]; then
    if [[ -z "${INTEL_ONEAPI_ROOT}" ]]; then
        if [[ -d $INTEL_ONEAPI_DEFAULT ]]; then
            export INTEL_ONEAPI_ROOT=$INTEL_ONEAPI_DEFAULT
        else
            >&2 echo "MKL not found. Set one of the environment variables 'MKLROOT' or 'INTEL_ONEAPI_ROOT' to the correct location."
            return
        fi
    fi
    export MKLROOT=$INTEL_ONEAPI_ROOT/mkl/latest
fi

# Set environment
source $MKLROOT/env/vars.sh
export CC=gcc
export CXX=g++
