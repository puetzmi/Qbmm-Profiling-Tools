# This script rebuilds all build-subdirectories in the given directory by running the script 'scripts/build.sh' on each one of them, provided that the compiler type 'intel' or 'gnu' is contained in the name of the build subdirectory.

#!/bin/env bash
BUILD_DIR=$1
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
CWD=$PWD
cd $BUILD_DIR
for build_dir in */; do
    env -i HOME=${HOME} bash -l -c "$SCRIPT_DIR/build.sh $build_dir"
done
cd $CWD
