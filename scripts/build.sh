# This script checks the given directory name for the compiler type 'intel' or 'gnu', sources the corresponding script to load environment variables and builds the project with the current configuration. The script is run by 'scripts/build_all.sh' and should not be executed by the user.

TARGET=$1
if [[ $TARGET == *"intel"* ]]; then
    source $(dirname $0)/load-intel.sh
elif [[ $TARGET == *"gnu"* ]]; then
    source $(dirname $0)/load-gnu.sh
fi
cd $TARGET
NPROC=$(nproc)
make -j$NPROC
cd - > /dev/null
