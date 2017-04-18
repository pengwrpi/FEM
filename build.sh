set -x
cd build
verbose=""
[ "$2" = "on" ] && verbose="VERBOSE=1"
make $verbose $1
cd -
set +x

