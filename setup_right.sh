export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/home/wei/Downloads/installer/scorec/install
set -x
mkdir -p build
cd build
flags='-g '
cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CXX_FLAGS=$flags \
    -DSCOREC_PREFIX="/home/wei/Downloads/installer/scorec/install"  ../
cd -
set +x
