source ~/.bashrc
cantera_gcc

BASE_DIR=/home/alberti3/WORKSPACE/cantera/plasma_cantera
INSTALL_DIR=$BASE_DIR/cantera_lib-with-$COMPILER-$COMPILER_VERSION_MICRO

scons clean

scons build optimize='False' prefix=$INSTALL_DIR \
#python_package = 'minimal' \
matlab_toolbox = 'n' \
f90_interface = 'y' \
use_sundials = 'n'
#blas_lapack_libs = 'blas,lapack'
#blas_lapack_dir = '/usr/lib'

scons install
