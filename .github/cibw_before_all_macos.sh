#!/bin/bash
# Adapted from h5py. Licensed under the BSD 3-Clause license.
# Copyright (c) 2008 Andrew Collette and contributors
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:

# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the
#    distribution.

# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

set -eo pipefail # Option pour que le script exit si une erreur est détecté dans l'un des pipeline 
set +x # set -x correspond au mode debug qui permet de print les commende avant de les exécuter, set -x permet d'éviter ça (utile pour éviter que les token github apparaissent dans les log) 

function setup_github_env {
    echo "HDF5_ROOT=${HDF5_DIR}" | tee -a $GITHUB_ENV
    echo "HighFive_ROOT=${HIGHFIVE_DIR}" | tee -a $GITHUB_ENV
    echo "SUNDIALS_ROOT=${SUNDIALS_DIR}" | tee -a $GITHUB_ENV
    echo "MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET}" | tee -a $GITHUB_ENV
    echo "DYLD_FALLBACK_LIBRARY_PATH=${HDF5_DIR}/lib:${YAML_CPP_DIR}/lib" | tee -a $GITHUB_ENV
    echo "yaml-cpp_ROOT=${YAML_CPP_DIR}" | tee -a $GITHUB_ENV
}

# Check if the current script ($0) is called with the correct option ($1) which should correspond to the project path 
if [[ "$1" == "" ]] ; then
    echo "Usage: $0 <PROJECT_PATH>"
    exit 1
fi

PROJECT_PATH="$1"
ARCH=$(uname -m) # Récupère l'architecture du système (arm64 / x86_64 ...)
GENERATOR="Ninja"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

source "${SCRIPT_DIR}/dependencies.sh"

HDF5_DIR="${PROJECT_PATH}/cache/hdf5/${HDF5_VERSION}-${ARCH}"
HIGHFIVE_DIR="${PROJECT_PATH}/cache/highfive/${HIGHFIVE_VERSION}-${ARCH}"
SUNDIALS_DIR="${PROJECT_PATH}/cache/sundials/${SUNDIALS_VERSION}-${ARCH}"
SUNDIALS_BUILD_OPTIONS=(
    "-DENABLE_LAPACK=ON"
    "-DBLA_VENDOR=Apple"
    "-DSUNDIALS_LAPACK_CASE=LOWER"
    "-DSUNDIALS_LAPACK_UNDERSCORES=NONE"
)
YAML_CPP_DIR="${PROJECT_PATH}/cache/yaml-cpp/${YAML_CPP_VERSION}-${ARCH}"

# When compiling HDF5, we should use the minimum across all Python versions for a given
# arch, for versions see for example a more updated version of the following:
# https://github.com/pypa/cibuildwheel/blob/9c75ea15c2f31a77e6043b80b1b7081372319d85/cibuildwheel/macos.py#L302-L315
if [[ "${ARCH}" == "arm64" ]]; then
    export MACOSX_DEPLOYMENT_TARGET="11.0"
else
    # This is the minimum version for Cantera
    export MACOSX_DEPLOYMENT_TARGET="10.15"
fi

lib_name=libhdf5.dylib
inc_name=highfive.hpp

if [ -f ${HDF5_DIR}/lib/${lib_name} ] && [ -f ${HIGHFIVE_DIR}/include/highfive/${inc_name} ]; then
    echo "using cached build"
    setup_github_env
    exit 0
else
    echo "building dependencies"
fi

brew install ninja cmake --formula

source "${SCRIPT_DIR}/build_dependencies.sh" # Build here hdf5, highFive, sundials, yaml_cpp

setup_github_env