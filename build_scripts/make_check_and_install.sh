#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)"
source ${SCRIPT_DIR}/setup_env_vars.sh


for BUILD_TYPE in debug release; do
    mkdir -p ${SCRIPT_DIR}/build/${BUILD_TYPE}
    cd ${SCRIPT_DIR}/build/${BUILD_TYPE}
    cmake -DCMAKE_INSTALL_PREFIX=${COELA_INSTALL} ${SCRIPT_DIR}
    make  check -j$NCPU
    make install -j$NCPU
done
