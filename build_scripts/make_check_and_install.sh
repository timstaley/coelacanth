#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)"
source ${SCRIPT_DIR}/setup_env_vars.sh

mkdir -p ${SCRIPT_DIR}/build
cd ${SCRIPT_DIR}/build
cmake -DCMAKE_INSTALL_PREFIX=${COELA_INSTALL} ${SCRIPT_DIR}
make  check -j$NCPU
make install -j$NCPU

