#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)"

if [ -z "$COELA_INSTALL" ]; then
	COELA_INSTALL=${SCRIPT_DIR}/install
fi

if [ -z "$NCPU" ]; then
	NCPU=1
fi

export C_INCLUDE_PATH=${C_INCLUDE_PATH:+${C_INCLUDE_PATH}:}${COELA_INSTALL}/include
export CPLUS_INCLUDE_PATH=${CPLUS_INCLUDE_PATH:+${CPLUS_INCLUDE_PATH}:}${COELA_INSTALL}/include
export LIBRARY_PATH=${LIBRARY_PATH:+${LIBRARY_PATH}:}${COELA_INSTALL}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:+${LD_LIBRARY_PATH}:}${COELA_INSTALL}/lib

mkdir -p ${SCRIPT_DIR}/build

echo $LIBRARY_PATH

cd ${SCRIPT_DIR}/build
cmake -DCMAKE_INSTALL_PREFIX=${COELA_INSTALL} ${SCRIPT_DIR}
make check -j$NCPU
make install -j$NCPU

