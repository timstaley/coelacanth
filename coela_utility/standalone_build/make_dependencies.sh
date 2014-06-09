#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)"

if [ -z "$COELA_DEPS" ]; then
	COELA_INSTALL=${SCRIPT_DIR}/install
fi

if [ -z "$NCPU" ]; then
	NCPU=1
fi

mkdir -p ${COELA_INSTALL}/lib
mkdir -p ${COELA_INSTALL}/include
export C_INCLUDE_PATH=${C_INCLUDE_PATH:+${C_INCLUDE_PATH}:}${COELA_INSTALL}/include
export CPLUS_INCLUDE_PATH=${CPLUS_INCLUDE_PATH:+${CPLUS_INCLUDE_PATH}:}${COELA_INSTALL}/include
export LIBRARY_PATH=${LIBRARY_PATH:+${LIBRARY_PATH}:}${COELA_INSTALL}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:+${LD_LIBRARY_PATH}:}${COELA_INSTALL}/lib


cd ${SCRIPT_DIR}/deps/src/
unzip unittest-cpp-1.4.zip

#Unittest++
cd ${SCRIPT_DIR}/deps/src/UnitTest++ 
make -j$NCPU
cp -v libUnitTest++.a ${COELA_INSTALL}/lib
cp -r src ${COELA_INSTALL}/include/UnitTest++

