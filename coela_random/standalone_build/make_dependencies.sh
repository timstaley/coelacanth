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
tar -xvf rngstreams-1.0.1.tar.gz
tar -xvf unuran-1.8.1.tar.gz


#Unittest++
cd ${SCRIPT_DIR}/deps/src/UnitTest++ 
make -j$NCPU
cp -v libUnitTest++.a ${COELA_INSTALL}/lib
cp -r src ${COELA_INSTALL}/include/UnitTest++

#RNGstreams
cd ${SCRIPT_DIR}/deps/src/rngstreams-1.0.1
./configure --enable-shared=yes --with-pic --prefix=${COELA_INSTALL}
make install -j$NCPU

#UNURan
cd ${SCRIPT_DIR}/deps/src/unuran-1.8.1
./configure --with-urng-rngstream --enable-shared=yes --prefix=${COELA_INSTALL} --with-pic --with-urng-default=rngstream
make install -j$NCPU

