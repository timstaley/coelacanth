#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)"

mkdir -p ${SCRIPT_DIR}/deps/src
cd ${SCRIPT_DIR}/deps/src

wget http://downloads.sourceforge.net/project/unittest-cpp/UnitTest%2B%2B/1.4/unittest-cpp-1.4.zip

wget http://www.cern.ch/mathlibs/sw/5_28_00/Minuit2/Minuit2-5.28.00.tar.gz

wget http://statmath.wu.ac.at/software/RngStreams/rngstreams-1.0.1.tar.gz
wget http://statmath.wu.ac.at/unuran/unuran-1.8.1.tar.gz
