#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)"

mkdir -p ${SCRIPT_DIR}/deps/src
cd ${SCRIPT_DIR}/deps/src

wget http://downloads.sourceforge.net/project/unittest-cpp/UnitTest%2B%2B/1.4/unittest-cpp-1.4.zip


