#!/bin/bash

## Sets up the default install / include / lib paths.
## Can be sourced to allow manual testing of build commands.

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)"


if [ -z "$COELA_INSTALL" ]; then
	COELA_INSTALL=${SCRIPT_DIR}/install
fi

if [ -z "$NCPU" ]; then
	NCPU=4
fi


export C_INCLUDE_PATH=${C_INCLUDE_PATH:+${C_INCLUDE_PATH}:}${COELA_INSTALL}/include
export CPLUS_INCLUDE_PATH=${CPLUS_INCLUDE_PATH:+${CPLUS_INCLUDE_PATH}:}${COELA_INSTALL}/include
export LIBRARY_PATH=${LIBRARY_PATH:+${LIBRARY_PATH}:}${COELA_INSTALL}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:+${LD_LIBRARY_PATH}:}${COELA_INSTALL}/lib
