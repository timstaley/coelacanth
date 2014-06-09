#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)"

export COELA_INSTALL=${SCRIPT_DIR}/install
export NCPU=5

bash get_dependencies.sh
bash make_dependencies.sh
bash make_standalone_lib.sh
