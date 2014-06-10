#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)"

bash ${SCRIPT_DIR}/get_dependencies.sh
bash ${SCRIPT_DIR}/make_dependencies.sh
bash ${SCRIPT_DIR}/make_check_and_install.sh
