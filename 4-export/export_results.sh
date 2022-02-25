#!/bin/bash

SET_NAME='oxides-verification-PBE-v1'

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PLUGIN_NAME=`cat "$SCRIPT_DIR"/../plugin_name.txt 2>/dev/null`

if [ "$PLUGIN_NAME" == "" ]
then
    echo 'ERROR: You need to define a file `../plugin_name.txt`, containing the name of your plugin (siesta, quantum_espresso, ...) in the format expected by the aiida-common-workflows project'
    exit 1
fi

verdi archive create -G acwf-verification/${SET_NAME}/structures/${PLUGIN_NAME} acwf-verification/${SET_NAME}/workflows/${PLUGIN_NAME} -- acwf-verification_${SET_NAME}_results_${PLUGIN_NAME}.aiida


