#!/bin/bash

SET_NAME='set2'

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PLUGIN_NAME=`cat "$SCRIPT_DIR"/../plugin_name.txt 2>/dev/null`

if [ "$PLUGIN_NAME" == "" ]
then
    echo 'ERROR: You need to define a file `../plugin_name.txt`, containing the name of your plugin (siesta, quantum_espresso, ...) in the format expected by the aiida-common-workflows project'
    exit 1
fi

verdi archive create -G commonwf-oxides/${SET_NAME}/structures/${PLUGIN_NAME} commonwf-oxides/${SET_NAME}/workflows/${PLUGIN_NAME} -- commonwf-oxides_${SET_NAME}_results_${PLUGIN_NAME}.aiida


