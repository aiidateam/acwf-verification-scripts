#!/bin/bash

PLUGIN_NAME=quantum_espresso

verdi archive create -G commonwf-oxides/set1/structures/${PLUGIN_NAME} commonwf-oxides/set1/workflows/${PLUGIN_NAME} -- commonwf-oxides_set1_results_${PLUGIN_NAME}.aiida


