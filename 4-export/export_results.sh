#!/bin/bash

CODE_NAME=quantum_espresso

verdi archive create -G commonwf-oxides/set1/structures/${CODE_NAME} commonwf-oxides/set1/workflows/${CODE_NAME} -- commonwf-oxides_set1_results_${CODE_NAME}.aiida


