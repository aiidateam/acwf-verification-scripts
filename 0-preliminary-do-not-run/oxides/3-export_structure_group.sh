#!/bin/bash

SET_NAME=set2

verdi archive create -G commonwf-oxides/${SET_NAME}/structures -- commonwf-oxides_${SET_NAME}_structures.aiida

## To check the nodes that were imported: you can e.g. use:
##
# for node in Group.get(label='commonwf-oxides/set1/structures').nodes:
#     print(node.extras)
