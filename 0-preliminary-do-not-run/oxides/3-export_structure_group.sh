#!/bin/bash

SET_NAME=oxides-verification-PBE-v1

verdi archive create -G acwf-verification/${SET_NAME}/structures -- acwf-verification_${SET_NAME}_structures.aiida

## To check the nodes that were imported: you can e.g. use:
##
# for node in Group.get(label='acwf-verification/oxides-verification-PBE-v1/structures').nodes:
#     print(node.extras)
