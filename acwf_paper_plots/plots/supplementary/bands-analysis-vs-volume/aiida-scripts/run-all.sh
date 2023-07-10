#!/bin/bash

for s in Cs2O CsO2 RbO3 ; do ./1-process-bands.py $s ; ./2-plot.py $s ; done
