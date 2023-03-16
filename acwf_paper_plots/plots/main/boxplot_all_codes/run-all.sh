#!/bin/bash
./plot_box_all.py together up-to-Bi-no-lanthanides
./plot_box_all.py together only-lanthanides "WIEN2k@(L)APW+lo+LO" "FLEUR@LAPW+LO" "CASTEP@PW|C19MK2" "VASP@PW|GW-PAW54*" "Quantum ESPRESSO@PW|SSSP-prec-1.2"
./plot_box_all.py together only-actinides "WIEN2k@(L)APW+lo+LO" "FLEUR@LAPW+LO" "CASTEP@PW|C19MK2" "VASP@PW|GW-PAW54*"
