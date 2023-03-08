#!/bin/bash
./plot_box_all.py together up-to-Bi-no-lanthanides
./plot_box_all.py together only-lanthanides "WIEN2k" "FLEUR" "CASTEP+..." "VASP+..." "Quantum ESPRESSO+SSSP"
./plot_box_all.py together only-actinides "WIEN2k" "FLEUR" "CASTEP+..." "VASP+..."
