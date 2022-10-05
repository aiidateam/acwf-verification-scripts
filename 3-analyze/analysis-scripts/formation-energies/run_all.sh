#!/bin/bash

REF_PLUGIN=fleur

set -ev

for plugin in abinit-JTH-1.1-PBE abinit-PseudoDojo-0.4-PBE-SR-standard-psp8 abinit-PseudoDojo-0.5b1-PBE-SR-standard-psp8 cp2k_DZVP gpaw mk2-castep quantum_espresso vasp_recpot vasp_recpot_e700 wien2k-dk_0.06 
do  
    ./compute_formation_energies.py $plugin $REF_PLUGIN
    ./plot_histo_formation_energies.py $plugin $REF_PLUGIN formation-energy -0.05
done
