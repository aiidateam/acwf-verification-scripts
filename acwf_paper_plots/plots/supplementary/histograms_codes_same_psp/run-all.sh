#!/bin/bash
./generate_paper_histos_all_materials.py all "SIRIUS/CP2K@PW|SSSP-prec-v1.2" "Quantum ESPRESSO@PW|SSSP-prec-v1.2" sssp-prec
./generate_paper_histos_all_materials.py all "ABINIT@PW|PseudoDojo-v0.4" "CASTEP@PW|PseudoDojo-v0.4" dojo-abinit-castep
./generate_paper_histos_all_materials.py all "ABINIT@PW|PseudoDojo-v0.4" "Quantum ESPRESSO@PW|PseudoDojo-v0.4" dojo-abinit-qe
./generate_paper_histos_all_materials.py all "CASTEP@PW|PseudoDojo-v0.4" "Quantum ESPRESSO@PW|PseudoDojo-v0.4" dojo-castep-qe