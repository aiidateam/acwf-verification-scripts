#!/bin/bash
./generate_paper_periodic_tables_same_pseudos.py "SIRIUS/CP2K@PW|SSSP-prec-v1.2" "Quantum ESPRESSO@PW|SSSP-prec-v1.2" sssp-prec
./generate_paper_periodic_tables_same_pseudos.py "ABINIT@PW|PseudoDojo-v0.4" "CASTEP@PW|PseudoDojo-v0.4" dojo
./generate_paper_periodic_tables_same_pseudos.py "ABINIT@PW|PseudoDojo-v0.4" "Quantum ESPRESSO@PW|PseudoDojo-v0.4" dojo
./generate_paper_periodic_tables_same_pseudos.py "CASTEP@PW|PseudoDojo-v0.4" "Quantum ESPRESSO@PW|PseudoDojo-v0.4" dojo