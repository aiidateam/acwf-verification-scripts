This script compare codes result data get from same pseudopotential 
and creates an histogram of the differences
in V0,B0 and B1.

A command argument can be used to create a single
histogram with both oxides and unaries or to
create a separete histogram for the two sets.

There are three codes using PseudoDojo-v04:
- abinit
- castep
- quantum_espresso

There are two codes using SSSP-v1.2:
- quantum_espresso
- cp2k-sirus

The script can be run with the following command
The last argument is the comment string appended to the output file name.

```bash
./generate_paper_histos_all_materials.py all "SIRIUS/CP2K@PW|SSSP-prec-v1.2" "Quantum ESPRESSO@PW|SSSP-prec-v1.2" sssp-prec
```

```bash
./generate_paper_histos_all_materials.py all "ABINIT@PW|PseudoDojo-v0.4" "CASTEP@PW|PseudoDojo-v0.4" dojo
./generate_paper_histos_all_materials.py all "ABINIT@PW|PseudoDojo-v0.4" "Quantum ESPRESSO@PW|PseudoDojo-v0.4" dojo
./generate_paper_histos_all_materials.py all "CASTEP@PW|PseudoDojo-v0.4" "Quantum ESPRESSO@PW|PseudoDojo-v0.4" dojo
```