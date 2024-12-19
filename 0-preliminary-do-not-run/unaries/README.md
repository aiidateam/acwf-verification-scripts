# Data generation

## New datasets for PBE Z=97-103 and Z=1-103 for PBEsol and LDA
Based on the minimum volumes obtained from EOS calculations by three different AE codes, i.e. `WIEN2k`, `FLEUR` and `SIRIUS+CP2k`, we have generated new datasets for the PBEsol and LDA functionals for `Z=1-103`. Moreover, we extended the previous PBE data by adding the elements `Z=97-103` (which can be accessed in the corresponding `PBE-missing-actinides` files).
We took the average of the three minimum volumes and converted the average central volumes into the corresponding lattice constant to create the `xsf` files. 
Due to technical issues, the `SIRIUS+CP2k` results for `Lr` are missing. Hence, the average for `Lr` was calculated using the results from `WIEN2k` and `FLEUR`.