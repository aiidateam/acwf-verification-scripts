# Verification of magnetic systems

After several tests, we decided to use the energy as a function of fixed total magnetization to verify the magnetic systems. 

## Summary of the procedure

* **The volume remains fixed at the magnetic minimum volume previously determined.**
* We evaluate the energy at 9 (symmetrically distributed) fixed total magnetization values in steps of 0.03 mu_B/atom around the magnetization corresponding to the magnetic minimum volume previously determined. 
    * In some cases, when the magnetization reaches an integer value that would correspond to a full shell configuration, we adjust the number of points to avoid that value and only consider values on one side of that integer value, at least 5 datapoints.

## Inputs

* We provide a set of structures with the magnetic minimum volume in `.xsf` format and a corresponding `AiiDA` archive. You can find them in the usual locations in the `0-preliminary-do-not-run/unaries/` and `1-preliminary/magnetic-verification/` folders of this repository.
* The total magnetizations (in mu_B/cell) to be checked for each structure are provided in the input JSON file `fixed_total_magnetizations-unaries-{SET_NAME}.json`. We explicitly define them as the number of points might change depending on the configuration

As some configurations are still refined, we start the verification for a subset. This list will be updated:
* `PBE-magnetic-d-block-OFClBr`: Only d-block elements and `O, F, Cl, Br`. Note, only magnetic configurations are included.

## Outputs

The expected format of the outputs is defined in the `expected-output-format-magnetic-verification-v001.json` file. Since this is still part of the development, this output format might be extended in the future.
* The output contains te following keys:
    * `script_version`: version of the script that generated the output
    * `set_name`: name of the set of structures. This is important since we only start with a subset that might be extended in the future.
    * `total_magnetization_energy_data`: a list of lists, where the first element of each sub-list contains the fixed total magnetization (in mu_B/cell) and the second element corresponding total energy (in eV/cell) for each configuration.
    * `fermi_energies`: a list of lists, with the same order as `total_magnetization_energy_data`, containing the Fermi energies (in eV) for each magnetization value. The first element of each sub-list corresponds to the `Fermi` energy of the spin-down channel, the second element to the spin-up channel.
    * `density_matrices`: a list of dictionaries, with the same order as `total_magnetization_energy_data`, containing the density matrices for each magnetization value. The dictionary has two keys: `spin-up` and `spin-down`. Each key contains the actual density matrices. For d-block elements, these are `5 x 5 x 2` list. The two elements of each inner list correspond to the real and imaginary part. The density matrices are given with respect to the complex spherical harmonics. Therefore, the `5` columns correspond to `-m, ..., +m`. \
    In case of f-block elements, the format will be `7 x 7 x 2`.