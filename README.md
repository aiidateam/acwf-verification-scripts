# acwf-verification-scripts
Scripts to facilitate the [AiiDA common workflows (ACWF)](https://github.com/aiidateam/aiida-common-workflows) sub-project focusing on verification of codes (on a set of unaries and oxides for all chemical elements between Z=1 and Z=96).

These scripts were used to manage all simulations, analyze the data and produce the plots and tables of the paper:

>  [E. Bosoni et al., *How to verify the precision of density-functional-theory implementations via reproducible and universal workflows*, **Nat. Rev. Phys. 6**, 45 (2024)](https://doi.org/10.1038/s42254-023-00655-3)

If you use these scripts, we would be grateful if you could cite the paper above.

Results from the paper above (whose raw files are available the subfolder `acwf_paper_plots`) can be visualized directly in your browser on this Materials Cloud Discover section: [https://acwf-verification.materialscloud.org](https://acwf-verification.materialscloud.org).

## Content of the folders

These scripts are a template (based on Quantum ESPRESSO) that can be adapted relatively simply for any other code.

- `0-preliminary-do-not-run`: ignore this folder - these are helper scripts to create the import the structures from CIF into AiiDA nodes and create initial AiiDA archive file. The scripts will be run only once by one of us.

- `1-preliminary`: scripts to import the group with the initial structures, and to create a subgroup containing only the nodes you actually want to run.

- `2-submit`: scripts to submit (in batches) the common EOS workflows for all systems

- `3-analyze`: scripts to analyze the results and create some output (plots, JSON with data and fitting parameters, ...)

- `4-export`: bash script to export the group with all results.


## Starting your project

- Create a new virtual environment (if you want) to have different packages than 
- `pip install aiida-core==1.6.*` (until further notice, do NOT use more modern versions)
- install your plugin package (e.g.: `pip install aiida_quantumespresso`)
- `verdi quicksetup` and create a new profile, e.g. `acwf-verification` or any name you want
- (optional) If you want to avoid having to specify the profile each time: `verdi profile setdefault acwf-verification`
- `pip install aiida-common-workflows==1.0.*` (until further notice, do NOT use more modern versions)
- `reentry scan`
- Configure your computer and code (ideally, prepare some .yml config files so it's easy to recreate them if needed)
- Install needed quantities (e.g. `aiida-pseudo install sssp -v 1.1 -x PBE -p precision`)

At this point, you can already submit some tests if you want to check if everything is correclty setup:
  - `aiida-common-workflows launch eos YOUR_CODE_TYPE -p fast -S Si -X YOUR_CODE_LABEL@YOUR_COMPUTER_LABEL -r none -m 1 -w 3600 -d`
    - Note: if you have to add options in the metadata of the calculations, you can use this command line flag: `--engine-options='{"relax": {"account": "mr0"}}'`


## Adapting your scripts and running

Install also the packages in the `requirements.txt`.

Then go through the folders, in order, and follow the instructions in the respective README files.
