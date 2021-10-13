# commonwf-oxides-scripts
Scripts to facilitate the common-workflow project on oxides

## Content of the folders

These scripts are a template (based on Quantum ESPRESSO) that can be adapted relatively simply for any other code.

- `0-preliminary-do-not-run`: ignore this folder - these are helper scripts to create the import the structures from CIF into AiiDA nodes and create initial AiiDA archive file. The scripts will be run only once by one of us.

- `1-preliminary`: scripts to import the group with the initial structures, and to create a subgroup containing only the nodes you actually want to run.

- `2-submit`: scripts to submit (in batches) the common EOS workflows for all systems

- `3-analyze`: scripts to analyze the results and create some output (plots, JSON with data and fitting parameters, ...)

- `4-export`: bash script to export the group with all results.


## Starting your project

- Create a new virtual environment (if you want) to have different packages than 
- `pip install aiida-core==1.6.3` (until further notice, do NOT use develop)
- install your plugin package (e.g.: `pip install aiida_quantumespresso`)
- `verdi quicksetup` and create a new profile, e.g. `commonwf` or any name you want
- (optional) If you want to avoid having to specify the profile each time: `verdi profile setdefault commonwf`
- `pip install aiida-common-workflows` (you can use the most common develop branch)
- `reentry scan`
- Configure your computer and code (ideally, prepare some .yml config files so it's easy to recreate them if needed)
- Install needed quantities (e.g. `aiida-pseudo install sssp -v 1.1 -x PBE -p precision`) and document them (if we need to report them in the paper, similarly to the supplementary of the first paper)

At this point, you can already submit some tests if you want to check if everything is correclty setup:
  - `aiida-common-workflows launch eos YOUR_CODE_TYPE -p precise -S Si -X YOUR_CODE_LABEL@YOUR_COMPUTER_LABEL -r none -m 1 -w 3600 -d`
    - Note: if you have to add options in the metadata of the calculations, you can use this command line flag that we recently introduced: `--engine-options='{"relax": {"account": "mr0"}}'`


## Adapting your scripts and running

Install also the packages in the `requirements.txt`.

Then go through the folders, in order, and follow the instructions in the respective README files.
