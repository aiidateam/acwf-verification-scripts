# Adapting the script to launch a single calculation

You now need to prepare the script of the common workflow to run the calculations with your code.

In principle, being the common interface generic, the amount of changes should be small: but there might still be code-specific changes to do (e.g. for optimal parallelization).

You can change the example file `first-tests/launch_example_one_calc_only_template.py`. 

This will keep running the same system (Ag2O), so you can check if the script is correct.

You need to change at least the initial `PLUGIN_NAME` and the `CODE_LABEL` at the beginning of the script.

Also adapt the part generating the `engines` dict, and finally either remove or adapt the `sub_process` input in the `inputs` dict.

When you are done, just execute the script with `verdi run` and adapt it until when you're happy with the result.

# Preparing the script to run all systems

You are now ready to create your own launch script.
We suggest that you copy the file `launch_calculations_qe.py`, rename it as `launch_calculations_<PLUGIN_NAME>.py` and modify as follows:

- `DRY_RUN = True`: Change it to False when you are ready to run. At the beginning you can keep it to True to check if the script is correct without running.
- `MAX_CONCURRENT = 200`: change to the maximum number of workchains that are allowed to run at any time (note! this is the number of *workchains* - the actual number of CalcJobs submitted to the queue will be 6 times larger at the end of the EOS as points are submitted in parallel)
- `PLUGIN_NAME = 'quantum_espresso'`: replace with your plugin name
- `CODE_LABEL = 'qe-6.7-pw@daint-mc'`: replace with your code name
-  adapt the content of the `get_inputs_and_processclass_from_extras` method of the EosSubmissionController class reusing/copying the code that you wrote in the previous step in the file `launch_example_one_calc_only.py`.

NOTE: the script uses classes of the `aiida-submission-controller` package that must be installed before running the script (see ../requirements.txt)

When you are ready, you can execute it to see if everything works as expected.
The output will contain information on how many simulations are still to run, how many
can actually be run before reaching the `MAX_CONCURRENT` limit, ...

When you are ready, just set `DRY_RUN` to False and then run the script in a loop, e.g. (in bash, maybe in a `screen` session):
```bash
while true ; do verdi run launch_calculations.py ; sleep 60 ; done
```
(you can adapt the sleep, probably you might even increase it to 5 minutes or so).

Continue with the next folder when all simulations are done.


## Note: failed simulations, or simulations to re-run
The `aiida-submission-controller` will put submitted workchains in a group
`commonwf-oxides/set2/workflows/PLUGIN_NAME` (the set name, here `set2`, could be different - use the most recent set available).   

If you want to rerun one or more workflows, just remove it from the group (I suggest to create another group `commonwf-oxides/set2/workflows/PLUGIN_NAME/failed` [remember to use the correct set name], and actually also add the nodes there, so you don't lose track of them, they might be useful for futher analysis). Or you can just delete the workflows.
