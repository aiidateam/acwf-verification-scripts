# Analyzing the results

In this folder you find scripts to help you analyze your results.

(Note: the folders `eos_utils` and `cottenier_data_set1` are needed by the scripts, you don't need to modify them).

# get_errors.py

After modifying your plugin name, run it to store some information on the EOS workchains that failed.

# get_results.py

The main script to get results from your calculations.
Also in this case, the only thing you should do is to adapt the `PLUGIN_NAME` variable at the beginning.

Then, run it and wait. It will fit all EOS, and create text files and PNG plots in the `outputs` folder, that you can inspect.

# Creating a single collated PDF
If you want to create a single PDF from all PNG plots:

- go in the folder `collate-plots`
- adapt the `create_latex_file.py` to use your plugin name
- run the `create_latex_file.py`
- enter the `tex-template` subfolder
- run `pdflatex results.tex`
- inspect the output `results.pdf`

