# Analyzing the results

In this folder you find scripts to help you analyze your results.

(Note: the folders `eos_utils` and `cottenier_data_set1` are needed by the scripts, you don't need to modify them).

# get_errors.py

Run it to store some information on the EOS workchains that failed.
It expects that you already created the file `../plugin_name.txt` (see README file in the folder `1-preliminary` for more details).

# get_results.py

The main script to get results from your calculations.
It expects that you already created the file `../plugin_name.txt` (see README file in the folder `1-preliminary` for more details).

After you run your simulations, just run this script it with `verdi run` and wait.
It will fit all EOS, and create text files and PNG plots in the `outputs` folder, that you can inspect.

# Creating a single collated PDF
If you want to create a single PDF from all PNG plots, to inspect graphically the results do the following
(note: currently, it will compare with preliminary Wien2K results - in the future it will
directly perform a cross-code validation):

- go in the folder `collate-plots`
- run the `create_latex_file.py`
- enter the `tex-template` subfolder
- run `pdflatex results.tex`
- inspect the output `results.pdf`

