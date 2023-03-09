#!/usr/bin/env python
import os
import glob
import json
import copy
from acwf_paper_plots.quantities_for_comparison import get_num_atoms_in_formula_unit

PREFIX = "nu-unaries-"
SUFFIX = "-vs-ae.json"
ALL_MEASURES = ["nu", "epsilon", "delta_per_formula_unit"]

all_methods = []
# assume that if the file for nu and unaries is there, all other files are there as well
for filename in glob.glob(f"{PREFIX}*{SUFFIX}"):
    method_short_name = filename[len(PREFIX):-len(SUFFIX)]
    all_methods.append(method_short_name)

data = {}
for measure in ALL_MEASURES:
    data[measure] = {}
    for method in all_methods:
        data[measure][method] = {}
        for set_name in ["unaries", "oxides"]:
            with open(f"{measure}-{set_name}-{method}{SUFFIX}") as fhandle:
                data_from_file = json.load(fhandle)
                if measure == "delta_per_formula_unit":
                    # Convert in delta per atom
                    for k in data_from_file:
                        conf = k.split('-')[1]
                        data_from_file[k] /= get_num_atoms_in_formula_unit(conf)

                # print('>', measure, method, set_name, len(data_from_file))
                data[measure][method].update(data_from_file)

# Check that, given a method, the keys are the same for all measures.
# This should always be the case! If I have data, I should have all measures
flat_data = {}
for method in all_methods:
    ref_keys = None    
    for measure in data:
        if ref_keys is None:
            ref_keys = set(data[measure][method].keys())
            #print(method, len(ref_keys))
        else:
            assert set(data[measure][method].keys()) == ref_keys, f"Different set of keys for method {method}"

    # Generate flat lists, ensuring the order is the same
    ref_keys = sorted(ref_keys)
    data_for_measures = {}
    for measure in data:
        data_for_measures[measure] = []
        for key in ref_keys:
            data_for_measures[measure].append(data[measure][method][key])
    flat_data[method] = copy.deepcopy(data_for_measures)

import pylab as pl

for (meas1, meas1name), (meas2, meas2name), file_basename in [
    [
        ("epsilon", r"$\varepsilon$"),
        ("nu", r"$\nu$"),
        "epsilon-vs-nu"
    ],
    [
        ("delta_per_formula_unit", r"$\Delta$ per atom"),
        ("epsilon", r"$\varepsilon$"),
        # Note: above I have converted it in delta per atom
        "delta-vs-epsilon"
    ],
    [
        # Note: above I have converted it in delta per atom
        ("delta_per_formula_unit", r"$\Delta$ per atom"),
        ("nu", r"$\nu$"),
        "delta-vs-nu"
    ],        
]:
    fig = pl.figure()
    for method in all_methods:
        pl.plot(flat_data[method][meas1], flat_data[method][meas2], '.', color='#2b8cbe', label=method)
    pl.xlabel(meas1name)
    pl.ylabel(meas2name)
    #pl.legend(loc='best')

    filename = f"comparison-{file_basename}.png"
    pl.savefig(filename)
    print(f"File '{filename}' written.")

    
