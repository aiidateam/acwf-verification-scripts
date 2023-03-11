#!/usr/bin/env python
import glob
import json
import copy

import numpy as np

from acwf_paper_plots.quantities_for_comparison import get_num_atoms_in_formula_unit

PREFIX = "nu-unaries-"
SUFFIX = "-vs-ae.json"
ALL_MEASURES = ["nu", "epsilon", "delta_per_formula_unit"]

# List of unaries that are from the Science 2016 paper that also in our set
overlapping_elements = [
    "Ag-X/FCC",
    "Al-X/FCC",
    "Ar-X/FCC",
    "Au-X/FCC",
    "Ca-X/FCC",
    "Cu-X/FCC",
    "Ir-X/FCC",
    "Kr-X/FCC",
    "Ne-X/FCC",
    "Pb-X/FCC",
    "Pd-X/FCC",
    "Pt-X/FCC",
    "Rh-X/FCC",
    "Rn-X/FCC",
    "Sr-X/FCC",
    "Xe-X/FCC",
    "Ba-X/BCC",
    "Cs-X/BCC",
    "K-X/BCC",
    "Mo-X/BCC",
    "Nb-X/BCC",
    "Rb-X/BCC",
    "Ta-X/BCC",
    "V-X/BCC",
    "W-X/BCC",
    "Po-X/SC",
    "Ge-X/Diamond",
    "Si-X/Diamond",
    "Sn-X/Diamond"
]

all_methods = []
# assume that if the file for nu and unaries is there, all other files are there as well
for filename in glob.glob(f"{PREFIX}*{SUFFIX}"):
    method_short_name = filename[len(PREFIX):-len(SUFFIX)]
    all_methods.append(method_short_name)
all_methods = sorted(all_methods)

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

## Generate pairwise comparison of metrics
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

## Compare delta on old set with nu and epsilon on new set
print("# METHOD EPS_AVERAGE NU_AVERAGE DELTA_AVERAGE_SCIENCE_SUBSET")
all_eps_average = []
all_nu_average = []
all_delta_average = []
all_delta_subset_average = []
for method in all_methods:
    eps_data = list(data['epsilon'][method].values())
    nu_data = list(data['nu'][method].values())
    delta_data = list(data['delta_per_formula_unit'][method].values())
    delta_subset_data = [v for k, v in data['delta_per_formula_unit'][method].items() if k in overlapping_elements]
    #print(f"# Method: {method} ({len(eps_data)}/960 systems, {len(delta_subset_data)}/{len(overlapping_elements)} in the Delta subset)")
    #print(f"#  - average epsilon          : {np.mean(eps_data)}")
    #print(f"#  - average nu               : {np.mean(nu_data)}")
    #print(f"#  - average delta (on subset): {np.mean(delta_subset_data)}")
    print(f"{method} {np.mean(eps_data)} {np.mean(nu_data)} {np.mean(delta_subset_data)}")
    all_eps_average.append(np.mean(eps_data))
    all_nu_average.append(np.mean(nu_data))
    all_delta_average.append(np.mean(delta_data))
    all_delta_subset_average.append(np.mean(delta_subset_data))

for xdata, xlabel, ydata, ylabel, filename in [
    [all_delta_subset_average, r"Average $\Delta$ per atom (on Science 2016 subset)", all_eps_average, r"Average $\varepsilon$", "average-delta-vs-eps-on-science-subset.png"],
    [all_delta_subset_average, r"Average $\Delta$ per atom (on Science 2016 subset)", all_nu_average, r"Average $\nu$", "average-delta-vs-nu-on-science-subset.png"],
    [all_eps_average, r"Average $\varepsilon$", all_nu_average, r"Average $\nu$", "average-eps-vs-nu-on-science-subset.png"],
    [all_delta_subset_average, r"Average $\Delta$ per atom (on Science 2016 subset)", all_delta_average, r"Average $\Delta$ per atom (full set)", "average-delta-subset-vs-full-delta-on-science-subset.png"],
]:
    pl.figure()
    print(xdata)
    print(ydata)
    pl.plot(xdata, ydata, 'o')
    for label, x, y in zip(all_methods, xdata, ydata):
        pl.annotate(label, (x, y))
    pl.xlabel(xlabel)
    pl.ylabel(ylabel)
    pl.savefig(filename)
    print(f"File '{filename}' written.")
