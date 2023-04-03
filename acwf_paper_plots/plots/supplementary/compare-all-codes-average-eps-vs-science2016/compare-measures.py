#!/usr/bin/env python
import json
import os
import pylab as pl
import numpy as np
#from adjustText import adjust_text
import acwf_paper_plots.quantities_for_comparison as qc

DO_ZOOM_PANEL = False

EXCELLENT_EPS_THR = 0.06
GOOD_EPS_THR = 0.1


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

DATA_FOLDER = "../../../code-data"
with open(os.path.join(DATA_FOLDER, "labels.json")) as fhandle:
    labels_data = json.load(fhandle)
all_methods = sorted(labels_data['methods-main'].keys())

raw_data = {}
for set_name in ['unaries', 'oxides']:
    raw_data[set_name] = {}
    for method in all_methods:
        with open(os.path.join(DATA_FOLDER, labels_data['methods-main'][method][set_name])) as fhandle:
            raw_data[set_name][method] = json.load(fhandle)

# measure = "epsilon"
data = {}
for i, method1 in enumerate(all_methods):
    for j, method2 in enumerate(all_methods):
        if j < i:
            # Compute only once per pair and avoid i==j that would give zero distance
            continue
        print(f"Computing {method1} vs {method2}...")
        data[(method1, method2)] = {}
        for set_name in ["unaries", "oxides"]:
            # Compute here epsilon between the two sets
            for structure in raw_data[set_name][method1]['BM_fit_data'].keys():
                element, configuration = structure.split('-')
                BM_fit_data_1 = raw_data[set_name][method1]['BM_fit_data'][structure]
                if BM_fit_data_1 is None:
                    continue
                scaling_factor_1 = qc.get_volume_scaling_to_formula_unit(
                    raw_data[set_name][method1]['num_atoms_in_sim_cell'][f'{element}-{configuration}'],
                    element, configuration
                )
                # Here I normalize quantities, so that they are now per atom and not per formula unit!
                # This does not change anything for epsilon and nu, but changes for delta
                V0=BM_fit_data_1['min_volume']/scaling_factor_1
                B0=BM_fit_data_1['bulk_modulus_ev_ang3']
                B01=BM_fit_data_1['bulk_deriv']

                # Get the data for the compare_with plugin, if specified (and if the EOS worked for the 
                # reference plugin, otherwise we don't know which E0 to use)
                try:
                    BM_fit_data_2 = raw_data[set_name][method2]['BM_fit_data'][structure]
                    if BM_fit_data_2 is None:
                        continue               
                except KeyError:
                    continue
                scaling_factor_2 = qc.get_volume_scaling_to_formula_unit(
                    raw_data[set_name][method1]['num_atoms_in_sim_cell'][f'{element}-{configuration}'],
                    element, configuration
                )
                # Here I normalize quantities, so that they are now per atom and not per formula unit!
                # This does not change anything for epsilon and nu, but changes for delta
                CV0=BM_fit_data_2['min_volume']/scaling_factor_2
                CB0=BM_fit_data_2['bulk_modulus_ev_ang3']
                CB01=BM_fit_data_2['bulk_deriv']

                prefactor = 1.
                # Next two not used
                DEFAULT_wb0 = 1.
                DEFAULT_wb1 = 1.
                data[(method1, method2)][structure]= qc.epsilon(V0,B0,B01,CV0,CB0,CB01,prefactor,DEFAULT_wb0,DEFAULT_wb1)

## Compare delta on old set with nu and epsilon on new set
print("# METHOD EPS_AVERAGE EPS_AVERAGE_SCIENCE_SUBSET")
all_eps_average = []
all_eps_subset_average = []
labels = []
for method_pair in data.keys():
    method1, method2 = method_pair
    eps_data = list(data[method_pair].values())
    eps_subset_data = [v for k, v in data[method_pair].items() if k in overlapping_elements]
    print(f"{method_pair} {np.mean(eps_data)} {np.mean(eps_subset_data)}")
    all_eps_average.append(np.mean(eps_data))
    all_eps_subset_average.append(np.mean(eps_subset_data))
    labels.append(f'{method1}-{method2}')

xdata = all_eps_subset_average
xlabel = r"Average $\varepsilon$ (on Science 2016 subset)"
ydata = all_eps_average
ylabel = r"Average $\varepsilon$"
filename = "average-eps-on-science-subset-vs-eps.pdf"
fig = pl.figure(figsize=(5,5))
#print(xdata)
#print(ydata)
pl.plot(xdata, ydata, 'o')
#for label, x, y in zip(labels, xdata, ydata):
#    pl.text(x, y, label)
pl.xlabel(xlabel)
pl.ylabel(ylabel)
#adjust_text(texts)

ZOOM_RANGE = 0.4
pl.plot([0, ZOOM_RANGE], [0, ZOOM_RANGE], '-k')
pl.xlim(0, ZOOM_RANGE)
pl.ylim(0, ZOOM_RANGE)

if DO_ZOOM_PANEL:
    left, bottom, width, height = [0.55, 0.2, 0.25, 0.25]
    ax2 = fig.add_axes([left, bottom, width, height])
    pl.plot(xdata, ydata, 'o')
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel)
    for item in (
        [ax2.title, ax2.xaxis.label, ax2.yaxis.label] +
            ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(8)
    pl.plot([0, 0.8], [0, 0.8], '-k')
    ax2.set_xlim([0,0.15])
    ax2.set_ylim([0,0.15])


pl.plot([EXCELLENT_EPS_THR, EXCELLENT_EPS_THR], [0., 0.3], '--', color='green')
pl.plot([0., EXCELLENT_EPS_THR], [EXCELLENT_EPS_THR, EXCELLENT_EPS_THR], '--', color='green')
pl.plot([GOOD_EPS_THR, GOOD_EPS_THR], [0., 0.3], '--', color='red')
pl.plot([0., GOOD_EPS_THR], [GOOD_EPS_THR, GOOD_EPS_THR], '--', color='red')

pl.savefig(filename)
print(f"File '{filename}' written.")
