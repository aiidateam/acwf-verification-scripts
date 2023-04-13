#!/usr/bin/env python
import glob
import json
import copy

import numpy as np
from adjustText import adjust_text

from acwf_paper_plots.quantities_for_comparison import get_num_atoms_in_formula_unit

PREFIX = "nu-unaries-"
SUFFIX = "-vs-ae.json"
ALL_MEASURES = ["nu", "epsilon", "delta_per_formula_unit", "delta_per_formula_unit_over_b0"]

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
                if measure.startswith("delta_"):
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
for LOGLOG in [True, False]:
    # factor1 and 2 are used to change units
    for (meas1, meas1name, factor1), (meas2, meas2name, factor2), file_basename in [
        [
            ("epsilon", r"$\varepsilon$", 1.),
            ("nu", r"$\nu$", 1.),
            "epsilon-vs-nu"
        ],
        [
            ("delta_per_formula_unit", r"$\Delta$ per atom (meV/atom)", 1.),
            ("epsilon", r"$\varepsilon$", 1.),
            # Note: above I have converted it in delta per atom
            "delta-vs-epsilon"
        ],
        [
            # Delta is in meV. B_0 is in eV/ang^3.
            # If I want to get Ang^-3, I need to divide by 1000
            ("delta_per_formula_unit_over_b0", r"$\Delta/B_0$ per atom ($\AA^{3}$)", 0.001),
            ("epsilon", r"$\varepsilon$", 1.),
            # Note: above I have converted it in delta per atom
            "deltaoverB0-vs-epsilon"
        ],    
        [
            # Note: above I have converted it in delta per atom
            ("delta_per_formula_unit", r"$\Delta$ per atom (meV/atom)", 1.),
            ("nu", r"$\nu$", 1.),
            "delta-vs-nu"
        ],        
    ]:
        
        if not LOGLOG and meas1 == 'epsilon' and meas2 == 'nu':
            fig = pl.figure(figsize=(8, 4))
            pl.subplot(121)
        else:
            fig = pl.figure(figsize=(4, 3))
            fig = pl.figure()
        all_data_x = []
        all_data_y = []
        for method in all_methods:
            all_data_x += (np.array(flat_data[method][meas1]) * factor1).tolist()
            all_data_y += (np.array(flat_data[method][meas2]) * factor2).tolist()
            if LOGLOG:
                pl.loglog(np.array(flat_data[method][meas1]) * factor1, np.array(flat_data[method][meas2]) * factor2, '.', color='#2b8cbe', label=method)
            else:
                pl.plot(np.array(flat_data[method][meas1]) * factor1, np.array(flat_data[method][meas2]) * factor2, '.', color='#2b8cbe', label=method)
        if LOGLOG:
            # Cache limits
            xlim = pl.xlim()
            ylim = pl.ylim()

            if meas1 == 'epsilon' and meas2 == 'nu':
                xmin = 1e-3; xmax = 1; slope = 1.673

                for m1_threshold, m2_threshold, color in [
                    (0.06, 0.1, 'green'),
                    (0.2, 0.35, 'red'),
                ]:
                    pl.plot([xlim[0], m2_threshold / slope], [m2_threshold, m2_threshold], '--', color=color, linewidth=0.5)
                    pl.plot([m1_threshold, m1_threshold], [ylim[0], m1_threshold * slope], '--', color=color, linewidth=0.5)
                    pl.annotate('(d)', (0.03, 0.95), xycoords='axes fraction')

            elif meas1 == 'delta_per_formula_unit' and meas2 == 'epsilon':
                xmin = 1e-3; xmax = 10; slope = 0.2

                for m1_threshold, m2_threshold, color in [
                    (0.3, 0.06, 'green'),
                    (0.95, 0.2, 'red'),
                ]:
                    pl.plot([xlim[0], m2_threshold / slope], [m2_threshold, m2_threshold], '--', color=color, linewidth=0.5)
                    pl.plot([m1_threshold, m1_threshold], [ylim[0], m1_threshold * slope], '--', color=color, linewidth=0.5)
                    pl.annotate('(a)', (0.03, 0.95), xycoords='axes fraction')

            elif meas1 == 'delta_per_formula_unit' and meas2 == 'nu':
                xmin = 2e-3; xmax = 30; slope = 0.33

                for m1_threshold, m2_threshold, color in [
                    (0.3, 0.1, 'green'),
                    (0.95, 0.35, 'red'),
                ]:
                    pl.plot([xlim[0], m2_threshold / slope], [m2_threshold, m2_threshold], '--', color=color, linewidth=0.5)
                    pl.plot([m1_threshold, m1_threshold], [ylim[0], m1_threshold * slope], '--', color=color, linewidth=0.5)
                    pl.annotate('(b)', (0.03, 0.95), xycoords='axes fraction')

            elif meas1 == 'delta_per_formula_unit_over_b0' and meas2 == 'epsilon':
                xmin = 1e-5; xmax = 1e-2; slope = 100.
                pl.annotate('(c)', (0.03, 0.95), xycoords='axes fraction')
            else:
                xmin = None
            if xmin is not None:
                # Plot line of slope = 1 (i.e. linear relationship)
                pl.loglog([xmin, xmax], [xmin * slope, xmax * slope], '-k')
            # Put back cached limits
            pl.xlim(xlim)
            pl.ylim(ylim)

        pl.xlabel(meas1name)
        pl.ylabel(meas2name)
        #pl.legend(loc='best')
        if not LOGLOG and meas1 == 'epsilon' and meas2 == 'nu':
            all_data_x = np.array(all_data_x)
            all_data_y = np.array(all_data_y)

            # Theoretical fit
            m_theor =6. / np.sqrt(15) #  ~1.549

            threshold = 1/m_theor # So that we go up to nu ~= 1
            filter = np.copy(all_data_x < threshold)
            all_data_x = all_data_x[filter]
            all_data_y = all_data_y[filter]
            pl.subplot(122)
            pl.plot(all_data_x, all_data_y, '.', color='#2b8cbe')
            pl.xlabel(meas1name)
            pl.ylabel(meas2name)


            # Linear fit, see https://numpy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html
            A = np.vstack([all_data_x, np.ones(len(all_data_x))]).T
            m, c = np.linalg.lstsq(A, all_data_y, rcond=None)[0]
            print(f">> eps-vs-nu fit (on linear scale): nu = {m} * eps + {c}")
            xmin, xmax = all_data_x.min(), all_data_x.max()
            pl.plot([xmin, xmax], [m*xmin + c, m*xmax + c], 'r', label='Linear fit') #label=f'Fit (m={m:.3f}, c={c:.3f})')
            pl.plot([xmin, xmax], [m_theor*xmin, m_theor*xmax], 'y', label='Theoretical relation')#label=f'Theoretical (m={m_theor:.3f}, c=0)')
            pl.xlim(0, threshold)
            pl.ylim(0, threshold * m_theor)
            pl.legend(loc='lower right')


        if not LOGLOG and meas1 == "delta_per_formula_unit_over_b0":
            xlim_max = 750

            outliers_count = 0
            for method in all_methods:
                outliers_count += sum(np.array(flat_data[method][meas1]) > xlim_max)

            pl.xlim(0, xlim_max)
            pl.title(f"{outliers_count} outliers on the right")

        filename = f"comparison-{file_basename}{'-loglog' if LOGLOG else ''}.png"
        pl.savefig(filename)
        print(f"File '{filename}' written.")
