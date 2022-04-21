#!/usr/bin/env python
import json
import os
import sys

import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
import tqdm

import quantities_for_comparison as qc

# Adapt this factor to change the zoom on the x axis
# The default zoom is obtained from the standard deviation of the data
# This is a multiplicative number; a number > 1 means zoom in, a number < 1 means zoom out


BINS = 100
DEFAULT_PREFACTOR = 100
DEFAULT_wb0 = 1.0/8.0
DEFAULT_wb1 = 1.0/64.0
EXPECTED_SCRIPT_VERSION = ["0.0.3","0.0.4"]
LIMITS = {"V0_rel_diff":1,"B0_rel_diff":3,"B1_rel_diff":15}


def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


quantity_for_comparison_map = {
    "delta_per_formula_unit": qc.delta,
    "B0_rel_diff": qc.B0_rel_diff,
    "V0_rel_diff": qc.V0_rel_diff,
    "B1_rel_diff": qc.B1_rel_diff,
    "rel_errors_vec_length": qc.rel_errors_vec_length,
    "epsilon": qc.epsilon
}


if __name__ == "__main__":
    try:
        SET_NAME = sys.argv[1]
    except IndexError:
        print(f"The first argument must be the set name.")
        sys.exit(1)

    compare_plugin_data = {}
    
    try:
        with open(f'results-{SET_NAME}-fleur.json') as fhandle:
            compare_plugin_data["fleur"] = json.load(fhandle)
            if not compare_plugin_data["fleur"]['script_version'] in EXPECTED_SCRIPT_VERSION:
                raise ValueError(
                    f"This script only works with data generated at version {EXPECTED_SCRIPT_VERSION}. "
                    f"Please re-run ./get_results.py to update the data format for results-{SET_NAME}-fleur.json!"
                    )
                sys.exit(1)
    except OSError:
        print(f"No data found for fleur (set '{SET_NAME}')")

    try:
        with open(f'results-{SET_NAME}-wien2k.json') as fhandle:
            compare_plugin_data["wien2k"] = json.load(fhandle)
            if not compare_plugin_data["wien2k"]['script_version'] in EXPECTED_SCRIPT_VERSION:
                raise ValueError(
                    f"This script only works with data generated at version {EXPECTED_SCRIPT_VERSION}. "
                    f"Please re-run ./get_results.py to update the data format for results-{SET_NAME}-wien2k.json!"
                    )
                sys.exit(1)
    except OSError:
        print(f"No data found for wien2k (set '{SET_NAME}')")

    if not compare_plugin_data:
        print(f"At least one file with name among [results-{SET_NAME}-fleur.json, results-{SET_NAME}-wien2k.json] should be present.")
        sys.exit(1)

    file_prefix = f'results-{SET_NAME}-'
    file_suffix = '.json'
    results_folder = os.curdir
    code_results = {}
    for fname in os.listdir(results_folder):
        if fname.startswith(file_prefix) and fname.endswith(file_suffix):
            label = fname[len(file_prefix):-len(file_suffix)]
            if label in ["fleur","wien2k"]:
                continue
            with open(os.path.join(results_folder, fname)) as fhandle:
                code_results[label] = json.load(fhandle)
                if not code_results[label]['script_version'] in EXPECTED_SCRIPT_VERSION:
                    raise ValueError(
                        f"This script only works with data generated at version {EXPECTED_SCRIPT_VERSION}. "
                        f"Please re-run ./get_results.py to update the data format for {fname}! Skipping {label}"
                        )
                    code_results.pop(label)
    
    for plugin, plugin_data in code_results.items():
 
        print(f"Analyzing {plugin}")

        name_file = f'histo-{SET_NAME}-{plugin}.pdf'

        all_systems = set(plugin_data['eos_data'].keys())
        all_systems = set(plugin_data['BM_fit_data'].keys())

        # Plotting
        fig, ax = pl.subplots(1, 3, figsize=(18,6))

        SMALL_SIZE = 22
        MEDIUM_SIZE = 24
        BIGGER_SIZE = 28

        pl.rc('font', size=SMALL_SIZE)# controls default text sizes
        pl.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
        pl.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
        pl.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        pl.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        pl.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        pl.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

        
        for indx, QUANTITY in enumerate(["V0_rel_diff","B0_rel_diff","B1_rel_diff"]):

            max_hist = 0
            index=-1
            for name, compare_plugin in compare_plugin_data.items():
                index=index+1
                collect = []

                print(f"  comparing with {name}")
                progress_bar = tqdm.tqdm(sorted(all_systems))
                for element_and_configuration in progress_bar:
                    progress_bar.set_description(f"{element_and_configuration:12s}")
                    progress_bar.refresh()

                    element, configuration = element_and_configuration.split('-')
                    # Get the data for the reference plugin
                    ref_BM_fit_data = plugin_data['BM_fit_data'][f'{element}-{configuration}']
                
                    if ref_BM_fit_data is None:
                        continue
               
                    scaling_factor_ref = qc.get_volume_scaling_to_formula_unit(
                            plugin_data['num_atoms_in_sim_cell'][f'{element}-{configuration}'],
                            element, configuration
                        )

                    V0=ref_BM_fit_data['min_volume']/scaling_factor_ref
                    B0=ref_BM_fit_data['bulk_modulus_ev_ang3']
                    B01=ref_BM_fit_data['bulk_deriv']

                    # Get the data for the compare_with plugin, if specified (and if the EOS worked for the 
                    # reference plugin, otherwise we don't know which E0 to use)
                    try:
                        compare_BM_fit_data = compare_plugin['BM_fit_data'][f'{element}-{configuration}']
                        if compare_BM_fit_data is None:
                            # No fitting data in the plugin to compare with.
                            # Raise this exception that is catched one line below, so
                            # it will set `compare_eos_fit_energy` to None.
                            raise KeyError                    
                    except KeyError:
                        # Set to None if fit data is missing (if we are here, the EOS points
                        # are there, so it means that the fit failed). I will still plot the
                        # points
                        continue

                    scaling_factor_comp = qc.get_volume_scaling_to_formula_unit(
                            compare_plugin['num_atoms_in_sim_cell'][f'{element}-{configuration}'],
                            element, configuration
                        )

                    CV0=compare_BM_fit_data['min_volume']/scaling_factor_comp
                    CB0=compare_BM_fit_data['bulk_modulus_ev_ang3']
                    CB01=compare_BM_fit_data['bulk_deriv']

                    quant = quantity_for_comparison_map[QUANTITY](V0,B0,B01,CV0,CB0,CB01,DEFAULT_PREFACTOR,DEFAULT_wb0,DEFAULT_wb1)

                    collect.append(quant)

                sta_dev = np.std(np.array(collect))
                lim = LIMITS[QUANTITY]
                mean = np.mean(collect)
                #std = np.std(collect)


                hist_y, bins, patches = ax[indx].hist(
                    collect,
                    bins=BINS,
                    range=[-lim,lim],
                    alpha=0.5,
                    label=f"{name} {round(mean,3)} {round(sta_dev,2)}"
                )
                countBig = 0
                countSmall = 0
                for alls in collect:
                    if alls > lim:
                        countBig = countBig + 1
                    if alls < -lim:
                        countSmall = countSmall + 1
                if countBig > 0:
                    ax[indx].annotate(f"+{countBig} {name}", xy=(lim, max(hist_y)/(2+index)), xytext=(0.1*lim, max(hist_y)/(2+index)), arrowprops=dict(facecolor='black', shrink=0.05))
                if countSmall:
                    ax[indx].annotate(f"+{countSmall} {name}", xy=(-lim, max(hist_y)/(2+index)), xytext=(-lim+lim/5, max(hist_y)/(2+index)), arrowprops=dict(facecolor='black', shrink=0.05))

                if index == 0:
                    ax[indx].axvline(mean, color='b', linestyle=':')#, label=f"mean {round(mean,3)}, std {round(sta_dev,3)}")
                else:
                    ax[indx].axvline(mean, color='r', linestyle=':')
           
                max_hist = max(max_hist,max(hist_y))
            # Reset the xlim
            ax[indx].set_xlim([-lim, lim])
            ax[indx].set_ylim([0,max_hist+max_hist/5])
            #ax[indx].annotate(f"mean {round(mean,3)}",xy=(-lim+lim/20,max(hist_y)-max(hist_y)/10)) 
            #ax[indx].annotate(f"std {round(sta_dev,2)}",xy=(-lim+lim/20,max(hist_y)-max(hist_y)/5))

            ax[indx].legend(loc='upper center')
            ax[indx].set_xlabel(f"% difference in {QUANTITY.strip('_rel_diff')}",fontsize=22)
            ax[indx].set_ylabel("Frequency",fontsize=22)
            ax[indx].tick_params(axis="x",labelsize=20)
            ax[indx].tick_params(axis="y",labelsize=20)
            #set(xlabel=f"{DEFAULT_PREFACTOR}*{QUANTITY}", ylabel='Frequency', Fontsize=30)
        if "unaries" in SET_NAME:
            fig.suptitle(f"{plugin} unaries")
        else:
            fig.suptitle(f"{plugin} oxides")
        fig.tight_layout()
        fig.savefig(f"{name_file}")
        pl.close(fig)
