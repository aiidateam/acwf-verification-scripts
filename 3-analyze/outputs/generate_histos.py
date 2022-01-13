#!/usr/bin/env python
import json
import os
import sys

import numpy as np
import pylab as pl
import tqdm

import quantities_for_comparison as qc

def get_plugin_name():
    file_name = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        os.pardir, os.pardir, 'plugin_name.txt'
    )
    try:
        with open(file_name) as fhandle:
            plugin_name = fhandle.read().strip()
            # Simple check e.g. to make sure there are no weird characters,
            # newlines, ... - one might still make a typo, but at least we
            # do a basic check
            assert plugin_name.isidentifier()
        return plugin_name
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "You need to define a file `../../plugin_name.txt`, containing the "
            "name of your plugin (siesta, quantum_espresso, ...) in the format "
            "expected by the aiida-common-workflows project"
        ) from exc

PLUGIN_NAME = get_plugin_name()

BINS = 100
DEFAULT_PREFACTOR = 100
DEFAULT_wb0 = 0.1
DEFAULT_wb1 = 0.01


quantity_for_comparison_map = {
    "delta": qc.delta,
    "delta_per_atom": qc.delta_per_atom,
    "B0_rel_diff": qc.B0_rel_diff,
    "V0_rel_diff": qc.V0_rel_diff,
    "B1_rel_diff": qc.B1_rel_diff,
    "rel_errors_vec_length": qc.rel_errors_vec_length,
    "epsilon2": qc.epsilon2
}


if __name__ == "__main__":
    try:
        QUANTITY = sys.argv[1]
    except IndexError:
        print(f"The first argument must be the quantity to use for comparison. Choose among {quantity_for_comparison_map.keys()}")
        sys.exit(1)

    if QUANTITY not in quantity_for_comparison_map.keys():
        print(f"The first argument must be the quantity to use for comparison. Choose among {quantity_for_comparison_map.keys()}")
        sys.exit(1)

    all_args = sys.argv[2:]

    if not all_args:
        print("The plugin's names whose results will be plotted must be listed explicitely as script arguments.")
        sys.exit(1)
    
    try:
        with open(f'results-{PLUGIN_NAME}.json') as fhandle:
            reference_plugin_data = json.load(fhandle)
    except OSError:
        print(f"No data found for your plugin '{PLUGIN_NAME}'. Did you run `./get_results.py` first?")
        sys.exit(1)
    
    print(f"Using data for plugin '{PLUGIN_NAME}' compared with {all_args}.")

    compare_plugin_data = []
    for compare_with in all_args:
        try:
            with open(f'results-{compare_with}.json') as fhandle:
                compare_plugin_data.append(json.load(fhandle))
        except OSError:
            print(f"No data found for the plugin '{compare_with}': you need the file results-{compare_with}.json.")
            sys.exit(1)

    name_file = f'histo-{QUANTITY}-{PLUGIN_NAME}'

    all_systems = set(reference_plugin_data['eos_data'].keys())
    all_systems = set(reference_plugin_data['BM_fit_data'].keys())
    #all_systems.update(compare_plugin_data['BM_fit_data'].keys())

    # Plotting
    fig = pl.figure(figsize=(18,6))

    SMALL_SIZE = 20
    MEDIUM_SIZE = 24
    BIGGER_SIZE = 28

    pl.rc('font', size=SMALL_SIZE)# controls default text sizes
    pl.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
    pl.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    pl.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    pl.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    pl.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    pl.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    
    for index, compare_plugin in enumerate(compare_plugin_data):
        
        collect = []

        print(f"comparing with {all_args[index]}")
        progress_bar = tqdm.tqdm(sorted(all_systems))
        for element_and_configuration in progress_bar:
            progress_bar.set_description(f"{element_and_configuration:12s}")
            progress_bar.refresh()

            element, configuration = element_and_configuration.split('-')
            # Get the data for the reference plugin
            ref_BM_fit_data = reference_plugin_data['BM_fit_data'][f'{element}-{configuration}']
        
            if ref_BM_fit_data is None:
                continue
        
            V0=ref_BM_fit_data['min_volume']
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

            CV0=compare_BM_fit_data['min_volume']
            CB0=compare_BM_fit_data['bulk_modulus_ev_ang3']
            CB01=compare_BM_fit_data['bulk_deriv']

            quant = quantity_for_comparison_map[QUANTITY](V0,B0,B01,CV0,CB0,CB01,configuration,DEFAULT_PREFACTOR,DEFAULT_wb0,DEFAULT_wb1)

            collect.append(quant)

        mini = min(collect)

        if mini >= 0:
            sta_dev=np.sqrt(np.mean(np.array(collect)**2))
            pl.hist(collect, bins=BINS, range=[0, sta_dev], label=f"{all_args[index]}", alpha=0.5)
            countBig = 0
            for alls in collect:
                if alls > sta_dev:
                    countBig = countBig + 1
            if countBig > 0:
                pl.annotate(f"{countBig} more for {all_args[index]}", xy=(pl.xlim()[1], (pl.ylim()[1]-pl.ylim()[0])/2/(index+1)), xytext=(pl.xlim()[1]-0.5*sta_dev, (pl.ylim()[1]-pl.ylim()[0])/2/(index+1)), arrowprops=dict(facecolor='black', shrink=0.05))

        else:
            sta_dev = np.std(np.array(collect))
            pl.hist(collect, bins=BINS, range=[-2*sta_dev, 2*sta_dev], label=f"{all_args[index]}", alpha=0.5)
            countBig = 0
            countSmall = 0
            for alls in collect:
                if alls > 2*sta_dev:
                    countBig = countBig + 1
                if alls < -2*sta_dev:
                    countSmall = countSmall + 1
            if countBig > 0:
                pl.annotate(f"{countBig} more for {all_args[index]}", xy=(pl.xlim()[1], (pl.ylim()[1]-pl.ylim()[0])/2/(index+1)), xytext=(pl.xlim()[1]-1.5*sta_dev, (pl.ylim()[1]-pl.ylim()[0])/2/(index+1)), arrowprops=dict(facecolor='black', shrink=0.05))
            if countSmall:
                pl.annotate(f"{countSmall} more for {all_args[index]}", xy=(pl.xlim()[0], (pl.ylim()[1]-pl.ylim()[0])/2/(index+1)), xytext=(pl.xlim()[0]+0.2*sta_dev, (pl.ylim()[1]-pl.ylim()[0])/2/(index+1)), arrowprops=dict(facecolor='black', shrink=0.05))
   
    
    pl.legend(loc='upper right')
    if QUANTITY in ["delta_per_atom","delta"]:
        pl.xlabel(f"{QUANTITY}")
    elif QUANTITY == "rel_errors_vec_length":
        pl.xlabel(f"{DEFAULT_PREFACTOR}*{QUANTITY}({DEFAULT_wb0},{DEFAULT_wb1})")
    else:
        pl.xlabel(f"{DEFAULT_PREFACTOR}*{QUANTITY}")
    pl.ylabel("Frequency")
    pl.title(f"{PLUGIN_NAME}")
    pl.tight_layout()
    pl.savefig(f"{name_file}")
    pl.close(fig)
