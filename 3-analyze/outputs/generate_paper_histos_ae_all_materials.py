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

PLUGIN_NAME = 'wien2k-dk_0.06'

BINS = 100
DEFAULT_PREFACTOR = 100
DEFAULT_wb0 = 1.0/8.0
DEFAULT_wb1 = 1.0/64.0
EXPECTED_SCRIPT_VERSION = ["0.0.3","0.0.4"]
LIMITS = {"V0_rel_diff":0.3,"B0_rel_diff":2,"B1_rel_diff":10}

QUANTITY_FANCY_NAMES = {
    'B0': "$B_0$",
    'V0': "$V_0$",
    'B1': "$B'$"
}


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
    SET_NAME_1 = 'unaries-verification-PBE-v1'
    SET_NAME_2 = 'oxides-verification-PBE-v1'

    reference_plugin_data = []

    try:
        with open(f'results-{SET_NAME_1}-{PLUGIN_NAME}.json') as fhandle:
            reference_plugin_data.append(json.load(fhandle))
    except OSError:
        print(f"No data found for your plugin '{PLUGIN_NAME}' (set '{SET_NAME_1}'). Did you run `./get_results.py` first?")
        sys.exit(1)
 
    if not reference_plugin_data[-1]['script_version'] in EXPECTED_SCRIPT_VERSION:
        raise ValueError(
            f"This script only works with data generated at version {EXPECTED_SCRIPT_VERSION}. "
            f"Please re-run ./get_results.py to update the data format for {PLUGIN_NAME}!"
            )

    try:
        with open(f'results-{SET_NAME_2}-{PLUGIN_NAME}.json') as fhandle:
            reference_plugin_data.append(json.load(fhandle))
    except OSError:
        print(f"No data found for your plugin '{PLUGIN_NAME}' (set '{SET_NAME_2}'). Did you run `./get_results.py` first?")
        sys.exit(1)

    if not reference_plugin_data[-1]['script_version'] in EXPECTED_SCRIPT_VERSION:
        raise ValueError(
            f"This script only works with data generated at version {EXPECTED_SCRIPT_VERSION}. "
            f"Please re-run ./get_results.py to update the data format for {PLUGIN_NAME}!"
            )

    compare_plugin_data = []
    compare_with = 'fleur'
    for SET_NAME in [SET_NAME_1,SET_NAME_2]:
        try:
            with open(f'results-{SET_NAME}-{compare_with}.json') as fhandle:
                compare_plugin_data.append(json.load(fhandle))
            if not compare_plugin_data[-1]['script_version'] in EXPECTED_SCRIPT_VERSION:
                raise ValueError(
                    f"This script only works with data generated at version {EXPECTED_SCRIPT_VERSION}. "
                    f"Please re-run ./get_results.py to update the data format for {compare_with}!"
                    )
        except OSError:
            print(f"No data found for the plugin '{compare_with}' (set '{SET_NAME}'): you need the file results-{SET_NAME}-{compare_with}.json.")
            sys.exit(1)

    name_file = f'histo-final.pdf'

        # Plotting
    #fig = pl.figure(figsize=(18,6))
    fig, ax = pl.subplots(1, 3, figsize=(18,6))

    TINY_SIZE = 18
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
       

        collect = []

        for ind in range(2):
        
            all_systems = set(reference_plugin_data[ind]['eos_data'].keys())
            all_systems = set(reference_plugin_data[ind]['BM_fit_data'].keys())
            #all_systems.update(compare_plugin_data['BM_fit_data'].keys())

       
            progress_bar = tqdm.tqdm(sorted(all_systems))
            for element_and_configuration in progress_bar:
                progress_bar.set_description(f"{element_and_configuration:12s}")
                progress_bar.refresh()

                element, configuration = element_and_configuration.split('-')
                # Get the data for the reference plugin
                ref_BM_fit_data = reference_plugin_data[ind]['BM_fit_data'][f'{element}-{configuration}']
            
                if ref_BM_fit_data is None:
                    continue
            
                scaling_factor_ref = qc.get_volume_scaling_to_formula_unit(
                        reference_plugin_data[ind]['num_atoms_in_sim_cell'][f'{element}-{configuration}'],
                        element, configuration
                    )

                V0=ref_BM_fit_data['min_volume']/scaling_factor_ref
                B0=ref_BM_fit_data['bulk_modulus_ev_ang3']
                B01=ref_BM_fit_data['bulk_deriv']

                # Get the data for the compare_with plugin, if specified (and if the EOS worked for the 
                # reference plugin, otherwise we don't know which E0 to use)
                try:
                    compare_BM_fit_data = compare_plugin_data[ind]['BM_fit_data'][f'{element}-{configuration}']
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
                        compare_plugin_data[ind]['num_atoms_in_sim_cell'][f'{element}-{configuration}'],
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


        hist_y, bins, patches = ax[indx].hist(collect, bins=BINS, range=[-lim,lim], alpha=0.5)
        countBig = 0
        countSmall = 0
        for alls in collect:
            if alls > lim:
                countBig = countBig + 1
            if alls < -lim:
                countSmall = countSmall + 1
        if countBig > 0:
            ax[indx].annotate(f"+{countBig}", xy=(lim, max(hist_y)/2), xytext=(0.5*lim, max(hist_y)/2), arrowprops=dict(color='#999999', shrink=0.05),va='center')
        if countSmall:
            ax[indx].annotate(f"+{countSmall}", xy=(-lim, max(hist_y)/2), xytext=(-lim*0.8, max(hist_y)/2), arrowprops=dict(color='#999999', shrink=0.05),va='center')

        # Fit Gaussian and plot it
        #hist_x = (bins[1:] + bins[:-1])/2
        #popt, pcov = curve_fit(gaussian, hist_x, hist_y, p0=[10., 0., 1.])
        #x = np.linspace(-lim, lim, 1000)
        #sigma = abs(popt[2])
        ## NOTES ON THE RELATION BETWEEN THE SIGMA OF THE GAUSSIAN AND THE FWHM
        #  np.exp(-HWHM**2/(2*sigma**2)) = 1/2
        #  -HWHM**2/(2*sigma**2) = ln(1/2)
        #  HWHM**2/(2*sigma**2) = ln(2)
        #  HWHM**2 = ln(2) * (2*sigma**2)
        #  HWHM = sqrt(ln(2)) * sqrt(2) * sigma
        #  FWHM = 2*HWHM = 2*sqrt(2)*sqrt(ln(2)) * sigma
        #ax[indx].plot(x,gaussian(x,*popt),'r:',label=rf'Gauss fit FWHM = {2*np.sqrt(2)*np.sqrt(np.log(2))*sigma:.2f}')
        ax[indx].axvline(mean, color='b', linestyle=':')#, label=f"mean {round(mean,3)}, std {round(sta_dev,3)}")
        # Reset the xlim
        ax[indx].set_xlim([-lim, lim])
        ax[indx].set_ylim([0,max(hist_y)+max(hist_y)/20])
        ax[indx].annotate(f"Mean: {round(mean,3)}",xy=(-lim+lim/20,max(hist_y)-max(hist_y)/10), fontsize=TINY_SIZE) 
        ax[indx].annotate(f"Stdev: {round(sta_dev,2)}",xy=(-lim+lim/20,max(hist_y)-max(hist_y)/5), fontsize=TINY_SIZE)

        #ax[indx].legend(loc='upper center')
        ax[indx].set_xlabel(f"{QUANTITY_FANCY_NAMES[QUANTITY.strip('_rel_diff')]} difference [%]",fontsize=22)
        ax[indx].set_ylabel("Count",fontsize=SMALL_SIZE)
        ax[indx].tick_params(axis="x",labelsize=SMALL_SIZE)
        ax[indx].tick_params(axis="y",labelsize=SMALL_SIZE)
        #set(xlabel=f"{DEFAULT_PREFACTOR}*{QUANTITY}", ylabel='Frequency', Fontsize=30)
        
    fig.suptitle(f"FLEUR vs WIEN2K")
    fig.tight_layout()
    fig.savefig(f"{name_file}")
    pl.close(fig)
