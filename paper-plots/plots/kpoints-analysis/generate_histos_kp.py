#!/usr/bin/env python
import json
import os
import sys

import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
import tqdm

import quantities_for_comparison as qc


BINS = 100
DEFAULT_PREFACTOR = 100
DEFAULT_wb0 = 1.0/8.0
DEFAULT_wb1 = 1.0/64.0
LIMITS = {"V0_rel_diff":0.03,"B0_rel_diff":0.25,"B1_rel_diff":2.5}


quantity_for_comparison_map = {
    "B0_rel_diff": qc.B0_rel_diff,
    "V0_rel_diff": qc.V0_rel_diff,
    "B1_rel_diff": qc.B1_rel_diff,
}


if __name__ == "__main__":
   

    for set_name in ['unaries','oxides']:
        
        with open(f'../../code-data/results-{set_name}-verification-PBE-v1-wien2k-dk_0.06.json') as fhandle:
            reference_plugin_data = json.load(fhandle)

        compare_plugin_data = []
        with open(f'results-{set_name}-0045.json') as fhandle:
            compare_plugin_data.append(json.load(fhandle))

        name_file = f'histo-{set_name}.pdf'

        all_systems = set(reference_plugin_data['eos_data'].keys())
        all_systems = set(reference_plugin_data['BM_fit_data'].keys())

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

            print(f'{set_name} {QUANTITY}')

            for index, compare_plugin in enumerate(compare_plugin_data):

                collect = []

                progress_bar = tqdm.tqdm(sorted(all_systems))
                
                for element_and_configuration in progress_bar:
                    #progress_bar.set_description(f"{element_and_configuration:12s}")
                    #progress_bar.refresh()

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

                    quant = quantity_for_comparison_map[QUANTITY](V0,B0,B01,CV0,CB0,CB01,DEFAULT_PREFACTOR,DEFAULT_wb0,DEFAULT_wb1)

                    if quant > LIMITS[QUANTITY]:
                        print(f'{element}-{configuration} {quant}')
                    if quant < -LIMITS[QUANTITY]:
                        print(f'{element}-{configuration} {quant}')

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
                    ax[indx].annotate(f"+{countBig}", xy=(lim, max(hist_y)/2), xytext=(0.5*lim, max(hist_y)/2), arrowprops=dict(facecolor='black', shrink=0.05))
                if countSmall:
                    ax[indx].annotate(f"+{countSmall}", xy=(-lim, max(hist_y)/2), xytext=(-lim*0.8, max(hist_y)/2), arrowprops=dict(facecolor='black', shrink=0.05))
                
                ax[indx].axvline(mean, color='b', linestyle=':')#, label=f"mean {round(mean,3)}, std {round(sta_dev,3)}")
                
                # Reset the xlim
                ax[indx].set_xlim([-lim, lim])
                ax[indx].set_ylim([0,max(hist_y)+max(hist_y)/20])
                ax[indx].annotate(f"mean {round(mean,3)}",xy=(-lim+lim/20,max(hist_y)-max(hist_y)/10)) 
                ax[indx].annotate(f"std {round(sta_dev,2)}",xy=(-lim+lim/20,max(hist_y)-max(hist_y)/5))

                #ax[indx].legend(loc='upper center')
                ax[indx].set_xlabel(f"% difference in {QUANTITY.strip('_rel_diff')}",fontsize=22)
                ax[indx].set_ylabel("Frequency",fontsize=22)
                ax[indx].tick_params(axis="x",labelsize=20)
                ax[indx].tick_params(axis="y",labelsize=20)
                #set(xlabel=f"{DEFAULT_PREFACTOR}*{QUANTITY}", ylabel='Frequency', Fontsize=30)
        
        fig.suptitle(f"WIEN2K kpoints mesh 0.06 vs 0.045 {set_name}")
        fig.tight_layout()
        fig.savefig(f"{name_file}")
        pl.close(fig)
