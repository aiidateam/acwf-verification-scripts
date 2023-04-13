#!/usr/bin/env python
import json
import os
import sys

import numpy as np
import pylab as pl
import tqdm

from acwf_paper_plots.eosfit_31_adapted import BM, echarge
from acwf_paper_plots.quantities_for_comparison import birch_murnaghan

EXPECTED_SCRIPT_VERSION = '0.0.3'
RESIDUALS_THRESHOLD = 1.e-3


DATA_FOLDER = "../../../code-data"
with open(os.path.join(DATA_FOLDER, "labels.json")) as fhandle:
    labels_data = json.load(fhandle)
FLEUR_LABEL = labels_data['all-electron-keys']['FLEUR']

# ONLY_THESE_SYSTEMS = None # Run for all
ONLY_THESE_SYSTEMS = [
    ('Er', 'X/Diamond'),    
    ('Cs', 'XO2'),
    ('Rb', 'XO3')
]


def fit_eos_data(eos_data):
    # I need to pass a numpy array
    min_volume, E0, bulk_modulus_internal, bulk_deriv, residuals = BM(np.array(eos_data))
    bulk_modulus_GPa = bulk_modulus_internal * echarge * 1.0e21
    #1 eV/Angstrom3 = 160.21766208 GPa
    bulk_modulus_ev_ang3 = bulk_modulus_GPa / 160.21766208
    BM_fit_data = {
        'min_volume': min_volume,
        'E0': E0,
        'bulk_modulus_ev_ang3': bulk_modulus_ev_ang3,
        'bulk_deriv': bulk_deriv,
        'residuals': residuals[0]
    }
    return BM_fit_data

def get_conf_nice(configuration_string):
    """Convert the configuration string to a nicely typeset string in LaTeX."""
    ret_pieces = []
    for char in configuration_string:
        if char in "0123456789":
            ret_pieces.append(f"$_{char}$")
        else:
            ret_pieces.append(char)
    return "".join(ret_pieces)


def plot(SET_NAME):
    # Load the TS data provided by FLEUR
    try:
        with open(f'TS_data_fleur/ts_contributions-{SET_NAME}-verification-PBE-v1-fleur.json') as fhandle:
            reference_plugin_data = json.load(fhandle)
    except OSError:
        print(f"No TS data found for FLEUR")
        sys.exit(1)
    
    if not reference_plugin_data['script_version'] == EXPECTED_SCRIPT_VERSION:
        raise ValueError(
            f"This script only works with data generated at version {EXPECTED_SCRIPT_VERSION}. "
            "Please re-run ./get_results.py to update the data format!")

    PLOT_FOLDER = f'TS-plots-{SET_NAME}-fleur'
    os.makedirs(PLOT_FOLDER, exist_ok=True)

    F_keys = set(reference_plugin_data['free_energy_data'].keys())
    TS_keys = set(reference_plugin_data['ts_contribution'].keys())
    all_systems = F_keys.intersection(TS_keys)

    progress_bar = tqdm.tqdm(sorted(all_systems))
    for element_and_configuration in progress_bar:
        progress_bar.set_description(f"{element_and_configuration:12s}")
        progress_bar.refresh()

        element, configuration = element_and_configuration.split('-')

        if (element, configuration) not in ONLY_THESE_SYSTEMS:
            continue

        F_data = np.array(reference_plugin_data['free_energy_data'][f'{element}-{configuration}'])
        TS_data = np.array(reference_plugin_data['ts_contribution'][f'{element}-{configuration}'])

        # Get the x axis for the plot
        volumes, free_energies = F_data.T
        volumes_TS, TS_contrib = TS_data.T

        assert all([diff < 1.e-6 for diff in np.abs(volumes - volumes_TS)]), 'Volumes from free energy and from TS do not match'

        dense_volumes = np.linspace(
            min(volumes),
            max(volumes),
            100
        )


        # FIT CURVES!
        BM_fit_data_free_energy = fit_eos_data(np.array([volumes, free_energies]).T)
        BM_fit_data_E = fit_eos_data(np.array([volumes, free_energies + TS_contrib]).T) # I cannot do F_data + TS_data, it would also sum the volumes
        BM_fit_data_E_minus_TS_half = fit_eos_data(np.array([volumes, free_energies + TS_contrib / 2]).T)

        fitted_free_energy = birch_murnaghan(
             V=dense_volumes,
             E0 =BM_fit_data_free_energy['E0'],
             V0 =BM_fit_data_free_energy['min_volume'],
             B0 =BM_fit_data_free_energy['bulk_modulus_ev_ang3'],
             B01=BM_fit_data_free_energy['bulk_deriv']
        )
        residuals_free_energy = BM_fit_data_free_energy['residuals']
        fitted_E = birch_murnaghan(
             V=dense_volumes,
             E0 =BM_fit_data_E['E0'],
             V0 =BM_fit_data_E['min_volume'],
             B0 =BM_fit_data_E['bulk_modulus_ev_ang3'],
             B01=BM_fit_data_E['bulk_deriv']
        )
        residuals_E = BM_fit_data_E['residuals']
        fitted_E_minus_TS_half = birch_murnaghan(
             V=dense_volumes,
             E0 =BM_fit_data_E_minus_TS_half['E0'],
             V0 =BM_fit_data_E_minus_TS_half['min_volume'],
             B0 =BM_fit_data_E_minus_TS_half['bulk_modulus_ev_ang3'],
             B01=BM_fit_data_E_minus_TS_half['bulk_deriv']
        )
        residuals_E_minus_TS_half = BM_fit_data_E_minus_TS_half['residuals']

        fig, eos_ax = pl.subplots(nrows=1, ncols=1, gridspec_kw={'left': 0.15, 'right': 0.95})

        # Plot EOS: this will be done anyway
        eos_ax.plot(volumes, free_energies - BM_fit_data_free_energy['E0'], 'ob', label=f'E-TS')
        eos_ax.plot(dense_volumes, fitted_free_energy - BM_fit_data_free_energy['E0'], '-b', label=f'E-TS fit')# (residuals: {residuals_free_energy:.3g})')
        eos_ax.axvline(BM_fit_data_free_energy['min_volume'], linestyle='--', color='b', alpha=0.5)

        eos_ax.plot(volumes, free_energies + TS_contrib / 2  - BM_fit_data_E_minus_TS_half['E0'], 'xg', label=f'E - TS/2')
        eos_ax.plot(dense_volumes, fitted_E_minus_TS_half - BM_fit_data_E_minus_TS_half['E0'], '-g', label=f'E - TS/2 fit')# (residuals: {residuals_E_minus_TS_half:.3g})')
        eos_ax.axvline(BM_fit_data_E_minus_TS_half['min_volume'], linestyle='--', color='g', alpha=0.5)

        eos_ax.plot(volumes, free_energies + TS_contrib - BM_fit_data_E['E0'], 'sr', label=f'E')
        eos_ax.plot(dense_volumes, fitted_E - BM_fit_data_E['E0'], '-r', label=f'E fit')# (residuals: {residuals_E:.3g})')
        eos_ax.axvline(BM_fit_data_E['min_volume'], linestyle='--', color='r', alpha=0.5)
        
        eos_ax.legend(loc='upper center')
        eos_ax.set_xlabel("Cell volume per simulation cell ($\\AA^3$)")
        eos_ax.set_ylabel("$E - TS$ per simulation cell (eV)")

        #LIGHTYELLOW = (255/255, 244/255, 214/255)
        #LIGHTORANGE = (255/255, 205/255, 171/255)
        #if residuals is None:
        #    eos_ax.set_facecolor(LIGHTYELLOW)
        #elif residuals > RESIDUALS_THRESHOLD:
        #    eos_ax.set_facecolor(LIGHTORANGE)

        conf_nice = get_conf_nice(configuration)
        fig.suptitle(f"FLEUR - {element} ({conf_nice})")

        pl.savefig(f"{PLOT_FOLDER}/{element}-{configuration.replace('/', '_')}.pdf")
        pl.close(fig)

    print(f"Plots written to: '{PLOT_FOLDER}'")

if __name__ == "__main__":
    plot("unaries")
    plot("oxides")
