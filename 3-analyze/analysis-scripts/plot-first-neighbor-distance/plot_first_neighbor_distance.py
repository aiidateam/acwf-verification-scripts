#!/usr/bin/env python
import json
import sys
from collections import defaultdict

from ase.data import chemical_symbols, atomic_numbers
import numpy as np
import pylab as pl

Z_max = 96
Pettifor_max = 103

# Precompute Pettifor values
from pymatgen.core.periodic_table import Element
pettifor_scale = {}
# I can create the Pettifor 
for Z in range(1, Pettifor_max + 1):
    pettifor = int(Element(chemical_symbols[Z]).mendeleev_no)
    pettifor_scale[chemical_symbols[Z]] = pettifor
assert len(pettifor_scale.values()) == len(set(pettifor_scale.values())), "Duplicate Pettifor values found!"

# Generate also reverse table, 
inverse_pettifor_to_Z = [None] * (max(pettifor_scale.values()) + 1)
for symbol, pettifor in pettifor_scale.items():
    inverse_pettifor_to_Z[pettifor] = atomic_numbers[symbol]

#print(inverse_pettifor_to_Z)

alat_to_first_neighbor_factor = {
    'SC': 1.,
    'FCC': np.sqrt(2)/2,
    'BCC': np.sqrt(3)/2,
    'Diamond': np.sqrt(3)/4,
    'XO': 1/2,
    'X2O': np.sqrt(3)/4,
    'XO3': 1/2,
    'XO2': np.sqrt(3)/4,
    'X2O3': np.sqrt(3)/4,
    'X2O5': np.sqrt(3)/4
}

MARKERS = {
    'SC': 's',
    'FCC': 'o',
    'BCC': 'X',
    'Diamond': '^',
    'XO': 's',
    'X2O': 'o',
    'XO3': 'X',
    'XO2': '^',
    'X2O3': 'D',
    'X2O5': 'P',
}

# IMPORTANT: update correct path!
print("*"*72)
print("*   IMPORTANT!! PLEASE UPDATE THE REFERENCE DATA FOR THE FINAL PLOT!   *")
print("*"*72)
FLEUR_UNARIES_FILE = "../temp-results/results-unaries-set2-fleur-k_0.06.json"
WIEN2K_UNARIES_FILE = "../temp-results/results-unaries-set2-wien2k-dk_0.06.json"
FLEUR_OXIDES_FILE = "../temp-results/results-set2plusoxygen-fleur-temp-v003.json"
WIEN2K_OXIDES_FILE = "../temp-results/results-set2plusoxygen-wien2k-temp-v003.json"

# Essentially, how many atoms there are in the conventional cube
volume_per_atom_to_cubic_volume = {
    'SC': 1.,
    'BCC': 2.,
    'FCC': 4.,
    'Diamond': 8.,
    'X2O': 12,
    'XO': 8,
    'XO2': 12,
    'X2O3': 10,
    'X2O5': 14,
    'XO3': 4,
}

def get_alat_from_raw_json(json_data):
    assert json_data['script_version'] == "0.0.3"

    data = defaultdict(dict)

    for key, values in json_data['BM_fit_data'].items():
        element, config = key.split('-')
        if config.startswith('X/'):
            config = config[len('X/'):]

        volume = values['min_volume']
        num_atoms_in_sim_cell = json_data['num_atoms_in_sim_cell'][key]
        volume_per_atom = volume / num_atoms_in_sim_cell
        cubic_volume = volume_per_atom * volume_per_atom_to_cubic_volume[config]
        data[config][element] = cubic_volume ** (1/3)

    return dict(data)

def generate_plots(fleur_alats, wien2k_alats, plot_vs_pettifor=False):

    # We compute the average between fleur and wien2k
    # First, check that there are all the same configurations and elements
    assert fleur_alats.keys() == wien2k_alats.keys()
    for key in fleur_alats:
        assert fleur_alats[key].keys() == wien2k_alats[key].keys()
    # Now we generate the average, and do some double checks
    average_alats = defaultdict(dict)
    for config in fleur_alats:
        for element in fleur_alats[config]:
            fleur = fleur_alats[config][element]
            wien2k = wien2k_alats[config][element]
            average_alats[config][element] = (fleur + wien2k) / 2
            assert abs((fleur - wien2k) / wien2k) < 0.01, f"Data for {element}-{config} has large error: {fleur} {wien2k} {abs((fleur - wien2k) / wien2k)}!"
            assert abs((fleur - wien2k) / wien2k) > 1.e-14, f"Data for {element}-{config} seem to be really identical! Maybe a copy-paste error?"

    if 'Diamond' in average_alats:
        set_type = 'unaries'
        valid_configurations = ['SC', 'BCC', 'FCC', 'Diamond'] # I hardcode them here to have a given order
    elif 'X2O5' in average_alats:
        set_type = 'oxides'
        valid_configurations = ['X2O', 'XO', 'X2O3', 'XO2', 'X2O5', 'XO3'] # I hardcode them here to have a given order
    else:
        raise ValueError("Unknown set!")
    assert set(valid_configurations) == set(average_alats.keys())

    if plot_vs_pettifor:
        pettifor_values = [True, False]
    else:
        pettifor_values = [False]
    for plot_vs_pettifor in pettifor_values:
        fig, ax = pl.subplots(figsize=(9, 3))
        pl.subplots_adjust(left=0.07, right=0.99, bottom=0.15)
        if plot_vs_pettifor:
            pl.xlabel("Mendeleev's number")
            x = np.arange(1, max(pettifor_scale.values()) + 1)
        else:
            pl.xlabel("Atomic number $Z$")
            x = np.arange(1, Z_max + 1)
            y = np.zeros(len(x))
            for conf in valid_configurations: # I use this dictionary so they are in the order I want
                for idx in range(len(x)):
                    if plot_vs_pettifor:
                        pettifor = x[idx]
                        Z = inverse_pettifor_to_Z[pettifor]
                    else:
                        Z = x[idx]

                    # In case Z is None, I pre-set to None so the point is not plotted
                    first_neighbor = None
                    if Z <= Z_max:
                        # alat * factor
                        first_neighbor = average_alats[conf][chemical_symbols[Z]] * alat_to_first_neighbor_factor[conf]
                    else:
                        print(idx, Z, 'SKIP')
                    y[idx] = first_neighbor

                marker = MARKERS[conf]

                latex_conf = "".join([f"$_{char}$" if char in "0123456789" else char for char in conf])
                pl.plot(x, y, f'{marker}-', markersize=3, linewidth=1, label=f'{latex_conf}')

            pl.ylabel("First-neighbor distance (Å)")

            from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                                        AutoMinorLocator)

            ax.xaxis.set_major_locator(MultipleLocator(10))
            ax.minorticks_on()
            ax.xaxis.set_minor_locator(MultipleLocator(1))
            pl.grid(which='major', axis='x', color='#ccc', linestyle='-')
            pl.grid(which='minor', axis='x', color='#eee', linestyle='-')
            pl.xlim(1,x.max())

            pl.grid(which='major', axis='y', color='#ccc', linestyle='-')
            pl.grid(which='minor', axis='y', color='#eee', linestyle='dotted')

            sec = ax.secondary_xaxis(location='top')
            if plot_vs_pettifor:
                ticks_x = np.arange(1, max(pettifor_scale.values()) + 1)
                sec.set_xticks(ticks_x)
                sec.set_xticklabels([chemical_symbols[inverse_pettifor_to_Z[pettifor]] for pettifor in ticks_x])
            else:
                ticks_x = np.arange(1, Z_max+1)
                sec.set_xticks(ticks_x)
                sec.set_xticklabels([chemical_symbols[Z] for Z in ticks_x], fontsize=7)

            sec.tick_params(rotation=90)

            if plot_vs_pettifor:
                pl.legend(loc='upper right', ncol=2)
            else:
                if set_type == 'oxides':
                    pl.ylim(1.4, 3.7)
                    pl.legend(loc='lower right', ncol=6)
                else:
                    pl.ylim(1, 5.8)
                    pl.legend(loc='lower right', ncol=4)

        fname_suffix = "-vs_mendeleev.pdf" if plot_vs_pettifor else ""
        fname = f'first-neighbor-distance-{set_type}{fname_suffix}.pdf'
        pl.savefig(fname)
        print(f"'{fname}' written.")


if __name__ == "__main__":
    if sys.argv[1:2] == ['unaries']:
        with open(FLEUR_UNARIES_FILE) as fhandle:
            fleur_data = json.load(fhandle)
        with open(WIEN2K_UNARIES_FILE) as fhandle:
            wien2k_data = json.load(fhandle)
    elif sys.argv[1:2] == ['oxides']:
        with open(FLEUR_OXIDES_FILE) as fhandle:
            fleur_data = json.load(fhandle)
        with open(WIEN2K_OXIDES_FILE) as fhandle:
            wien2k_data = json.load(fhandle)
    else:
        print("Pass either 'oxides' or 'unaries' on the command line.")
        sys.exit(2)

    fleur_alats = get_alat_from_raw_json(fleur_data)
    wien2k_alats = get_alat_from_raw_json(wien2k_data)

    generate_plots(fleur_alats, wien2k_alats, plot_vs_pettifor=False)

