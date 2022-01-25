#!/usr/bin/env python
import json
import os
import sys

# UPDATE if you change the output data format
DATA_VERSION = '0.0.1'

# TODO: replace with final list of elements
ALL_ELEMENTS = ['Si', 'Ba', 'Au']

def get_plugin_name():
    file_name = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        os.pardir, os.pardir, os.pardir, 'plugin_name.txt'
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

def get_num_atoms(key):
    conf_natoms_map = {'XO': 2, 'XO2': 3, 'XO3': 4, 'X2O': 3, 'X2O3': 10, 'X2O5': 14}
    # NOTE! This for now just works for the few elements in the monoelemental test set
    # Needs to be udpated, or even better queried from the AiiDA group (and/or exported
    # to JSON?)
    monoelemental_natoms_map = {'O': 8, 'Si': 2, 'Ba': 1, 'Au': 1}

    element, configuration = key.split('-')
    if element == 'O' or configuration == 'X':
        return monoelemental_natoms_map[element]
    else:
        return conf_natoms_map[configuration]

def get_O_percentage(key):
    conf_percentage_map = {'XO': 1/2, 'XO2': 2/3, 'XO3': 3/4, 'X2O': 1/3, 'X2O3': 3/5, 'X2O5': 5/7}
    element, configuration = key.split('-')

    if element == 'O':
        return 1.
    if configuration == 'X':
        return 0.
    return conf_percentage_map[configuration] 

def get_formation_energy(energy0, energy1, percentage, energy):
    # Linear interpolation of energy0 and enery1 at the `percentage` concentration
    # NOTE: no corrections applied!
    ref_energy = energy0 + percentage * (energy1 - energy0)
    return energy - ref_energy

if __name__ == "__main__":
    try:
        with open(f'../monoelemental-results-{PLUGIN_NAME}.json') as fhandle:
            monoelemental = json.load(fhandle)

        with open(f'../results-{PLUGIN_NAME}.json') as fhandle:
            oxides = json.load(fhandle)
    except OSError as exc:
        print(f"Error, probably you didn't generate the file yet! Original error: {exc}")
        sys.exit(1)

    all_data = {}
    for element in ALL_ELEMENTS:
        key = f'{element}-X'
        E0_element = monoelemental['BM_fit_data'][key]['E0'] / get_num_atoms(key)

        # TODO: replace this if the configuration for oxygen is not the other element,
        # but another string
        oxygen_key = f"O-{element}"
        E0_oxygen = monoelemental['BM_fit_data'][oxygen_key]['E0'] / get_num_atoms("O-X")

        all_data[element] = {
            'O': E0_oxygen,
            'element': E0_element,
            'oxides': {}
        }

        for key in oxides['BM_fit_data'].keys():
            if not key.startswith(f"{element}-"):
                continue

            E0_oxide = oxides['BM_fit_data'][key]['E0'] / get_num_atoms(key)
            percentage_oxide = get_O_percentage(key)
            element, configuration = key.split('-')
            all_data[element]['oxides'][configuration] = [percentage_oxide, E0_oxide]

    all_data_flat = {}
    all_convex_hull = {}
    for element in all_data:
        all_data_flat[element] = [
            [0, 0, element],
            [1, 0, 'O']
        ]
        for configuration, (percentage_oxide, E0_oxide) in all_data[element]['oxides'].items():
            all_data_flat[element].append([
                percentage_oxide,
                get_formation_energy(all_data[element]['element'], all_data[element]['O'], percentage_oxide, E0_oxide),
                configuration.replace('X', element)
            ])
        all_data_flat[element].sort()

        # Check the endpoints of the convex hull
        last_idx = 0
        on_convex_hull = [last_idx]
        while last_idx < len(all_data_flat[element]) - 1:
            tmp_slopes = []
            lastx, lasty, lastlabel = all_data_flat[element][last_idx]
            for idx in range(last_idx + 1, len(all_data_flat[element])):
                x, y, label = all_data_flat[element][idx]
                tmp_slopes.append(((y - lasty) / (x - lastx), idx))
            tmp_slopes.sort()
            # Append the idx of the first element, i.e. the one with the smallest slope
            newidx = tmp_slopes[0][1]
            on_convex_hull.append(newidx)
            last_idx = newidx

        above_hull = {}
        # There are always at least two (the two endpoints),
        # and the first is always the first point (idx=0)
        previous_on_hull = 0
        next_convex_hull_idx = 1
        #next_on_hull = on_convex_hull[next_convex_hull_idx]
        for idx in range(1, len(all_data_flat[element])):
            if idx in on_convex_hull:
                previous_on_hull = idx
                next_convex_hull_idx += 1
                # I don't do it here as it would fail for the very last point
                #next_on_hull = on_convex_hull[next_convex_hull_idx]
            else:
                # add a new entry, the key is the label, the values are: the
                # idx, and the labels of the previous and next
                next_on_hull = on_convex_hull[next_convex_hull_idx]
                previous_on_hull_label = all_data_flat[element][previous_on_hull][2]
                next_on_hull_label = all_data_flat[element][next_on_hull][2]
                above_hull[all_data_flat[element][idx][2]] = (
                    idx,previous_on_hull_label, next_on_hull_label)

        all_convex_hull[element] = {
            'on_hull': {all_data_flat[element][idx][2]: idx for idx in on_convex_hull},
            'above_hull': above_hull,
            }

    # Plotting now
    OUT_FOLDER = 'convex-hulls'
    os.makedirs(OUT_FOLDER, exist_ok=True)

    fname = f'{OUT_FOLDER}/convex-hull-data-{PLUGIN_NAME}.json'
    with open(fname, 'w') as fhandle:
        json.dump({
            'all_data': all_data,
            'all_data_flat': all_data_flat,
            'all_convex_hull': all_convex_hull,
            'data_version': DATA_VERSION
        }, fhandle)
    print(f"File '{fname}' written.")
