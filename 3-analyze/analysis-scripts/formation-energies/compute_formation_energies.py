#!/usr/bin/env python
import json
import os
import re
import sys

import ase.data
import numpy as np

# TODO: add back oxygen
ALL_ELEMENTS = [ase.data.chemical_symbols[Z] for Z in range(1, 96+1) if Z != 8]

def get_num_atoms(configuration):
    conf_natoms_map = {'XO': 2, 'XO2': 3, 'XO3': 4, 'X2O': 3, 'X2O3': 10, 'X2O5': 14}
    return conf_natoms_map[configuration]

def get_ratio_to_formula_unit(configuration):
    conf_map = {'XO': 1, 'XO2': 1, 'XO3': 1, 'X2O': 1, 'X2O3': 2, 'X2O5': 2}
    return conf_map[configuration]

def get_num_atoms_in_formula_unit(configuration):
    num_atoms_in_cell = get_num_atoms(configuration)
    ratio = get_ratio_to_formula_unit(configuration)
    assert num_atoms_in_cell % ratio == 0
    return num_atoms_in_cell // ratio

def get_X_O_in_formula_unit(configuration):
    map = {'XO': (1, 1), 'XO2': (1, 2), 'XO3': (1, 3), 'X2O': (2, 1), 'X2O3': (2, 3), 'X2O5': (2, 5)}
    return map[configuration] 

def get_O_percentage(configuration):
    X, O = get_X_O_in_formula_unit(configuration)
    return O / (X + O)

def get_formation_energy(element, configuration_triple, plugin):
    """
    Given an element and three configurations to compare, compute the formation energy 
    of the central configuration w.r.t. to the two extremal ones.

    You need to pass the plugin data, from which E0 is obtained.
    """    
    assert len(configuration_triple) == 3

    O_percentages = [get_O_percentage(conf) for conf in configuration_triple]
    assert O_percentages[0] < O_percentages[2], "percentage_left must be < percentage_right"
    assert O_percentages[1] >= O_percentages[0], "percentage must be >= percentage_left"
    assert O_percentages[1] <= O_percentages[2], "percentage must be <= percentage_right"

    try:
        energies_per_formula_unit = [
            plugin['BM_fit_data'][f'{element}-{configuration}']['E0'] / get_ratio_to_formula_unit(configuration) for configuration in configuration_triple
        ]
    except (KeyError, TypeError): # Something is None because of a fit problem, of wasn't even computed
        return None

    X_L, O_L = get_X_O_in_formula_unit(configuration_triple[0])
    X_C, O_C = get_X_O_in_formula_unit(configuration_triple[1])
    X_R, O_R = get_X_O_in_formula_unit(configuration_triple[2])

    # I want to find the linear combination of the two endpoints that gives me the central point
    # so e.g. if I have X2O, XO and X2O3 as the three points, for X I want to find alpha (to multiply to X2O)
    # and beta (to multiply to X2O3) that gives me XO, i.e. for X:
    # alpha * 2 + beta * 2 = 1
    # and for O:
    # alpha * 1 + beta * 3 = 1
    # So we need to solve this linear system (_L, _R and _C for left, right and center):
    # 
    #  X_L  X_R  |  X_C
    #  O_L  O_R  |  O_C
    # I call A the 2x2 matrix, b the vector
    # The solution is given by (A^-1) @ b
    A = np.array([[X_L, X_R], [O_L, O_R]])
    b = np.array([X_C, O_C]).T
    alpha, beta = np.linalg.inv(A) @ b

    # Now we can do a linear interpolation of the energies and compare with the central energy
    # Note that this is per formula unit, so we also need to divide by the number of atoms per formula unit in the center
    reference_energy = alpha * energies_per_formula_unit[0] + beta * energies_per_formula_unit[2]
    formation_per_formula_unit = energies_per_formula_unit[1] - reference_energy
    return formation_per_formula_unit / get_num_atoms_in_formula_unit(configuration_triple[1])

def generate_json_data(plugin_name_1, plugin_name_2):

    try:
        with open(f'../results-{plugin_name_1}.json') as fhandle:
            plugin1 = json.load(fhandle)
        with open(f'../results-{plugin_name_2}.json') as fhandle:
            plugin2 = json.load(fhandle)
    except OSError as exc:
        print(f"Error, probably you didn't generate the file yet! Original error: {exc}")
        sys.exit(1)

    all_data = {}

    for element in ALL_ELEMENTS:
        for configuration_triple in [
                ('X2O', 'XO', 'XO3'),
                ('X2O', 'XO2', 'XO3'),
                ('X2O', 'X2O3', 'XO3'),
                ('X2O', 'X2O5', 'XO3'),
            ]:
            #percentage_left = get_O_percentage(configuration_triple[0])
            #percentage = get_O_percentage(configuration_triple[1])
            #percentage_right = get_O_percentage(configuration_triple[2])

            plugin1_formation_energy = get_formation_energy(element, configuration_triple, plugin1)
            plugin2_formation_energy = get_formation_energy(element, configuration_triple, plugin2)

            all_data[f"{element}-{configuration_triple[0]}/{configuration_triple[1]}/{configuration_triple[2]}"] = {
                plugin_name_1: plugin1_formation_energy,
                plugin_name_2: plugin2_formation_energy
            }

    # Plotting now
    OUT_FOLDER = 'formation-energies-output'
    os.makedirs(OUT_FOLDER, exist_ok=True)

    fname = f'{OUT_FOLDER}/formation-energies-{plugin_name_1}-VS-{plugin_name_2}.json'
    with open(fname, 'w') as fhandle:
        json.dump({'formation_energies': all_data}, fhandle, sort_keys=True, indent=2)
    print(f"File '{fname}' written.")

if __name__ == "__main__":
    try:
        plugin1 = sys.argv[1]
        plugin2 = sys.argv[2]
    except IndexError:
        print("Pass as two parameters the two plugins to compare")
        sys.exit(1)
    generate_json_data(plugin1, plugin2)
