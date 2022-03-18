#!/usr/bin/env python
import json
import os
import re
import sys

import ase.data
import numpy as np

UNARIES_CONFIGURATIONS = ['X/BCC', 'X/SC', 'X/FCC', 'X/Diamond']
OXIDES_CONFIGURATIONS = ['XO', 'XO2', 'XO3', 'X2O', 'X2O3', 'X2O5']
ALL_ELEMENTS = [ase.data.chemical_symbols[Z] for Z in range(1, 96+1)]
VERBOSE = False

def get_num_atoms_in_formula_unit(configuration):
    return sum(get_X_O_in_formula_unit(configuration))

def get_X_O_in_formula_unit(configuration):
    map = {
        'XO': (1, 1), 'XO2': (1, 2), 'XO3': (1, 3),
        'X2O': (2, 1), 'X2O3': (2, 3), 'X2O5': (2, 5),
        # For the unaries I could put any value, even (2, 0) or similar; anyway in the end
        # everything is renormalized per atom. It is however important that they are all the same,
        # to properly deal with the cases below ('X' meaning lowest-in-energy X/..., and O meaning
        # lowest-in-energy X/... for element 'X=oxygen')
        'X/BCC': (1, 0), 'X/SC': (1, 0), 'X/FCC': (1, 0), 'X/Diamond': (1, 0),
        # Deal with the special cases
        'X': (1, 0), 'O': (0, 1)
        }
    return map[configuration]

def get_O_percentage(configuration):
    X, O = get_X_O_in_formula_unit(configuration)
    return O / (X + O)

def _get_energy_internal(plugin, element, configuration):
    system = f'{element}-{configuration}'
    if plugin['BM_fit_data'][system] is None:
        return None
    return plugin['BM_fit_data'][system]['E0'] / plugin['num_atoms_in_sim_cell'][system] * get_num_atoms_in_formula_unit(configuration)


def get_energy_per_formula_unit(plugin_unaries, plugin_oxides, element, configuration, unaries_energy_diff=None):
    if configuration in UNARIES_CONFIGURATIONS + ['X', 'O']:
        plugin = plugin_unaries
    elif configuration in OXIDES_CONFIGURATIONS:
        plugin = plugin_oxides
    else:
        raise ValueError("Unknown configuration string")
    
    if configuration == 'X':
        if unaries_energy_diff is None:
            raise ValueError("You need to pass `unaries_energy_diff` if one of the configurations is 'X'")
        # Enforce comparing with the energy of the reference element
        return _get_energy_internal(plugin, element, unaries_energy_diff[element]['reference_configuration'])
    elif configuration == 'O':
        if unaries_energy_diff is None:
            raise ValueError("You need to pass `unaries_energy_diff` if one of the configurations is 'X'")
        if unaries_energy_diff is None:
            raise ValueError("You need to pass `unaries_energy_diff` if one of the configurations is 'X'")
        # Enforce comparing with the energy of oxygen, independent of the element
        return _get_energy_internal(plugin, 'O', unaries_energy_diff['O']['reference_configuration'])

    else:
        return _get_energy_internal(plugin, element, configuration)


def get_formation_energy(element, configuration_triple, plugin_unaries, plugin_oxides, unaries_energy_diff):
    """
    Given an element and three configurations to compare, compute the formation energy per atom
    of the central configuration w.r.t. to the two extremal ones.

    You need to pass the two JSON-loaded data from the plugin (for unaries and for oxides),
    from which E0 is extracted.
    """    
    assert len(configuration_triple) == 3

    assert plugin_unaries['script_version'] == '0.0.3'
    assert plugin_oxides['script_version'] == '0.0.3'

    O_percentages = [get_O_percentage(conf) for conf in configuration_triple]
    assert O_percentages[0] < O_percentages[2], "percentage_left must be <= percentage_right"
    assert O_percentages[1] >= O_percentages[0], "percentage must be >= percentage_left"
    assert O_percentages[1] <= O_percentages[2], "percentage must be <= percentage_right"

    #try:
    energies_per_formula_unit = [
        get_energy_per_formula_unit(plugin_unaries, plugin_oxides, element, configuration, unaries_energy_diff)
        for configuration in configuration_triple
    ]
    #except (KeyError, TypeError): # Something is None because of a fit problem, of wasn't even computed
    #    return None

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
    if any(energy is None for energy in energies_per_formula_unit):
        return None
    reference_energy = alpha * energies_per_formula_unit[0] + beta * energies_per_formula_unit[2]
    formation_per_formula_unit = energies_per_formula_unit[1] - reference_energy
    return formation_per_formula_unit / get_num_atoms_in_formula_unit(configuration_triple[1])

def get_unaries_energy_difference(
    plugin_name_1, plugin_name_2, plugin1_unaries, plugin2_unaries):
    unaries_data = {}
    for element in ALL_ELEMENTS:
        plugin1_energies = []
        plugin2_energies = []
        # TODO: fix when there are missing elements
        for configuration in UNARIES_CONFIGURATIONS:
            # I pass None for the plugin_oxides since they should not be used
            plugin1_energies.append(
                get_energy_per_formula_unit(plugin1_unaries, None, element, configuration)
            )
            plugin2_energies.append(
                get_energy_per_formula_unit(plugin2_unaries, None, element, configuration)
            )
        plugin1_min_energy_conf = UNARIES_CONFIGURATIONS[np.argmin(plugin1_energies)]
        plugin2_min_energy_conf = UNARIES_CONFIGURATIONS[np.argmin(plugin2_energies)]
        # I use the first plugin as reference
        ref_min_energy_conf = plugin1_min_energy_conf

        plugin1_energy_dict = dict(zip(UNARIES_CONFIGURATIONS, plugin1_energies))
        plugin2_energy_dict = dict(zip(UNARIES_CONFIGURATIONS, plugin2_energies))
        if VERBOSE:
            print(f"{element:2s}: {ref_min_energy_conf}")      
        if plugin1_min_energy_conf != plugin2_min_energy_conf:
            print(f"  >> WARNING! Different minimum-energy configuration for {element}! {plugin1_min_energy_conf} (plugin1) vs {plugin2_min_energy_conf} (plugin 2)")

        unaries_data[element] = {'configurations': {}, 'reference_configuration': ref_min_energy_conf}
        for configuration in UNARIES_CONFIGURATIONS:
            unaries_data[element]['configurations'][configuration] = {}

        for plugin_energy_dict, plugin_name in [
            (plugin1_energy_dict, plugin_name_1),
            (plugin2_energy_dict, plugin_name_2),
        ]:
            if VERBOSE:
                print(f"  - {plugin_name:10s} [eV/atom]:", end="")
                print("    ", end="")
            for configuration in UNARIES_CONFIGURATIONS:
                en_diff = plugin_energy_dict[configuration] - plugin_energy_dict[ref_min_energy_conf]
                unaries_data[element]['configurations'][configuration][plugin_name] = en_diff
                if configuration != ref_min_energy_conf and VERBOSE:
                    print(f"({configuration:9s}: {en_diff:+4.3f})  ", end="")
            if VERBOSE:
                print() # Newline

    return unaries_data

def generate_json_data(plugin_name_1, plugin_name_2, unaries_set="unaries-verification-PBE-v1", oxides_set="oxides-verification-PBE-v1"):

    try:
        with open(f'../temp-results/results-{unaries_set}-{plugin_name_1}.json') as fhandle:
            plugin1_unaries = json.load(fhandle)
        with open(f'../temp-results/results-{unaries_set}-{plugin_name_2}.json') as fhandle:
            plugin2_unaries = json.load(fhandle)
        with open(f'../temp-results/results-{oxides_set}-{plugin_name_1}.json') as fhandle:
            plugin1_oxides = json.load(fhandle)
        with open(f'../temp-results/results-{oxides_set}-{plugin_name_2}.json') as fhandle:
            plugin2_oxides = json.load(fhandle)
    except OSError as exc:
        print(f"Error, probably you didn't generate the file yet! Original error: {exc}")
        sys.exit(1)

    all_data = {
        'formation_energies': {},
        'unaries_energy_difference': get_unaries_energy_difference(
            plugin_name_1, plugin_name_2, plugin1_unaries, plugin2_unaries)
        }
    for element in ALL_ELEMENTS:
        for configuration_triple in [
                #('X2O', 'XO', 'XO3'),
                #('X2O', 'XO2', 'XO3'),
                #('X2O', 'X2O3', 'XO3'),
                #('X2O', 'X2O5', 'XO3'),
                # 'X' is a placeholder for the lowest-energy 'X/...' phase, and 'O' for the
                # lowest-energy 'X/...' phase with X=oxygen
                ('X', 'X2O', 'O'),
                ('X', 'XO', 'O'),
                ('X', 'XO2', 'O'),
                ('X', 'X2O3', 'O'),
                ('X', 'X2O5', 'O'),
                ('X', 'XO3', 'O'),
            ]:

            plugin1_formation_energy = get_formation_energy(element, configuration_triple, plugin1_unaries, plugin1_oxides, all_data['unaries_energy_difference'])
            plugin2_formation_energy = get_formation_energy(element, configuration_triple, plugin2_unaries, plugin2_oxides, all_data['unaries_energy_difference'])

            all_data['formation_energies'][f"{element}-{configuration_triple[0]}|{configuration_triple[1]}|{configuration_triple[2]}"] = {
                plugin_name_1: plugin1_formation_energy,
                plugin_name_2: plugin2_formation_energy
            }

    # Plotting now
    OUT_FOLDER = 'formation-energies-output'
    os.makedirs(OUT_FOLDER, exist_ok=True)

    fname = f'{OUT_FOLDER}/formation-energies-{plugin_name_1}-VS-{plugin_name_2}.json'
    with open(fname, 'w') as fhandle:
        json.dump(all_data, fhandle, sort_keys=True, indent=2)
    print(f"File '{fname}' written.")

if __name__ == "__main__":
    try:
        plugin1_name = sys.argv[1]
        plugin2_name = sys.argv[2]
    except IndexError:
        print("Pass as two parameters the two plugins to compare")
        sys.exit(1)
    generate_json_data(plugin1_name, plugin2_name)
