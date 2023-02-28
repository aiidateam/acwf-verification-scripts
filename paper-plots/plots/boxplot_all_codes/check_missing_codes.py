#!/usr/bin/env python
import json
import pathlib as pl
import os
import ase.data



if __name__ == "__main__":
    DATA_FOLDER = "../../code-data"
    with open(os.path.join(DATA_FOLDER, "labels.json")) as fhandle:
        labels_data = json.load(fhandle)
    code_labels = labels_data['methods-main'].keys()

    SET_1 = ('unaries',  ['X/SC', 'X/FCC', 'X/BCC', 'X/Diamond'])
    SET_2 = ('oxides', ['XO', 'XO2', 'X2O', 'X2O3', 'X2O5', 'XO3'])
    #chemical_numbers = list(range(1, 96+1)) # ALL
    chemical_numbers = list(range(1, 56+1)) + list(range(72, 83+1)) # Minimal set: H to Bi except Lanthanides (La to Lu)

    SKIP_BIGDFT = False

    plugin_cache = {}
    def get_plugin_fit_data(set_name, plugin_name):
        global plugin_cache, labels_data

        key = f'{set_name}-{plugin_name}'

        try:
            return plugin_cache[key]
        except KeyError:
            with open(os.path.join(DATA_FOLDER, labels_data['methods-main'][plugin_name][set_name])) as fhandle:
                plugin_cache[key] = json.load(fhandle)['BM_fit_data']
            return plugin_cache[key]


    for Z in chemical_numbers:
        for set_name, variants in (SET_1, SET_2):
            symbol = ase.data.chemical_symbols[Z]    
            for variant in variants:
                configuration = f'{symbol}-{variant}'
                for code_label in code_labels:
                    # Not ideal as it is opening a lot of times, but this gives the order I want
                    plugin_fit_data = get_plugin_fit_data(set_name, code_label)

                    if SKIP_BIGDFT and 'bigdft' in code_label.lower():
                        continue
                    if plugin_fit_data.get(configuration) is None:
                        print(f">> {Z}: {configuration} ({code_label})")

