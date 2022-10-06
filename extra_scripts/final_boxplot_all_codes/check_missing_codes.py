#!/usr/bin/env python
import json
import pathlib as pl

plugin_names = [
        'wien2k',
        'vasp',
        'quantum_espresso',
        'siesta',
        'gpaw',
        'fleur',
        'cp2k', 
        'castep', 
        'bigdft', 
        'abinit-PseudoDojo0.4', 
        'abinit-PseudoDojo0.5b1'
        ] # need them in reverse order for up to bottom in plot


if __name__ == "__main__":
    import ase.data
    SET_1 = ('unaries-verification-PBE-v1',  ['X/SC', 'X/FCC', 'X/BCC', 'X/Diamond'])
    SET_2 = ('oxides-verification-PBE-v1', ['XO', 'XO2', 'X2O', 'X2O3', 'X2O5', 'XO3'])
    #chemical_numbers = list(range(1, 96+1)) # ALL
    chemical_numbers = list(range(1, 56+1)) + list(range(72, 83+1)) # Minimal set: H to Pb except Lanthanides

    SKIP_BIGDFT = True

    plugin_cache = {}
    def get_plugin_fit_data(set_name, plugin_name):
        key = f'{set_name}-{plugin_name}'
        global plugin_cache

        try:
            return plugin_cache[key]
        except KeyError:
            with open(f'{data_dir}/results-{set_name}-{plugin_name}.json') as fhandle:
                plugin_cache[key] = json.load(fhandle)['BM_fit_data']
            return plugin_cache[key]


    for Z in chemical_numbers:
        for set_name, variants in (SET_1, SET_2):
            symbol = ase.data.chemical_symbols[Z]    
            for variant in variants:
                configuration = f'{symbol}-{variant}'
                for plugin_name in plugin_names:
                    # Not ideal as it is opening a lot of times, but this gives the order I want
                    data_dir = pl.Path(f'data')
                    plugin_fit_data = get_plugin_fit_data(set_name, plugin_name)

                    if SKIP_BIGDFT and plugin_name == 'bigdft':
                        continue
                    if plugin_fit_data.get(configuration) is None:
                        print(f">> {Z}: {configuration} ({plugin_name})")

