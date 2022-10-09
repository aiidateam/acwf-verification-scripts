# %%
import json
import pathlib as pl
import numpy as np
import matplotlib.pyplot as plt
import tabulate
import sys

REF_PLUGIN = "ae"

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
labels = [
        'wien2k',
        'vasp',
        'quantum-espresso',
        'siesta',
        'gpaw',
        'fleur',
        'cp2k',
        'castep',
        'bigdft',
        'abinit (PseudoDojo0.4)',
        'abinit (PseudoDojo0.5b1)'
        ]

def get_list(set_name, file_name):
    """
    """
    
    for SET_NAME in set_name:
        try:
            data_dir = pl.Path(f'../data')
            with open(f'{data_dir}/results-{SET_NAME}-{REF_PLUGIN}.json') as fhandle:
                ref_plugin_data = json.load(fhandle)
        except OSError:
            #print(f"No data found for {REF_PLUGIN} (set '{SET_NAME}'), {REF_PLUGIN} is the reference and must be present")
            sys.exit(1)
        ref_BM_fit_data = ref_plugin_data['BM_fit_data']
        ref_systems = set(
            el_conf for el_conf in ref_BM_fit_data.keys() if ref_BM_fit_data[el_conf] is not None
        )
        common_data = ref_systems
        print(f'{SET_NAME}  {REF_PLUGIN}', len(ref_systems))

        for plugin_name in plugin_names:
            plugin_values = []
            try:
                data_dir = pl.Path(f'../data')
                with open(f'{data_dir}/results-{SET_NAME}-{plugin_name}.json') as fhandle:
                    plugin_data = json.load(fhandle)
            except OSError:
                print(f"No data found for {plugin_name} (set '{SET_NAME}')")
                sys.exit(1)

            plugin_BM_fit_data = plugin_data['BM_fit_data']

            plugin_systems = set(
                el_conf for el_conf in plugin_BM_fit_data.keys()
                if plugin_BM_fit_data[el_conf] is not None
            )
            # Take the systems that are both in the reference and plugin sets
            common_data = common_data.intersection(plugin_systems)
            print(f'{SET_NAME}  {plugin_name}', len(plugin_systems))

        print(common_data)

        


if __name__ == "__main__":
    SET_NAME_1 = 'unaries-verification-PBE-v1'
    SET_NAME_2 = 'oxides-verification-PBE-v1'

    get_list([SET_NAME_1, SET_NAME_2], 'box_plot_all')

