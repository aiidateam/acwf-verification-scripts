# %%
import json
import pathlib as pl
from more_itertools import only
import numpy as np
import matplotlib.pyplot as plt
import tabulate
import sys

import quantities_for_comparison as qc

EXPECTED_SCRIPT_VERSION = ['0.0.3','0.0.4']
DEFAULT_PREFACTOR = 100
DEFAULT_WB0 = 1.0 / 8.0
DEFAULT_WB01 = 1.0 / 64.0
REF_PLUGIN = "ae"

CAPPROPS = {
    'linewidth': 1,
    'color': 'tab:green'
}
BOXPROPS = {
    'linewidth': 1,
    'color': 'tab:green'
}
WHISKERPROPS = {
    'linewidth': 1,
    'color': 'tab:green'
}
MEDIANPROPS = {
    'linewidth': 1,
    'color': 'tab:blue'
}
FLIERPROPS = {
    'marker': '.',
    'linestyle': '',
    'markerfacecolor': 'tab:grey',
    'markeredgewidth': 0,
    'markersize': 3
}
# %%

quantity_for_comparison_map = {
    "delta_per_formula_unit": qc.delta,
    "% difference in B0": qc.B0_rel_diff,
    "% difference in V0": qc.V0_rel_diff,
    "% difference in B1": qc.B1_rel_diff,
    "rel_errors_vec_length": qc.rel_errors_vec_length,
    "epsilon": qc.epsilon
}

xlims = {
    "delta_per_formula_unit": [0,3],
    "% difference in B0": [-20,20],
    "% difference in V0": [-8.0,8.0],
    "% difference in B1": [-50,50],
    "rel_errors_vec_length": [0,15],
    "epsilon": [0,15]
}

quantity_names = ["% difference in V0","% difference in B0","% difference in B1"]

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

def generate_box_plt(set_name, file_name, only_must_have_elements=None, skip_codes=None):
    if skip_codes is not None:
        for skip_code in skip_codes:
            try:
                position = plugin_names.index(skip_code)
                plugin_names.pop(position)
                labels.pop(position)
                print(f">>>>>> SKIPPING CODE '{skip_code}'")
            except IndexError:
                raise ValueError(f"Code '{skip_code}' asked to be skipped but does not exist!")
    all_data = {}
    for quantity_name in quantity_names:
        out_data = {}
        print(f"#### {quantity_name} ####")
        for plugin_name in plugin_names:
            plugin_values = []
            plugin_big = 0
            plugin_small = 0
            out_data[plugin_name] = {}
            for SET_NAME in set_name:
                try:
                    data_dir = pl.Path(f'data')
                    with open(f'{data_dir}/results-{SET_NAME}-{REF_PLUGIN}.json') as fhandle:
                        ref_plugin_data = json.load(fhandle)
                except OSError:
                    print(f"No data found for {REF_PLUGIN} (set '{SET_NAME}'), {REF_PLUGIN} is the reference and must be present")
                    sys.exit(1)


                try:
                    data_dir = pl.Path(f'data')
                    with open(f'{data_dir}/results-{SET_NAME}-{plugin_name}.json') as fhandle:
                        plugin_data = json.load(fhandle)
                except OSError:
                    print(f"No data found for {REF_PLUGIN} (set '{SET_NAME}'), {REF_PLUGIN} is the reference and must be present")
                    sys.exit(1)

                ref_BM_fit_data = ref_plugin_data['BM_fit_data']
                # List the reference systems that have BM fit data
                ref_systems = set(
                    el_conf for el_conf in ref_BM_fit_data.keys()
                    if ref_BM_fit_data[el_conf] is not None
                )

                plugin_BM_fit_data = plugin_data['BM_fit_data']

                plugin_systems = set(
                    el_conf for el_conf in plugin_BM_fit_data.keys()
                    if plugin_BM_fit_data[el_conf] is not None
                )
                # Take the systems that are both in the reference and plugin sets
                plot_systems = plugin_systems.intersection(ref_systems)
                set_type = SET_NAME.partition('-')[0] # oxides or unaries
                print(f"All systems: {len(plot_systems)} (set: {set_type}; full: {len(plugin_systems)}) for '{plugin_name}'")
                if only_must_have_elements is not None:
                    if set_type == 'oxides':
                        variants = ['XO', 'XO2', 'X2O', 'X2O3', 'X2O5', 'XO3']
                    elif set_type == 'unaries':
                        variants = ['X/SC', 'X/FCC', 'X/BCC', 'X/Diamond']
                    else:
                        raise ValueError("Unrecognized set type!")
                    missing = []
                    new_plot_systems = set()
                    for element in only_must_have_elements:
                        for variant in variants:
                            expected_key = f"{element}-{variant}"
                            if expected_key not in plot_systems:
                                missing.append(expected_key)
                            else:
                                new_plot_systems.add(expected_key)
                    if missing:
                        raise ValueError(f"{plugin_name} ({set_type}) misses the following keys: {missing}")
                    plot_systems = new_plot_systems
                print(f"   -> Plotting: {len(plot_systems)}")


                for element_and_configuration in plot_systems:
                    element, configuration = element_and_configuration.split('-')
                
                    ref_n_atoms = ref_plugin_data['num_atoms_in_sim_cell'][f'{element}-{configuration}']
                    ref_scaling_factor = qc.get_volume_scaling_to_formula_unit(
                        ref_n_atoms, element, configuration
                    )
                    ref_V0 = ref_BM_fit_data[f'{element}-{configuration}']['min_volume'] / ref_scaling_factor
                    ref_B0 = ref_BM_fit_data[f'{element}-{configuration}']['bulk_modulus_ev_ang3']
                    ref_B01 = ref_BM_fit_data[f'{element}-{configuration}']['bulk_deriv']

                    plugin_n_atoms = plugin_data['num_atoms_in_sim_cell'][f'{element}-{configuration}']
                    plugin_scaling_factor = qc.get_volume_scaling_to_formula_unit(
                        plugin_n_atoms, element, configuration
                    )
                    plugin_V0 = plugin_BM_fit_data[f'{element}-{configuration}']['min_volume'] / plugin_scaling_factor
                    plugin_B0 = plugin_BM_fit_data[f'{element}-{configuration}']['bulk_modulus_ev_ang3']
                    plugin_B01 = plugin_BM_fit_data[f'{element}-{configuration}']['bulk_deriv']

                    quantity_value = quantity_for_comparison_map[quantity_name](
                        ref_V0, ref_B0, ref_B01,
                        plugin_V0, plugin_B0, plugin_B01,
                        DEFAULT_PREFACTOR, DEFAULT_WB0, DEFAULT_WB01
                    )

                    plugin_values.append(quantity_value)
                    if quantity_value < xlims[quantity_name][0]:
                        plugin_big = plugin_big+1  
                    if quantity_value > xlims[quantity_name][1]:
                        plugin_small = plugin_small+1

            out_data[plugin_name]['values'] = plugin_values
            out_data[plugin_name]['big'] = plugin_big
            out_data[plugin_name]['small'] = plugin_small
        all_data[quantity_name] = out_data
    # %%
    # Set up the plot axes for each quantity
    for yy in quantity_names:
        for plugin_name in plugin_names:
            if all_data[yy][plugin_name]['small'] != 0:
                print(f'small {yy} {plugin_name}', all_data[yy][plugin_name]['small'])
    for yy in quantity_names:
        for plugin_name in plugin_names:
            if all_data[yy][plugin_name]['big'] != 0:
                print(f'big {yy} {plugin_name}', all_data[yy][plugin_name]['big'])
    
    n_quantities = len(quantity_names)
    fig, axes = plt.subplots(1, n_quantities, dpi=300, figsize=(4 * n_quantities, 4), sharey=True)
    axes = axes.flatten()

    for quantity_name, ax in zip(quantity_names, axes):
        quantity_values = [all_data[quantity_name][plugin_name]['values'] for plugin_name in plugin_names]
    
        #violin_parts = ax.violinplot(
        #    quantity_values,
        #    widths=0.75,
        #    points=1000,
        #    showextrema=False, showmedians=False,  # These are shown by the box plots
        #    vert=False
        #)
        #for part in violin_parts['bodies']:
        #    part.set_facecolor('tab:grey')
        #    part.set_alpha(0.50)
    
        ax.boxplot(
            quantity_values,
            flierprops=FLIERPROPS,
            boxprops=BOXPROPS,
            whiskerprops=WHISKERPROPS,
            medianprops=MEDIANPROPS,
            capprops=CAPPROPS,
            vert=False,
            labels=labels
        )
        ax.set_xlabel(quantity_name, fontsize=14)
        ax.set_xlim(xlims[quantity_name][0],xlims[quantity_name][1])

    fig.suptitle(f'Reference: all electron average.')
    fig.tight_layout()

    suffix = "_SUBSET" if only_must_have_elements is not None else ""
    fig.savefig(f'{file_name}{suffix}.pdf')


if __name__ == "__main__":
    import ase.data
    SET_NAME_1 = 'unaries-verification-PBE-v1'
    SET_NAME_2 = 'oxides-verification-PBE-v1'
    ONLY_SUBSET = False

    try:
        mode = sys.argv[1]
    except IndexError:
        print("Pass as first parameter 'separate' (separate analysis for unaries and oxides) or 'all'.")
        sys.exit(1)

    only_must_have_elements = None
    skip_codes = None


    if ONLY_SUBSET:
        # ONLY PLOT A SUBSET, AND SKIP BIGDFT
        
        # except f
        chemical_numbers = list(range(1, 56+1)) + list(range(72, 86+1))
        # VASP
        chemical_numbers.remove(1) #H
        chemical_numbers.remove(2) #He
        chemical_numbers.remove(4) #Be
        chemical_numbers.remove(45) #Rh

        # QE
        chemical_numbers.remove(85) #At
        chemical_numbers.remove(18) #Ar
        chemical_numbers.remove(53) #I

        # Siesta
        chemical_numbers.remove(80) #Hg
        chemical_numbers.remove(37) #Rb

        # GPAW
        chemical_numbers.remove(83) #Bi
        chemical_numbers.remove(43) #Tc
        chemical_numbers.remove(86) #Bi
        chemical_numbers.remove(84) #Po

        # CP2K
        chemical_numbers.remove(11) #Na

        only_must_have_elements = [ase.data.chemical_symbols[i] for i in
            chemical_numbers]
        skip_codes = ['bigdft']

    if mode == 'all':
        generate_box_plt([SET_NAME_1, SET_NAME_2], 'box_plot_all',
            only_must_have_elements=only_must_have_elements, skip_codes=skip_codes)
    elif mode == 'separate':
        generate_box_plt([SET_NAME_1], 'box_plot_unaries',
            only_must_have_elements=only_must_have_elements, skip_codes=skip_codes)
        generate_box_plt([SET_NAME_2], 'box_plot_oxides',    
            only_must_have_elements=only_must_have_elements, skip_codes=skip_codes)
    else:
        print("Pass as first parameter 'separate' (separate analysis for unaries and oxides) or 'all'.")
        sys.exit(1)

