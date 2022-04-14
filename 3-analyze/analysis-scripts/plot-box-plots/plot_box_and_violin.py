# %%
import json
import pathlib as pl
import numpy as np
import matplotlib.pyplot as plt
import tabulate

import quantities_for_comparison as qc

EXPECTED_SCRIPT_VERSION = '0.0.3'
DEFAULT_PREFACTOR = 100
DEFAULT_WB0 = 1.0 / 8.0
DEFAULT_WB01 = 1.0 / 64.0

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
set_name = 'oxides-verification-PBE-v1'
reference_plugin_name = 'wien2k_dk_0.06'
data_dir = pl.Path(f'../../outputs')
quantity_names = ['min_volume', 'bulk_modulus_ev_ang3', 'bulk_deriv', 'delta_per_formula_unit', 'epsilon']
# %%
data = {}
for filename in data_dir.glob(f'results-{set_name}*json'):
    plugin_name = filename.stem.replace(f'results-{set_name}-', '')
    with open(filename, 'r') as fp:
        plugin_data = json.load(fp)
        if plugin_data['script_version'] == EXPECTED_SCRIPT_VERSION:
            data[plugin_name] = plugin_data

pretty_plugin_name_map = {
    'quantum_espresso': 'Quantum ESPRESSO',
    'mk2-castep': 'CASTEP v2',
    'abinit-PseudoDojo-0.4-PBE-SR-standard-psp8': 'ABINIT NC v0.4',
    'cp2k_DZVP': 'CP2K DZVP',
    'fleur': 'Fleur',
    'castep': 'CASTEP',
    'gpaw': 'gpaw',
    'abinit-PseudoDojo-0.5b1-PBE-SR-standard-psp8': 'ABINIT NC v0.5b1',
    'vasp_recpot_e700': 'VASP recpot e700',
    'wien2k_dk_0.06': 'Wien2k',
    'vasp_recpot': 'VASP recpot',
    'abinit-JTH-1.1-PBE': 'ABINIT PAW v1.1'
}

quantity_name_map = {
    'bulk_modulus_ev_ang3': '$\Delta B_0$ [%]',
    'E0': '$E_0$ [eV]',
    'bulk_deriv': '$\Delta B_1$ [%]',
    'min_volume': '$\Delta V_0$ [%]',
    'delta_per_formula_unit': '$\Delta$',
    'epsilon': r'$\varepsilon$'
}

text_quantity_name_map = {
    'bulk_modulus_ev_ang3': 'dB0 [%]',
    'bulk_deriv': 'dB1 [%]',
    'min_volume': 'dV0 [%]',
    'delta_per_formula_unit': 'Delta',
    'epsilon': 'Epsilon'
}

quantity_for_comparison_map = {
    "delta_per_formula_unit": qc.delta,
    "bulk_modulus_ev_ang3": qc.B0_rel_diff,
    "min_volume": qc.V0_rel_diff,
    "bulk_deriv": qc.B1_rel_diff,
    "rel_errors_vec_length": qc.rel_errors_vec_length,
    "epsilon": qc.epsilon
}
# %%
# Sort in alphabetical order
plugin_names = sorted(list(data.keys()))[::-1]
# Remove the reference plugin from the list
plugin_names = [plugin_name for plugin_name in plugin_names if plugin_name != reference_plugin_name]

ref_data = data[reference_plugin_name]
ref_BM_fit_data = ref_data['BM_fit_data']
# List the reference systems that have BM fit data
ref_systems = set(
    el_conf for el_conf in ref_BM_fit_data.keys()
    if ref_BM_fit_data[el_conf] is not None
)

# Set up a place to store all the plotting data
quantity_data = {plugin_name: {} for plugin_name in plugin_names}
for quantity_name in quantity_names:
    for plugin_name in plugin_names:
        plugin_values = []

        plugin_data = data[plugin_name]
        plugin_BM_fit_data = data[plugin_name]['BM_fit_data']

        plugin_systems = set(
            el_conf for el_conf in plugin_BM_fit_data.keys()
            if plugin_BM_fit_data[el_conf] is not None
        )
        # Take the systems that are both in the reference and plugin sets
        plot_systems = plugin_systems.intersection(ref_systems)

        for element_and_configuration in plot_systems:
            element, configuration = element_and_configuration.split('-')
            
            ref_n_atoms = ref_data['num_atoms_in_sim_cell'][f'{element}-{configuration}']
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
            quantity_data[plugin_name][quantity_name] = plugin_values

# %%
# Set up the plot axes for each quantity
n_quantities = len(quantity_names)
fig, axes = plt.subplots(1, n_quantities, dpi=300, figsize=(4 * n_quantities, 4), sharey=True)
axes = axes.flatten()

for quantity_name, ax in zip(quantity_names, axes):
    quantity_values = [quantity_data[plugin_name][quantity_name] for plugin_name in plugin_names]
    
    violin_parts = ax.violinplot(
        quantity_values,
        widths=0.75,
        points=1000,
        showextrema=False, showmedians=False,  # These are shown by the box plots
        vert=False
    )
    for part in violin_parts['bodies']:
        part.set_facecolor('tab:grey')
        part.set_alpha(0.50)
    
    ax.boxplot(
        quantity_values,
        flierprops=FLIERPROPS,
        boxprops=BOXPROPS,
        whiskerprops=WHISKERPROPS,
        medianprops=MEDIANPROPS,
        capprops=CAPPROPS,
        vert=False,
        labels=[pretty_plugin_name_map[plugin_name] for plugin_name in plugin_names]
    )
    ax.set_xlabel(quantity_name_map[quantity_name], fontsize=14)

fig.suptitle(f'{pretty_plugin_name_map[reference_plugin_name]} Reference')
fig.tight_layout()

fig.savefig('boxplots.pdf')

# %%
def make_quantity_table(function, quantity_data):
    table = []
    for plugin_name in plugin_names:
        row = [pretty_plugin_name_map[plugin_name]]
        for quantity_name in quantity_names:
            val = function(quantity_data[plugin_name][quantity_name])
            row.append(val)
        table.append(row)
    table = table[::-1]

    return tabulate.tabulate(
        table, headers=['Plugin'] + [text_quantity_name_map[quantity_name] for quantity_name in quantity_names]
    )

# %%
table = make_quantity_table(lambda x: np.median(x), quantity_data)
print('MEDIAN\n' + table)
# %%
table = make_quantity_table(lambda x: np.mean(x), quantity_data)
print('MEAN\n' + table)
# %%
table = make_quantity_table(lambda x: np.quantile(x, 0.75) - np.quantile(x, 0.25), quantity_data)
print('IQR\n' + table)

# %%
