#!/usr/bin/env python
import json
import matplotlib.pyplot as plt
import os
import sys
import copy
import acwf_paper_plots.quantities_for_comparison as qc

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})

EXPECTED_SCRIPT_VERSION = ['0.0.3','0.0.4']
DEFAULT_PREFACTOR = 100 # to convert from relative to % errors
DEFAULT_WB0 = 0.
DEFAULT_WB01 = 0.

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
}

xlims = {
    "delta_per_formula_unit": [0,3],
    "% difference in B0": [-20,20],
    "% difference in V0": [-7,7],
    "% difference in B1": [-50,50],
}

quantity_names = ["% difference in V0","% difference in B0","% difference in B1"]

QUANTITY_FANCY_NAMES = {
    'B0': "$B_0$",
    'V0': "$V_0$",
    'B1': "$B_1$"
}

DATA_FOLDER = "../../../code-data"
with open(os.path.join(DATA_FOLDER, "labels.json")) as fhandle:
    labels_data = json.load(fhandle)
code_labels = list(labels_data['methods-main'].keys())[::-1] # invert order because they are plot bottom to top, so we keep them alphabetical

ALL_ELECTRON_CODES_SHORT = ["FLEUR", "WIEN2k"][::-1] # Revert order as they are printed from top to bottom
ALL_ELECTRON_CODES = [labels_data['all-electron-keys'][short_label] for short_label in ALL_ELECTRON_CODES_SHORT]


def generate_box_plt(set_names, file_name, material_set_label, file_suffix, only_must_have_elements=None, keep_only_codes=None):
    """
    Generate the box plot
    """
    # Double check that there are no mistakes
    if keep_only_codes is not None:
        for code in keep_only_codes:
            if code not in code_labels:
                raise ValueError(f"Asking to keep code '{code}' but it does not exist a code with such a label")

    for code in ALL_ELECTRON_CODES:
        if code not in code_labels:
            raise ValueError(f"All electron code '{code}' expected, but no data exists for a code with such a label")

    # Get the filtered list
    used_code_labels = []
    for code_label in code_labels:
        if keep_only_codes is not None and code_label not in keep_only_codes:
            print(f">> Skipping code {code_label}")
            continue
        # Skip the AE, will add at the end (so they are on the top)
        if code_label in ALL_ELECTRON_CODES:
            continue
        used_code_labels.append(code_label)
    for code_label in ALL_ELECTRON_CODES:
        used_code_labels.append(code_label)

    all_data = {}
    print()
    print('#####################################################################')
    print('#      Statistics on the number of elements for each code           #')
    print('#####################################################################')
    for quantity_name in quantity_names:
        out_data = {}
        for code_label in used_code_labels:
            plugin_values = []
            plugin_big = 0
            plugin_small = 0
            out_data[code_label] = {}
            for set_name in set_names:
                reference_data_files = labels_data['references']['all-electron average']
                with open(os.path.join(DATA_FOLDER, reference_data_files[set_name])) as fhandle:
                    ref_plugin_data = json.load(fhandle)
                with open(os.path.join(DATA_FOLDER, labels_data['methods-main'][code_label][set_name])) as fhandle:
                    plugin_data = json.load(fhandle)

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
                               
                if set_name == 'oxides':
                    variants = ['XO', 'XO2', 'X2O', 'X2O3', 'X2O5', 'XO3']
                elif set_name == 'unaries':
                    variants = ['X/SC', 'X/FCC', 'X/BCC', 'X/Diamond']
                else:
                    raise ValueError("Unrecognized set name!")
                missing = []
                new_plot_systems = set()
                for element in only_must_have_elements:
                    for variant in variants:
                        expected_key = f"{element}-{variant}"
                        if expected_key not in plot_systems:
                            missing.append(expected_key)
                        else:
                            new_plot_systems.add(expected_key)
                if quantity_name == '% difference in V0':
                    if missing:
                        print(f"{code_label} ({set_name}) misses the following keys: {missing}")
                        print(f"   -> Plotting: {len(new_plot_systems)}")
                    else:
                        print(f"{code_label} ({set_name}) is complete")
                        print(f"   -> Plotting: {len(new_plot_systems)}")
                plot_systems = new_plot_systems

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

            out_data[code_label]['values'] = plugin_values
            out_data[code_label]['big'] = plugin_big
            out_data[code_label]['small'] = plugin_small
        all_data[quantity_name] = out_data
    # %%
    # Set up the plot axes for each quantity
    print()
    print('#####################################################################')
    print('#   Description of the outliers (out of picture) for each code      #')
    print('#####################################################################')
    for yy in quantity_names:
        for code_label in used_code_labels:
            if all_data[yy][code_label]['small'] != 0:
                print(f'small {yy} {code_label}', all_data[yy][code_label]['small'])
    for yy in quantity_names:
        for code_label in used_code_labels:
            if all_data[yy][code_label]['big'] != 0:
                print(f'big {yy} {code_label}', all_data[yy][code_label]['big'])
    
    n_quantities = len(quantity_names)
    # + 2 is to keep space for header and footer
    fig_height = max((len(used_code_labels) + 2) * 0.45, 3) * 0.7 # Set min size of 3; rescaling factor to make it less wide
    fig, axes = plt.subplots(1, n_quantities, dpi=300, figsize=(4 * n_quantities, fig_height), sharey=True)
    axes = axes.flatten()

    for quantity_name, ax in zip(quantity_names, axes):
        quantity_values = [all_data[quantity_name][plugin_name]['values'] for plugin_name in used_code_labels]
    
        fancy_labels = []
        for label in used_code_labels:
            code_label, sep, rest = label.partition('@')
            # Use bold code name, and replace pipe as with OT1 font enc it is replaced by a dash
            fancy_labels.append(rf'\textbf{{{code_label}}}{sep}{rest}'.replace('|', '$|$'))

        ax.boxplot(
            quantity_values,
            flierprops=FLIERPROPS,
            boxprops=BOXPROPS,
            whiskerprops=WHISKERPROPS,
            medianprops=MEDIANPROPS,
            capprops=CAPPROPS,
            vert=False,
            labels=fancy_labels
        )
        ax.set_xlabel(rf"{QUANTITY_FANCY_NAMES[quantity_name[-2:]]} difference [\%]",fontsize=14)
        ax.tick_params(axis='x', which='major', labelsize=11)
        ax.tick_params(axis='y', which='major', labelsize=9)
        ax.set_xlim(xlims[quantity_name][0],xlims[quantity_name][1])

        # Put a horizontal dashed line to separate all-electron codes from the rest
        ax.axhline(len(used_code_labels) - len(ALL_ELECTRON_CODES) + 0.5, linestyle='--', color='gray', linewidth=1)

    fig.suptitle(f"Materials set: {material_set_label}",fontsize=14, y=0.98)
    #fig.tight_layout()
    fig.subplots_adjust(left=0.25, right=0.99, top=1., bottom=(1.8/(len(used_code_labels) + 2)), wspace=0.05)
    def make_space_above(axes, topmargin=1):
        """ increase figure size to make topmargin (in inches) space for 
            titles, without changing the axes sizes"""
        fig = axes.flatten()[0].figure
        s = fig.subplotpars
        w, h = fig.get_size_inches()

        figh = h - (1-s.top)*h  + topmargin
        fig.subplots_adjust(bottom=s.bottom*h/figh, top=1-topmargin/figh)
        fig.set_figheight(figh)
    make_space_above(axes, topmargin=0.35)
    fig.savefig(f'{file_name}{file_suffix}.pdf')


if __name__ == "__main__":
    import ase.data
    SET_NAME_1 = 'unaries'
    SET_NAME_2 = 'oxides'

    try:
        mode = sys.argv[1]
    except IndexError:
        print("Pass as first parameter 'separate' (separate analysis for unaries and oxides) or 'together'.")
        sys.exit(1)

    try:
        elements = sys.argv[2]
    except IndexError:
        print(
            "Pass as second parameter 'all' (atomic number 1-96), 'up-to-Bi-no-lanthanides' (1-56,71-83), "
            "'delta-set' (1-56,71-84+86), 'no-actinides' (1-86), 'only-actinides' (84-96), 'only-lanthanides'(57-71)."
        )
        sys.exit(1)

    only_must_have_elements = None
    if elements == 'all':
        chemical_numbers = list(range(1, 96+1))
        material_set_label = "Z=1-96"
    elif elements == 'delta-set':
        chemical_numbers = list(range(1, 56+1)) +  list(range(71, 84+1)) + [86]
        material_set_label = "Delta set (Science 2016)"
    elif elements == 'up-to-Bi-no-lanthanides':
        chemical_numbers = list(range(1, 56+1)) +  list(range(72, 83+1))
        material_set_label = "Z=1-56,72-83"
    elif elements == 'no-actinides':
        chemical_numbers = list(range(1, 88+1))
        material_set_label = "Z=1-88"
    elif elements == 'only-actinides':
        chemical_numbers = list(range(84, 96+1))
        material_set_label = "Z=84-96"
    elif elements == 'only-lanthanides':
        chemical_numbers = list(range(57, 71+1))
        material_set_label = "Z=57-71"
    else:
        print(
            "Pass as second parameter 'all' (atomic number 1-96), 'up-to-Bi-no-lanthanides' (1-56,71-83), "
            "'delta-set' (1-56,71-84+86), 'no-actinides' (1-86), 'only-actinides' (84-96), 'only-lanthanides'(57-71)."
        )
        sys.exit(1)
    # VASP  chemical_numbers.remove(1) #H
    #       chemical_numbers.remove(2) #He
    #       chemical_numbers.remove(4) #Be
    #       chemical_numbers.remove(45) #Rh
    # QE  chemical_numbers.remove(85) #At
    #     chemical_numbers.remove(18) #Ar
    #     chemical_numbers.remove(53) #I
    # Siesta     chemical_numbers.remove(80) #Hg
    #            chemical_numbers.remove(37) #Rb
    # GPAW     chemical_numbers.remove(83) #Bi
    #          chemical_numbers.remove(43) #Tc
    #          chemical_numbers.remove(86) #Bi
    #          chemical_numbers.remove(84) #Po
    # CP2K    chemical_numbers.remove(11) #Na
    only_must_have_elements = [ase.data.chemical_symbols[i] for i in chemical_numbers]
    
    keep_only_codes = sys.argv[3:]
    if not keep_only_codes:
        keep_only_codes = None

    if mode == 'together':
        generate_box_plt([SET_NAME_1, SET_NAME_2], 'box_plot_',
            only_must_have_elements=only_must_have_elements, material_set_label=material_set_label, file_suffix=elements, keep_only_codes=keep_only_codes)
    elif mode == 'separate':
        generate_box_plt([SET_NAME_1], 'box_plot_unaries_',
            only_must_have_elements=only_must_have_elements, material_set_label=material_set_label, file_suffix=elements, keep_only_codes=keep_only_codes)
        generate_box_plt([SET_NAME_2], 'box_plot_oxides_',    
            only_must_have_elements=only_must_have_elements, material_set_label=material_set_label, file_suffix=elements, keep_only_codes=keep_only_codes)
    else:
        print("Pass as first parameter 'separate' (separate analysis for unaries and oxides) or 'together'.")
        sys.exit(1)

