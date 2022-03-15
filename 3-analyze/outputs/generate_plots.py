#!/usr/bin/env python
import json
import os
import sys

import numpy as np
import pylab as pl
import tqdm

from quantities_for_comparison import birch_murnaghan, get_volume_scaling_to_formula_unit

def get_plugin_name():
    file_name = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        os.pardir, os.pardir, 'plugin_name.txt'
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

EXPECTED_SCRIPT_VERSION = '0.0.3'
RESIDUALS_THRESHOLD = 1.e-3

def get_conf_nice(configuration_string):
    """Convert the configuration string to a nicely typeset string in LaTeX."""
    ret_pieces = []
    for char in configuration_string:
        if char in "0123456789":
            ret_pieces.append(f"$_{char}$")
        else:
            ret_pieces.append(char)
    return "".join(ret_pieces)


if __name__ == "__main__":
    try:
        SET_NAME = sys.argv[1]
    except IndexError:
        print("Pass as first parameter the set name, e.g. oxides-verification-PBE-v1 or unaries-verification-PBE-v1")
        sys.exit(1)

    try:
        compare_with = sys.argv[2]
    except IndexError:
        compare_with = None

    try:
        with open(f'results-{SET_NAME}-{PLUGIN_NAME}.json') as fhandle:
            reference_plugin_data = json.load(fhandle)
    except OSError:
        print(f"No data found for your plugin '{PLUGIN_NAME}' (set '{SET_NAME}'). Did you run `./get_results.py` first?")
        sys.exit(1)
    
    if not reference_plugin_data['script_version'] == EXPECTED_SCRIPT_VERSION:
        raise ValueError(
            f"This script only works with data generated at version {EXPECTED_SCRIPT_VERSION}. "
            "Please re-run ./get_results.py to update the data format!")

    if compare_with is None:
        print(f"Plotting data for plugin '{PLUGIN_NAME}' only (set '{SET_NAME}').")
        print("If you want to compare, pass another parameter with the plugin to compare with.")
        compare_plugin_data = None
    else:
        print(f"Plotting data for plugin '{PLUGIN_NAME}' (set '{SET_NAME}') compared with '{compare_with}'.")
        try:
             with open(f'results-{SET_NAME}-{compare_with}.json') as fhandle:
                compare_plugin_data = json.load(fhandle)
        except OSError:
            print(f"No data found for the reference plugin '{compare_with}': you need the file results-{SET_NAME}-{compare_with}.json.")
            sys.exit(1)
        if not compare_plugin_data['script_version'] == EXPECTED_SCRIPT_VERSION:
            raise ValueError(
                f"This script only works with data generated at version {EXPECTED_SCRIPT_VERSION}. "
                "Please ask the other plugin you want to compare with to re-run ./get_results.py "
                "to update the data format!"
            )


    if compare_with is None:
        PLOT_FOLDER = f'plots-{SET_NAME}-{PLUGIN_NAME}'
    else:
        PLOT_FOLDER = f'plots-{SET_NAME}-{PLUGIN_NAME}-vs-{compare_with}'
    os.makedirs(PLOT_FOLDER, exist_ok=True)

    all_systems = set(reference_plugin_data['BM_fit_data'].keys())
    if compare_with:
        all_systems.update(compare_plugin_data['BM_fit_data'].keys())

    progress_bar = tqdm.tqdm(sorted(all_systems))
    for element_and_configuration in progress_bar:
        progress_bar.set_description(f"{element_and_configuration:12s}")
        progress_bar.refresh()

        element, configuration = element_and_configuration.split('-')
        try:
            eos_data = reference_plugin_data['eos_data'][f'{element}-{configuration}']
        except KeyError:
            # If this system does not exist in the reference data, skip it
            continue
        if eos_data is None:
            # If there is no data, I skip this material
            continue
        scaling_ref_plugin = get_volume_scaling_to_formula_unit(
            reference_plugin_data['num_atoms_in_sim_cell'][f'{element}-{configuration}'],
            element, configuration
        )

        # Get the x axis for the plot
        volumes, energies = (np.array(eos_data).T / scaling_ref_plugin).tolist()
        dense_volumes = np.linspace(
            min(volumes),
            max(volumes),
            100
        )

        # Get the data for the reference plugin
        try:
            ref_BM_fit_data = reference_plugin_data['BM_fit_data'][f'{element}-{configuration}']
            if ref_BM_fit_data is None:
                # No fitting data: data was there but was not fitted.
                # Raise this exception that is catched one line below, so
                # there is nothing plotted for the fit but just the data.
                raise KeyError
        except KeyError:
            # Set to None if fit data is missing (if we are here, the EOS points
            # are there, so it means that the fit failed). I will still plot the
            # points
            reference_eos_fit_energy = None
            residuals = None
        else:
            reference_eos_fit_energy = birch_murnaghan(
                V=dense_volumes,
                E0=ref_BM_fit_data['E0'] / scaling_ref_plugin,
                V0=ref_BM_fit_data['min_volume'] / scaling_ref_plugin,
                B0=ref_BM_fit_data['bulk_modulus_ev_ang3'],
                B01=ref_BM_fit_data['bulk_deriv']
            )
            residuals = ref_BM_fit_data['residuals']

            # Get the data for the compare_with plugin, if specified (and if the EOS worked for the 
            # reference plugin, otherwise we don't know which E0 to use)
            if compare_with is not None:
                try:
                    compare_BM_fit_data = compare_plugin_data['BM_fit_data'][f'{element}-{configuration}']
                    if compare_BM_fit_data is None:
                        # No fitting data in the plugin to compare with.
                        # Raise this exception that is catched one line below, so
                        # it will set `compare_eos_fit_energy` to None.
                        raise KeyError                    
                except KeyError:
                    # Set to None if fit data is missing (if we are here, the EOS points
                    # are there, so it means that the fit failed). I will still plot the
                    # points
                    compare_eos_fit_energy = None
                else:
                    scaling_compare_plugin = get_volume_scaling_to_formula_unit(
                        compare_plugin_data['num_atoms_in_sim_cell'][f'{element}-{configuration}'],
                        element, configuration
                    )

                    compare_eos_fit_energy = birch_murnaghan(
                        V=dense_volumes,
                        E0=ref_BM_fit_data['E0'] / scaling_ref_plugin, ## IMPORTANT! here we use the E0 of the reference plugin
                        V0=compare_BM_fit_data['min_volume'] / scaling_compare_plugin,
                        B0=compare_BM_fit_data['bulk_modulus_ev_ang3'],
                        B01=compare_BM_fit_data['bulk_deriv']
                    )
            else:
                # No compare_with plugin
                compare_eos_fit_energy = None

        # Fetch stress data, so I know if I need to do two panels or only one
        stress_data = reference_plugin_data['stress_data'][f'{element}-{configuration}']
        stress_volumes = []
        hydro_stresses_GPa = []

        # After this, `volumes` and `hydro_stresses_GPa`` are empty lists if all stresses are None
        for stress_volume, stress_tensor in stress_data:
            if stress_tensor is not None:
                stress_volumes.append(stress_volume / scaling_ref_plugin)
                #1 eV/Angstrom3 = 160.21766208 GPa
                hydro_stresses_GPa.append(
                    160.21766208 * (stress_tensor[0][0] + stress_tensor[1][1] + stress_tensor[2][2])/3
                    )

        #### START Plotting ####
        if hydro_stresses_GPa:
            fig, (stress_ax, eos_ax) = pl.subplots(nrows=2, ncols=1, gridspec_kw={'height_ratios': [1, 2], 'left': 0.15, 'right': 0.95}, sharex=True)
        else:
            # Only EOS panel
            fig, eos_ax = pl.subplots(nrows=1, ncols=1, gridspec_kw={'left': 0.15, 'right': 0.95})

        # Plot EOS: this will be done anyway
        eos_ax.plot(volumes, energies, 'ob', label=f'{PLUGIN_NAME} EOS data')
        if reference_eos_fit_energy is not None:
            eos_ax.plot(dense_volumes, reference_eos_fit_energy, '-b', label=f'{PLUGIN_NAME} fit (residuals: {residuals:.3g})')
            eos_ax.axvline(ref_BM_fit_data['min_volume'] / scaling_ref_plugin, linestyle='--', color='gray')
            if compare_eos_fit_energy is not None:
                eos_ax.plot(dense_volumes, compare_eos_fit_energy, '-r', label=f'{compare_with} fit')
                eos_ax.fill_between(dense_volumes, reference_eos_fit_energy, compare_eos_fit_energy, alpha=0.5, color='red')
        
        eos_ax.legend(loc='upper center')
        eos_ax.set_xlabel("Cell volume per formula unit ($\\AA^3$)")
        eos_ax.set_ylabel("$E - TS$ per formula unit (eV)")

        LIGHTYELLOW = (255/255, 244/255, 214/255)
        LIGHTORANGE = (255/255, 205/255, 171/255)
        if residuals is None:
            eos_ax.set_facecolor(LIGHTYELLOW)
        elif residuals > RESIDUALS_THRESHOLD:
            eos_ax.set_facecolor(LIGHTORANGE)

        conf_nice = get_conf_nice(configuration)
        fig.suptitle(f"{element} ({conf_nice})")

        # Plot stress, but only if there is data! (otherwise stress_ax is not even defined)
        if hydro_stresses_GPa:
            stress_ax.axhline(0.)
            stress_ax.plot(stress_volumes, hydro_stresses_GPa, 'o')

            # Quadratic fit (the linear one is typically not enough);
            a, b, c = np.polyfit(stress_volumes, hydro_stresses_GPa, 2)
            stress_ax.plot(dense_volumes, a * dense_volumes**2 + b * dense_volumes + c)
            # The quadratic fit leads to two solutions for zero stress, we choose the one within the volume range
            zero_stress_sol_1 = (-b - np.sqrt(b**2 - 4 * a * c))/2/a
            if zero_stress_sol_1 < max(stress_volumes) and zero_stress_sol_1 > min(stress_volumes):
                stress_ax.axvline((-b - np.sqrt(b**2 - 4 * a * c))/2/a, linestyle='--', color='gray')
            else:
                 stress_ax.axvline((-b + np.sqrt(b**2 - 4 * a * c))/2/a, linestyle='--', color='gray')
 
            stress_ax.set_ylabel("Volumetric stress (GPa)")

        pl.savefig(f"{PLOT_FOLDER}/{element}-{configuration.replace('/', '_')}.pdf")
        pl.close(fig)

    print(f"Plots written to: '{PLOT_FOLDER}'")
