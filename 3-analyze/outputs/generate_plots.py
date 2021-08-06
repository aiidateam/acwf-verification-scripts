#!/usr/bin/env python
import json
import os
import sys

import numpy as np
import pylab as pl
import tqdm

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


def get_conf_nice(configuration_string):
    """Convert the configuration string to a nicely typeset string in LaTeX."""
    ret_pieces = []
    for char in configuration_string:
        if char in "0123456789":
            ret_pieces.append(f"$_{char}$")
        else:
            ret_pieces.append(char)
    return "".join(ret_pieces)


def birch_murnaghan(V,E0,V0,B0,B01):
    r = (V0/V)**(2./3.)
    return (E0 +
            9./16. * B0 * V0 * (
            (r-1.)**3 * B01 + 
            (r-1.)**2 * (6. - 4.* r)))
   

if __name__ == "__main__":
    try:
        compare_with = sys.argv[1]
    except IndexError:
        compare_with = None

    try:
        with open(f'results-{PLUGIN_NAME}.json') as fhandle:
            reference_plugin_data = json.load(fhandle)
    except OSError:
        print(f"No data found for your plugin '{PLUGIN_NAME}'. Did you run `./get_results.py` first?")
        sys.exit(1)
    
    if compare_with is None:
        print(f"Plotting data for plugin '{PLUGIN_NAME}' only.")
        print("If you want to compare, pass another parameter with the plugin to compare with.")
        compare_plugin_data = None
    else:
        print(f"Plotting data for plugin '{PLUGIN_NAME}' compared with '{compare_with}'.")
        try:
             with open(f'results-{compare_with}.json') as fhandle:
                compare_plugin_data = json.load(fhandle)
        except OSError:
            print(f"No data found for the reference plugin '{compare_with}': you need the file results-{compare_with}.json.")
            sys.exit(1)

    if compare_with is None:
        PLOT_FOLDER = f'plots-{PLUGIN_NAME}'
    else:
        PLOT_FOLDER = f'plots-{PLUGIN_NAME}-vs-{compare_with}'
    os.makedirs(PLOT_FOLDER, exist_ok=True)

    all_systems = set(reference_plugin_data['eos_data'].keys())
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

        # Get the x axis for the plot
        volumes, energies = (np.array(eos_data).T).tolist()
        dense_volumes = np.linspace(min(volumes), max(volumes), 100)

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
        else:
            reference_eos_fit_energy = birch_murnaghan(
                V=dense_volumes,
                E0=ref_BM_fit_data['E0'],
                V0=ref_BM_fit_data['min_volume'],
                B0=ref_BM_fit_data['bulk_modulus_ev_ang3'],
                B01=ref_BM_fit_data['bulk_deriv']
            )

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
                    compare_eos_fit_energy = birch_murnaghan(
                        V=dense_volumes,
                        E0=ref_BM_fit_data['E0'], ## IMPORTANT! here we use the E0 of the reference plugin
                        V0=compare_BM_fit_data['min_volume'],
                        B0=compare_BM_fit_data['bulk_modulus_ev_ang3'],
                        B01=compare_BM_fit_data['bulk_deriv']
                    )
            else:
                # No compare_with plugin
                compare_eos_fit_energy = None

        # Plotting
        fig = pl.figure()

        # If we are here, 
        pl.plot(volumes, energies, 'ob', label=f'{PLUGIN_NAME} EOS data')
        
        if reference_eos_fit_energy is not None:
            pl.plot(dense_volumes, reference_eos_fit_energy, '-b', label=f'{PLUGIN_NAME} fit')
            if compare_eos_fit_energy is not None:
                pl.plot(dense_volumes, compare_eos_fit_energy, '-r', label=f'{compare_with} fit')
                pl.fill_between(dense_volumes, reference_eos_fit_energy, compare_eos_fit_energy, alpha=0.5, color='red')
        
        pl.legend(loc='upper center')
        pl.xlabel("Cell volume ($\\AA^2$)")
        pl.ylabel("$E_{tot}$ (eV)")

        conf_nice = get_conf_nice(configuration)
        pl.title(f"{element} ({conf_nice})")
        pl.savefig(f"{PLOT_FOLDER}/{element}-{configuration}.png")
        pl.close(fig)

    print(f"Plots written to: '{PLOT_FOLDER}'")
