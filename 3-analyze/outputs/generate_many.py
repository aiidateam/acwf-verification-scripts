#!/usr/bin/env python
import json
import os
import sys

import numpy as np
import pylab as pl
import tqdm


colors_list = ["b", "g", "r", "c", "m", "y", "k", "darkorange", "dimgrey", "sienna", "deeppink", "gold", "lime"]

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

    PLOT_FOLDER = 'plots-many'
    os.makedirs(PLOT_FOLDER, exist_ok=True)

    all_args = sys.argv[1:]

    if not all_args:
        print("The plugin's names whose results will be plotted must be listed explicitely as script arguments.")
        sys.exit(1)
    
    data={}
    colors={}
    i=0
    for PLUGIN_NAME in all_args:
        try:
            with open(f'results-{PLUGIN_NAME}.json') as fhandle:
                data[PLUGIN_NAME] = json.load(fhandle)
        except OSError:
            print(f"No data found for plugin '{PLUGIN_NAME}'. A file 'results-{PLUGIN_NAME}' must be in this folder.")
            sys.exit(1)
        colors[PLUGIN_NAME] = colors_list[i]
        i=i+1
 
    #Create a set with all the systems calculated
    all_systems = set([])
    for plug_name, plug_data in data.items():
        all_systems.update(set(plug_data['BM_fit_data'].keys()))
        try:
            all_systems.update(set(plug_data['eos_data'].keys()))
        except:
            pass

    progress_bar = tqdm.tqdm(sorted(all_systems))
    
    #Loop over systems
    for element_and_configuration in progress_bar:
        progress_bar.set_description(f"{element_and_configuration:12s}")
        progress_bar.refresh()
    
        fig = pl.figure()    

        element, configuration = element_and_configuration.split('-')

        #loop over plugins
        for plug_name, plug_data in data.items():
            try:
                eos_data = plug_data['eos_data'][f'{element}-{configuration}']
            except KeyError:
                eos_data = None
            try:
                BM_fit_data = plug_data['BM_fit_data'][f'{element}-{configuration}']
            except KeyError:
                BM_fit_data = None
            
            if eos_data is None and BM_fit_data is None:
                # If there is no data, I skip this material
                continue

            #Case when eos_data are present, the x range is set according to data
            if eos_data is not None:
                volumes, energies = (np.array(eos_data).T).tolist()
                dense_volumes = np.linspace(min(volumes), max(volumes), 100)
                #The data are printed with E min set to zero! In case E min was not obtained, skip the plot
                if BM_fit_data is not None:
                    pl.plot(volumes, np.array(energies)-BM_fit_data["E0"], 'o', label=f'{plug_name} EOS data', color=colors[plug_name]) 
                else:
                    message = "EoS data present but no fit. Curve will not be displayed since not possible to get E0"
                    print(plug_name, element_and_configuration, message)
            #Case when eos_data are not present (only BM_fit_data availables), the x range is set around V0
            else:
                dense_volumes = np.linspace(BM_fit_data["min_volume"]*0.94, BM_fit_data["min_volume"]*1.06, 100)

            #Plot the fit curve. Note: we could have the fit parameters but no eos data (e.g. cottonier-wien2k)
            if BM_fit_data is not None:
                eos_fit_energy = birch_murnaghan(
                    V=dense_volumes,
                    E0=0, #Set to zero!
                    V0=BM_fit_data['min_volume'],
                    B0=BM_fit_data['bulk_modulus_ev_ang3'],
                    B01=BM_fit_data['bulk_deriv']
                )
                pl.plot(dense_volumes, eos_fit_energy, label=f'{plug_name} fit', color=colors[plug_name])


        pl.legend(loc='upper center', fontsize=8)
        pl.xlabel("Cell volume ($\\AA^2$)")
        pl.ylabel("$E_{tot}$ (eV)")

        conf_nice = get_conf_nice(configuration)
        pl.title(f"{element} ({conf_nice})")
        pl.savefig(f"{PLOT_FOLDER}/{element}-{configuration}.png")
        pl.close(fig)


