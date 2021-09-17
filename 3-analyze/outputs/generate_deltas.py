#!/usr/bin/env python
import json
import os
import sys

import numpy as np
import pylab as pl
import tqdm
#import matplotlib
import matplotlib.pyplot as plt


def get_conf_nice(configuration_string):
    """Convert the configuration string to a nicely typeset string in LaTeX."""
    ret_pieces = []
    for char in configuration_string:
        if char in "0123456789":
            ret_pieces.append(f"$_{char}$")
        else:
            ret_pieces.append(char)
    return "".join(ret_pieces)


def calcDelta(v0w, b0w, b1w, v0f, b0f, b1f):
    """
    Calculate the Delta value, function copied from the official DeltaTest repository.
    I don't understand what it does, but it works.
    """

    Vi = 0.94 * (v0w + v0f) / 2.
    Vf = 1.06 * (v0w + v0f) / 2.

    a3f = 9. * v0f**3. * b0f / 16. * (b1f - 4.)
    a2f = 9. * v0f**(7. / 3.) * b0f / 16. * (14. - 3. * b1f)
    a1f = 9. * v0f**(5. / 3.) * b0f / 16. * (3. * b1f - 16.)
    a0f = 9. * v0f * b0f / 16. * (6. - b1f)

    a3w = 9. * v0w**3. * b0w / 16. * (b1w - 4.)
    a2w = 9. * v0w**(7. / 3.) * b0w / 16. * (14. - 3. * b1w)
    a1w = 9. * v0w**(5. / 3.) * b0w / 16. * (3. * b1w - 16.)
    a0w = 9. * v0w * b0w / 16. * (6. - b1w)

    x = [0, 0, 0, 0, 0, 0, 0]

    x[0] = (a0f - a0w)**2
    x[1] = 6. * (a1f - a1w) * (a0f - a0w)
    x[2] = -3. * (2. * (a2f - a2w) * (a0f - a0w) + (a1f - a1w)**2.)
    x[3] = -2. * (a3f - a3w) * (a0f - a0w) - 2. * (a2f - a2w) * (a1f - a1w)
    x[4] = -3. / 5. * (2. * (a3f - a3w) * (a1f - a1w) + (a2f - a2w)**2.)
    x[5] = -6. / 7. * (a3f - a3w) * (a2f - a2w)
    x[6] = -1. / 3. * (a3f - a3w)**2.

    y = [0, 0, 0, 0, 0, 0, 0]

    y[0] = (a0f + a0w)**2 / 4.
    y[1] = 3. * (a1f + a1w) * (a0f + a0w) / 2.
    y[2] = -3. * (2. * (a2f + a2w) * (a0f + a0w) + (a1f + a1w)**2.) / 4.
    y[3] = -(a3f + a3w) * (a0f + a0w) / 2. - (a2f + a2w) * (a1f + a1w) / 2.
    y[4] = -3. / 20. * (2. * (a3f + a3w) * (a1f + a1w) + (a2f + a2w)**2.)
    y[5] = -3. / 14. * (a3f + a3w) * (a2f + a2w)
    y[6] = -1. / 12. * (a3f + a3w)**2.

    Fi = np.zeros_like(Vi)
    Ff = np.zeros_like(Vf)

    Gi = np.zeros_like(Vi)
    Gf = np.zeros_like(Vf)

    for n in range(7):
        Fi = Fi + x[n] * Vi**(-(2. * n - 3.) / 3.)
        Ff = Ff + x[n] * Vf**(-(2. * n - 3.) / 3.)

        Gi = Gi + y[n] * Vi**(-(2. * n - 3.) / 3.)
        Gf = Gf + y[n] * Vf**(-(2. * n - 3.) / 3.)

    Delta = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi))
    #Deltarel = 100. * np.sqrt((Ff - Fi) / (Gf - Gi))
    #vref = 30.
    #bref = 100. * 10.**9. / 1.602176565e-19 / 10.**30. #100 GPa in ev_ang3
    #Delta1 = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi)) \
    #    / (v0w + v0f) / (b0w + b0f) * 4. * vref * bref

    return Delta  #, Deltarel, Delta1



if __name__ == "__main__":

    DELTA_FOLDER = 'deltas'
    os.makedirs(DELTA_FOLDER, exist_ok=True)

    all_args = sys.argv[1:]

    if not all_args:
        print("The plugin's names whose deltas will be compared must be listed explicitely as script arguments.")
        sys.exit(1)
    
    data={}
    for PLUGIN_NAME in all_args:
        try:
            with open(f'results-{PLUGIN_NAME}.json') as fhandle:
                data[PLUGIN_NAME] = json.load(fhandle)
        except OSError:
            print(f"No data found for plugin '{PLUGIN_NAME}'. A file 'results-{PLUGIN_NAME}' must be in this folder.")
            sys.exit(1)
 
    #Create a set with all the systems calculated, only the materials whose birch murnaghan fit is reported
    all_systems = set([])
    for plug_name, plug_data in data.items():
        all_systems.update(set(plug_data['BM_fit_data'].keys()))

    progress_bar = tqdm.tqdm(sorted(all_systems))
    
    #Loop over systems
    for element_and_configuration in progress_bar:
        progress_bar.set_description(f"{element_and_configuration:12s}")
        progress_bar.refresh()
    
        element, configuration = element_and_configuration.split('-')
        list_plugins_1 = []
        ind_1=0
        #This matrix will gather the results, I initialize it as a big matrix
        collect = [ [ 0 for i in range(50) ] for j in range(50) ]

        #loop over plugin 1
        for plug_name_1, plug_data_1 in data.items():
            try:
                BM_fit_data_1 = plug_data_1['BM_fit_data'][f'{element}-{configuration}']
            except KeyError:
                continue

            if not BM_fit_data_1:
                continue

            list_plugins_1.append(plug_name_1)
            ind_2 = 0
            #list_plugins_2 = []

            #loop over plugin 2
            for plug_name_2, plug_data_2 in data.items():
                
                if plug_name_2 == plug_name_1:
                    #list_plugins_2.append(plug_name_2)
                    collect[ind_1][ind_2] = 0.00
                    ind_2 = ind_2 + 1
                    continue

                try:
                    BM_fit_data_2 = plug_data_2['BM_fit_data'][f'{element}-{configuration}']
                except KeyError:
                    continue

                if not BM_fit_data_2:
                    continue

                V0_1 = BM_fit_data_1['min_volume']
                B0_1 = BM_fit_data_1['bulk_modulus_ev_ang3']
                B01_1 = BM_fit_data_1['bulk_deriv']
                V0_2 = BM_fit_data_2['min_volume']
                B0_2 = BM_fit_data_2['bulk_modulus_ev_ang3']
                B01_2 = BM_fit_data_2['bulk_deriv']

                delta = calcDelta(V0_1, B0_1, B01_1, V0_2, B0_2, B01_2)

                #list_plugins_2.append(plug_name_2)
                collect[ind_1][ind_2] = delta.round(2)
                ind_2 = ind_2 + 1

            ind_1 = ind_1 + 1
                
        to_plot = []
        for i in range(len(list_plugins_1)):
            to_plot.append([collect[i][j] for j in range(len(list_plugins_1))])

        to_plot = np.array(to_plot)

        fig, ax = plt.subplots()
        im = ax.imshow(to_plot,cmap="Reds")

        # We want to show all ticks...
        ax.set_xticks(np.arange(len(list_plugins_1)))
        ax.set_yticks(np.arange(len(list_plugins_1)))
        # ... and label them with the respective list entries
        ax.set_xticklabels(list_plugins_1)
        ax.set_yticklabels(list_plugins_1)

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(list_plugins_1)):
            for j in range(len(list_plugins_1)):
                text = ax.text(j, i, to_plot[i, j], ha="center", va="center", color="black")

        conf_nice = get_conf_nice(configuration)
        ax.set_title(f"{element} ({conf_nice})")
        fig.tight_layout()
        fig.savefig(f"{DELTA_FOLDER}/{element}-{configuration}.png")
        plt.close()

