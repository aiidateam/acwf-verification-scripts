#!/usr/bin/env python
import json
import os
import sys

import numpy as np
import pylab as pl
from scipy.optimize import curve_fit

BINS = 100
PRINT_THRESHOLD = 0.01 # eV/atom
VERBOSE = True

OUT_FOLDER = 'formation-energies-output'

def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


def get_dissimilarities(plugin1, plugin2, what):

    fname = f'formation-energies-{plugin1}-VS-{plugin2}.json'
    try:
        with open(f'{OUT_FOLDER}/{fname}') as fhandle:
            raw_data = json.load(fhandle)
    except OSError:
        print(f"No file '{fname}' found in the '{OUT_FOLDER}' subfolder! Run ./compute_formation_energies.py first.")
        sys.exit(1)

    dissimilarities = []
    missing_count = 0

    if what == 'formation-energy':
        formation_data = raw_data['formation_energies']
        for system, data in formation_data.items():
            # System is something like "Ac-X2O|X2O3|XO3"
            plugin1_formation_energy = data[plugin1]
            plugin2_formation_energy = data[plugin2]
            if plugin1_formation_energy is None or plugin2_formation_energy is None:
                # Deal with missing data for at least one plugin
                if VERBOSE:
                    print(f"WARNING: missing {system}")
                missing_count += 1                
                continue
            dissimilarities.append((plugin2_formation_energy - plugin1_formation_energy, system, plugin1_formation_energy, plugin2_formation_energy))
    elif what == 'unaries':
        unaries_data = raw_data['unaries_energy_difference']
        for element, data in unaries_data.items():
            reference_configuration = data['reference_configuration']
            for configuration, plugin_data in data['configurations'].items():
                if configuration == reference_configuration:
                    # Avoid to store zeros for the reference configuration, that would
                    # bias the final plot
                    continue
                plugin1_energy_difference = plugin_data[plugin1]
                plugin2_energy_difference = plugin_data[plugin2]
                if plugin1_energy_difference is None or plugin2_energy_difference is None:
                    if VERBOSE:
                        print(f"WARNING: missing {element} - {configuration}")
                    # Deal with missing data for at least one plugin
                    missing_count += 1
                    continue
                dissimilarities.append(
                    (plugin2_energy_difference - plugin1_energy_difference,
                    f"{element}-{configuration}-wrt-{reference_configuration}",
                    plugin1_energy_difference,
                    plugin2_energy_difference))
    else:
        raise ValueError(f"Unknown value of 'what': '{what}'")

    if missing_count:
        print(f"WARNING: {missing_count} systems missing when checking data for '{what}'")

    return dissimilarities

def generate_plots(plugin1, plugin2, what, x_zoom_factor=1., abs_x_range=None):

    # x_zoom_factor: Adapt this factor to change the zoom on the x axis
    # The default zoom is obtained from the standard deviation of the data
    # This is a multiplicative number; a number > 1 means zoom in, a number < 1 means zoom out

    # if abs_x_range is passed, x_zoom_factor is ignored and data is plotted in the range
    # [-abs_x_range, +abs_x_range]


    # Plotting
    fig = pl.figure(figsize=(18,6))

    TINY_SIZE = 18
    SMALL_SIZE = 20
    MEDIUM_SIZE = 24
    BIGGER_SIZE = 28

    pl.rc('font', size=SMALL_SIZE)# controls default text sizes
    pl.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
    pl.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    pl.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    pl.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    pl.rc('legend', fontsize=TINY_SIZE)    # legend fontsize
    pl.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    dissimilarities = get_dissimilarities(plugin1, plugin2, what=what)

    flat_data = np.array([_[0] for _ in dissimilarities])
    
    if abs_x_range is None:
        half_range = np.sqrt(np.mean(np.array(flat_data)**2)) / x_zoom_factor
    else:
        half_range = abs_x_range

    still_on_left = len(flat_data[flat_data < -half_range])
    still_on_right = len(flat_data[flat_data > half_range])

    hist_y, bins, patches = pl.hist(
        flat_data, bins=BINS, range=[-half_range, half_range],
        label=f"{still_on_left} more on the left, {still_on_right} more on the right"
    )

    # Fit Gaussian and plot it
    hist_x = (bins[1:] + bins[:-1])/2
    popt, pcov = curve_fit(gaussian, hist_x, hist_y, p0=[10., 0., 0.001])
    x = np.linspace(pl.xlim()[0], pl.xlim()[1], 1000)
    sigma = abs(popt[2])
    ## NOTES ON THE RELATION BETWEEN THE SIGMA OF THE GAUSSIAN AND THE FWHM
    #  np.exp(-HWHM**2/(2*sigma**2)) = 1/2
    #  -HWHM**2/(2*sigma**2) = ln(1/2)
    #  HWHM**2/(2*sigma**2) = ln(2)
    #  HWHM**2 = ln(2) * (2*sigma**2)
    #  HWHM = sqrt(ln(2)) * sqrt(2) * sigma
    #  FWHM = 2*HWHM = 2*sqrt(2)*sqrt(ln(2)) * sigma
    pl.plot(x, gaussian(x, *popt), 'r:', label=rf'Gaussian fit (FWHM = {2*np.sqrt(2)*np.sqrt(np.log(2))*sigma:.5f})')
    pl.axvline(popt[1], color='r', linestyle=':')
    # Reset the xlim
    pl.xlim(x[0], x[-1])

    pl.legend(loc='upper right')
    pl.xlabel("Formation energy dissimilarity (eV/atom)")
    pl.ylabel("Frequency")
    pl.title(f"{plugin1} VS {plugin2} ({what})")
    pl.xlim(-half_range, half_range)
    pl.tight_layout()
    pl.savefig(f"{OUT_FOLDER}/histogram-{what}-{plugin1}-VS-{plugin2}.png")
    pl.close(fig)

    abs_dissimilarities = [
        (abs(_[0]), _[1], _[2], _[3]) for _ in dissimilarities]
    
    # Sort in-place from the largest dissimilarity (in abs value)
    abs_dissimilarities.sort()
    abs_dissimilarities.reverse()

    with open(f"{OUT_FOLDER}/discrepancies-{what}-{plugin1}-VS-{plugin2}.txt", "w") as fhandle:
        fhandle.write(f"## {plugin1} VS {plugin2}\n")
        if what == 'formation-energy':
            fhandle.write(f"## Cases with abs(formation energies) > {PRINT_THRESHOLD} eV/atom:\n")
        elif what == 'unaries':
            fhandle.write(f"## Cases with abs(energy difference) > {PRINT_THRESHOLD} eV/atom:\n")
        else:
            raise ValueError(f"Unknown value of 'what': '{what}'")
        for dissimilarity, system, data_plugin1, data_plugin2 in abs_dissimilarities:
            if dissimilarity < PRINT_THRESHOLD:
                # Here I assume I already sorted them
                break
            fhandle.write(f"{system:30s}: {dissimilarity:.6f} ({data_plugin1} vs {data_plugin2})\n")

if __name__ == "__main__":
    try:
        plugin1_name = sys.argv[1]
        plugin2_name = sys.argv[2]
        what = sys.argv[3]
        zoom_value = sys.argv[4:5]
        if zoom_value == []:
            zoom_value = 1.
        else:
            zoom_value = zoom_value[0]
    except IndexError:
        print("Pass as two parameters the two plugins to compare, a third parameter with 'what' to compare, and a fourth (optional) parameter as the zoom factor")
        sys.exit(1)
    try:
        zoom_value = float(zoom_value)
    except IndexError:
        print("The zoom factor must be a float number")
        sys.exit(1)

    generate_plots(plugin1_name, plugin2_name, what, zoom_value)
