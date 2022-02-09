#!/usr/bin/env python
import json
import os
import sys

import numpy as np
import pylab as pl
from scipy.optimize import curve_fit


BINS = 100
FORM_ENERGY_PRINT_THRESHOLD = 0.03 # eV/atom

OUT_FOLDER = 'formation_energies_output'

def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


def generate_plots(plugin1, plugin2, x_zoom_factor=1., abs_x_range=None):

    # x_zoom_factor: Adapt this factor to change the zoom on the x axis
    # The default zoom is obtained from the standard deviation of the data
    # This is a multiplicative number; a number > 1 means zoom in, a number < 1 means zoom out

    # if abs_x_range is passed, x_zoom_factor is ignored and data is plotted in the range
    # [-abs_x_range, +abs_x_range]

    fname = f'formation-energies-{plugin1}-VS-{plugin2}.json'
    
    try:
        with open(f'{OUT_FOLDER}/{fname}') as fhandle:
            formation_data = json.load(fhandle)['formation_energies']
    except OSError:
        print(f"No file '{fname}' found in the '{OUT_FOLDER}' subfolder! Run ./compute_formation_energies.py first.")
        sys.exit(1)

    formation_energy_dissimilarities = []
    for system, data in formation_data.items():
        # System is something like "Ac-X2O/X2O3/XO3"
        plugin1_formation_energy = data[plugin1]
        plugin2_formation_energy = data[plugin2]
        if plugin1_formation_energy is None or plugin2_formation_energy is None:
            # Deal with missing data for at least one plugin
            continue
        formation_energy_dissimilarities.append((plugin2_formation_energy - plugin1_formation_energy, system, plugin1_formation_energy, plugin2_formation_energy))

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

    flat_data = np.array([_[0] for _ in formation_energy_dissimilarities])
    
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
    pl.title(f"{plugin1} VS {plugin2}")
    pl.xlim(-half_range, half_range)
    pl.tight_layout()
    pl.savefig(f"{OUT_FOLDER}/histogram-{plugin1}-VS-{plugin2}.png")
    pl.close(fig)

    abs_formation_energy_dissimilarities = [
        (abs(_[0]), _[1], _[2], _[3]) for _ in formation_energy_dissimilarities]
    
    # Sort in-place from the largest dissimilarity (in abs value)
    abs_formation_energy_dissimilarities.sort()
    abs_formation_energy_dissimilarities.reverse()

    with open(f"{OUT_FOLDER}/discrepancies-{plugin1}-VS-{plugin2}.txt", "w") as fhandle:
        fhandle.write(f"## {plugin1} VS {plugin2}\n")
        fhandle.write(f"## Cases with abs(Formation energies) > {FORM_ENERGY_PRINT_THRESHOLD} eV/atom:\n")
        for form_energy_dissimilarity, system, form_energy_plugin1, form_energy_plugin2 in abs_formation_energy_dissimilarities:
            if form_energy_dissimilarity < FORM_ENERGY_PRINT_THRESHOLD:
                # Here I assume I already sorted them
                break
            fhandle.write(f"{system:30s}: {form_energy_dissimilarity:.6f} ({form_energy_plugin1} vs {form_energy_plugin2})\n")

if __name__ == "__main__":
    try:
        plugin1 = sys.argv[1]
        plugin2 = sys.argv[2]
        zoom = sys.argv[3]
    except IndexError:
        print("Pass as two parameters the two plugins to compare, and a third parameter with the zoom factor")
        sys.exit(1)
    try:
        zoom = float(zoom)
    except IndexError:
        print("The zoom factor must be a float number")
        sys.exit(1)

    generate_plots(plugin1, plugin2, zoom)
    print("*** IMPORTANT WARNING!!! THIS SCRIPT USES (STILL) THE EOS (ENERGIES AND VOLUMES) FROM THE SIMULATION CELL, AND NOT FROM PER FORMULA UNIT!! ***")