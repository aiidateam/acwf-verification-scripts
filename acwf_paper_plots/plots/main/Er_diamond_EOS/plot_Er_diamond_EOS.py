#!/usr/bin/env python

# The version of the script will be placed in the json file containing the results.
# We should change this number anytime this script or `eos_utils.eosfit_31_adapted` is modified.
import json
import numpy as np
import pylab as pl
from scipy.interpolate import interp1d
from scipy.optimize import minimize

smearing_name_mapping = {
    'fermi-dirac': 'FD',
    'cold': 'Cold',
    'gaussian': 'Gaussian'
}

markers_dict = {
    'fermi-dirac': '-',
    'cold': ':',
    'gaussian': '--',
}

class ColorManager:
    RED_SCALE = ['#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000'][::-1]
    BLUE_SCALE = ['#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b'][::-1]
    GREEN_SCALE = ['#f7fcf5','#e5f5e0','#c7e9c0','#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b'][::-1]
    def __init__(self):
        self.red_cnt = 0
        self.blue_cnt = 0
        self.green_cnt = 0

    def get_next_color(self, smearing_type):
        if smearing_type == 'cold':
            scale = self.BLUE_SCALE
            color = scale[self.blue_cnt]
            self.blue_cnt += 1
        elif smearing_type == 'gaussian':
            scale = self.GREEN_SCALE
            color = scale[self.green_cnt]
            self.green_cnt += 1
        elif smearing_type == 'fermi-dirac':
            scale = self.RED_SCALE
            color = scale[self.red_cnt]
            self.red_cnt += 1
        else:
            raise ValueError("Unknown smearing type")
        return color


if __name__ == "__main__":
    num_atoms_per_cell = 2 # Hardcoded as I know this is Er diamond

    fig = pl.figure()
    color_manager = ColorManager()

    with open('Er-diamond-EOS-data.json') as fhandle:
        data = json.load(fhandle)

    min_y = None
    all_plots = []
    warning_lines = []
    for curve in data:
        smearing_type = curve['smearing_type']
        degauss_meV = curve['degauss_meV']
        volumes = curve['volumes']
        energies = curve['energies']

        marker = markers_dict[smearing_type]
        color = color_manager.get_next_color(smearing_type)

        # smooth interpolation between points
        interpolation_f = interp1d(volumes, energies, kind='cubic')

        # Left minimum
        bounds = [(25 * num_atoms_per_cell, 45 * num_atoms_per_cell)]
        epsilon = 1.
        minimization = minimize(interpolation_f, (30. * num_atoms_per_cell,), bounds = bounds, method='Nelder-Mead')
        if minimization.success and minimization.x[0] > bounds[0][0] + epsilon and minimization.x[0] < bounds[0][1] - epsilon:
            left_min = minimization.x[0]
            left_min_y = interpolation_f(left_min)
        else:
            left_min = None
            left_min_y = None
        
        # Right minimum
        bounds = [(50 * num_atoms_per_cell, 95 * num_atoms_per_cell)]
        minimization = minimize(interpolation_f, (90. * num_atoms_per_cell,), bounds = bounds, method='Nelder-Mead')
        if minimization.success and minimization.x[0] > bounds[0][0] + epsilon and minimization.x[0] < bounds[0][1] - epsilon:
            right_min = minimization.x[0]
            right_min_y = interpolation_f(right_min)
        else:
            right_min = None
            right_min_y = None

        # Get data to plot on dense grid
        volumes_dense = np.linspace(min(volumes), max(volumes), 1000)
        energies_dense = interpolation_f(volumes_dense)

        further_args = (marker,)
        kwargs = {
            'color': color,
            'label': f"{smearing_name_mapping[smearing_type]} ({int(degauss_meV)} meV)"
        }

        # IMPORTANT! These quantities are stil *per simulation cell*
        all_plots.append(
            (volumes_dense, energies_dense, further_args, kwargs, left_min, left_min_y, right_min, right_min_y))
        if min_y is None:
            min_y = np.array(energies_dense).min()
        else:
            min_y = min(min_y, np.array(energies_dense).min())

    # Plot all, shifting on y AND PLOTTING QUANTITIES PER ATOM AND NOT PER UNIT CELL
    for x, y, further_args, kwargs, left_min, left_min_y, right_min, right_min_y in all_plots:
        pl.plot(x / num_atoms_per_cell, (y - min_y)  / num_atoms_per_cell, *further_args, **kwargs)
        if left_min is not None:
            pl.plot([left_min / num_atoms_per_cell], [(left_min_y - min_y) / num_atoms_per_cell], 'x', color=kwargs['color'])
        if right_min is not None:
            pl.plot([right_min / num_atoms_per_cell], [(right_min_y - min_y) / num_atoms_per_cell], 'x', color=kwargs['color'])

    pl.xlabel('Volume per atom ($\AA^3$)')
    pl.ylabel('Free energy per atom (eV)')
    pl.xlim(volumes_dense[0]/ num_atoms_per_cell, volumes_dense[-1]/ num_atoms_per_cell)

    pl.ylim(0, 0.7)
    pl.legend(loc='upper center', ncol=3, fontsize=9)

    pl.savefig('Er-diamond-vs-smearing.pdf')
