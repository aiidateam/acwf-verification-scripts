#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 14:49:59 2022

@author: wolloch
"""
import json
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import tqdm
from eosfit_31_adapted import BM
from quantities_for_comparison import birch_murnaghan

eV_over_ang3_to_GPa = 160.21766208


def make_histograms(stats, bins, title, outlier_range=None, outlier_what=None, raw_data=None):
    """
    outlier_what should a valid key inside the raw data 'BM_fit_data' dictionaries, e.g. 'bulk_deriv'
    """
    B0V0 = []
    B1B0 = []
    B1V0 = []
    for k, v in stats.items():
        if outlier_range is not None and (
            raw_data['BM_fit_data'][k][outlier_what] < outlier_range[0] or
            raw_data['BM_fit_data'][k][outlier_what] > outlier_range[1]):
                continue
        B0V0.append(v['mean_B0']/v['mean_V0'])
        B1B0.append(v['mean_B1']/v['mean_B0'])
        B1V0.append(v['mean_B1']/v['mean_V0'])
    B0V0_hist, B0V0_bin_edges = np.histogram(B0V0,bins=bins, range=(0, 100))
    B1V0_hist, B1V0_bin_edges = np.histogram(B1V0,bins=bins, range=(0, 100))
    B1B0_hist, B1B0_bin_edges = np.histogram(B1B0,bins=bins, range=(0, 100))
    B0V0_top_bin = []
    B1V0_top_bin = []
    B1B0_top_bin = []
    for k, v in stats.items():
        if v['mean_B0']/v['mean_V0'] > B0V0_bin_edges[-2]:
            B0V0_top_bin.append(k)
        if v['mean_B1']/v['mean_V0'] > B1V0_bin_edges[-2]:
            B1V0_top_bin.append(k)
        if v['mean_B1']/v['mean_B0'] > B1B0_bin_edges[-2]:
            B1B0_top_bin.append(k)
    print(f'Systems with B0:V0 in largest bin: {B0V0_top_bin}\n'
          f'Systems with B1:V0 in largest bin: {B1V0_top_bin}\n'
          f'Systems with B1:B0 in largest bin: {B1B0_top_bin}\n')

    outlier_str = ""
    if outlier_range is not None:
        outlier_str = "-outliers-{outlier_what}-{outlier_range[0]}-{outlier_range[1]}"

    plt.hist(B0V0,bins=bins, range=(0, 100))
    plt.xlabel(r'absDevB0/absDevV0')
    plt.savefig(f'B0V0_{title}_hist{outlier_str}.png',
                dpi=300,
                orientation='landscape',
                bbox_inches='tight')
    plt.clf()
    plt.hist(B1V0,bins=bins, range=(0, 1000))
    plt.xlabel(r'absDevB1/absDevV0')
    plt.savefig(f'B1V0_{title}_hist{outlier_str}.png',
                dpi=300,
                orientation='landscape',
                bbox_inches='tight')
    plt.clf()
    plt.hist(B1B0,bins=bins, range=(0, 100))
    plt.xlabel(r'absDevB1/absDevB0')
    plt.savefig(f'B1B0_{title}_hist{outlier_str}.png',
                dpi=300,
                orientation='landscape',
                bbox_inches='tight')
    plt.clf()
    out_dict = {'B0V0': {'hist': B0V0_hist,
                         'bin_edges': B0V0_bin_edges,
                         'in_top_bin': B0V0_top_bin},
                'B1V0': {'hist': B1V0_hist,
                         'bin_edges': B1V0_bin_edges,
                         'in_top_bin': B1V0_top_bin},
                'B1B0': {'hist': B1B0_hist,
                         'bin_edges': B1B0_bin_edges,
                         'in_top_bin': B1B0_top_bin},}
    return out_dict


def guess_initial(EV):
    en=EV[:,1]
    vol=EV[:,0]
    a, b, c = np.polyfit(vol, en, 2)
    v0 = -b / (2 * a)
    e0 = a * (v0 ** 2) + b * v0 + c
    b0 = 2 * a * v0
    b1 = 4  # b1 is usually a small number like 4
    return [e0, v0, b0, b1]


def fit_eos(EV, eos_type='BM'):
    en=EV[:,1]
    vol=EV[:,0]
    initial_guess = guess_initial(EV)
    #initial_guess = [min(en), min(vol), 1, 1]
    if eos_type == 'AiiDA':
        V0, E0, B0, B1, res = BM(EV)
        out_dict = {'energy_vs_volume_array': EV,
                    'fitted_parameters': {'total_energy': E0,
                                          'equilibrium_volume': V0,
                                          'bulk_modulus': B0,
                                          'bulk_modulus_derivative': B1},
                    'residuals': res,
                    'bulk_modulus_in_GPa': B0*eV_over_ang3_to_GPa,
                    'equilibrium_vol_in_A3': V0}
        return out_dict
    elif eos_type == 'BM':
        popt, pcov = curve_fit(birch_murnaghan_EOS, vol, en, p0=initial_guess)
    elif eos_type == 'vignet':
        popt, pcov = curve_fit(vignet_EOS, vol, en, p0=initial_guess)
    else:
        print('incorrect eos_type')
        return
    if len(EV) == 4:
        perr = [0.0, 0.0, 0.0, 0.0] # only four inputs make fit exact!
    else:
        perr = np.sqrt(np.diag(pcov)) # one standard deviation errors of fitted params
    out_dict = {'energy_vs_volume_array': EV,
                'fitted_parameters': {'total_energy': popt[0],
                                      'equilibrium_volume': popt[1],
                                      'bulk_modulus': popt[2],
                                      'bulk_modulus_derivative': popt[3]},
                'estimated_errors': {'total_energy': perr[0],
                                     'equilibrium_volume': perr[1],
                                     'bulk_modulus': perr[2],
                                     'bulk_modulus_derivative': perr[3]},
                'bulk_modulus_in_GPa': popt[2]*eV_over_ang3_to_GPa,
                'equilibrium_vol_in_A3': popt[1]}
    return out_dict

def perturb_en(data, noise_sigma):
    l = len(data)
    noise = np.random.normal(0, noise_sigma, l)
    noisy_data = np.array(data) + noise
    return noisy_data

def get_statistics(dataset, noise_sigma, volumes_percents, nr_of_samples=10000):
    deviations = {}
    stats = {}
    progress_bar = tqdm.tqdm(sorted(dataset['BM_fit_data'].items()))
    for key, val in progress_bar:
        deviations[key] = {'bulk_modulus': [],
                           'bulk_modulus_derivative': [],
                           'equilibrium_volume': [],
                           'total_energy': [],
                           'failed_runs': 0}
        volumes_list = [i*val['min_volume'] for i in volumes_percents]
        en_ok_murn = birch_murnaghan(np.array(volumes_list), 0, val['min_volume'], val["bulk_modulus_ev_ang3"], val['bulk_deriv'])
        data_diff_vols = np.array([[volumes_list[i], en_ok_murn[i]] for i in range(len(en_ok_murn))])
        new_fit_with_different_vols = fit_eos(data_diff_vols, eos_type='AiiDA')
        print(val['min_volume'], new_fit_with_different_vols['fitted_parameters']['equilibrium_volume'])
        for i in range(nr_of_samples):
            new_en = perturb_en(en_ok_murn, noise_sigma)
            data = np.array([[volumes_list[i], new_en[i]] for i in range(len(new_en))])
            try:
                eos = fit_eos(data, eos_type='AiiDA')
                for k, v in eos['fitted_parameters'].items():
                    dev = 100*(new_fit_with_different_vols['fitted_parameters'][k]-eos['fitted_parameters'][k])/eos['fitted_parameters'][k]
                    deviations[key][k].append(dev)
            except:
                deviations[key]['failed_runs'] += 1
        stats[key] = {'mean_V0': np.mean(abs(np.asarray(deviations[key]['equilibrium_volume']))),
                      'mean_B0': np.mean(abs(np.asarray(deviations[key]['bulk_modulus']))),
                      'mean_B1': np.mean(abs(np.asarray(deviations[key]['bulk_modulus_derivative']))),
                      'mean_E0': np.mean(abs(np.asarray(deviations[key]['total_energy'])))}
    return stats, deviations


if __name__ == '__main__':
    with open('results-oxides-verification-PBE-v1-ae.json') as fhandle:
        oxides = json.load(fhandle)
    with open('results-unaries-verification-PBE-v1-ae.json') as fhandle:
        unaries = json.load(fhandle)
            
    noise_sigma = 1E-4
    nr_of_samples = 100
    nr_volume_points = 7
    max_and_min_percentage = [0.94, 1.06]

    interval = (max_and_min_percentage[1] - max_and_min_percentage[0]) / (nr_volume_points-1)
    set_vols_perc = [max_and_min_percentage[0]+interval*i for i in range(nr_volume_points)]

    print(set_vols_perc)

    unaries_title = f'unaries_ae_{noise_sigma}'
    oxides_title = f'oxides_ae_{noise_sigma}'

    stats_unaries, deviations_unaries = get_statistics(unaries, noise_sigma, set_vols_perc, nr_of_samples)
    stats_oxides, deviations_oxides = get_statistics(oxides, noise_sigma, set_vols_perc, nr_of_samples)

    print('make histograms:')
    make_histograms(stats=stats_unaries, bins=50, title=unaries_title)
    make_histograms(stats=stats_oxides, bins=50, title=oxides_title)

