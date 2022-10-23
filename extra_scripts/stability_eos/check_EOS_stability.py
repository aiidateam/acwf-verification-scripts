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

eV_over_ang3_to_GPa = 160.21766208

def birch_murnaghan_EOS(V, E0, V0, B0, dB0):
    return E0+9*V0*B0/16*(((V0/V)**(2./3.)-1)**3*dB0+((V0/V)**(2./3.)-1)**2*(6-4*(V0/V)**(2./3.)))

def fit_eos(EV):
    en=EV[:,1]
    vol=EV[:,0]
    initial_guess = [min(en), min(vol), 1, 1]
    popt, pcov = curve_fit(birch_murnaghan_EOS, vol, en, p0=initial_guess)
                           #bounds=([-np.inf, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]))
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

def return_noisy_eos(data, noise_sigma):
    l = len(data)
    noise = np.transpose(np.vstack((np.zeros(l),
                                    np.random.normal(0, noise_sigma, l))))
    noisy_data = data + noise
    return fit_eos(noisy_data) 

def get_statistics(dataset, noise_sigma, nr_of_samples=10000):
    deviations = {}
    stats = {}
    progress_bar = tqdm.tqdm(sorted(dataset['eos_data'].items()))
    for key, val in progress_bar:
        data = np.asarray(val)
        eos = fit_eos(data)
        deviations[key] = {'bulk_modulus': [],
                           'bulk_modulus_derivative': [],
                           'equilibrium_volume': [],
                           'total_energy': [],
                           'failed_runs': 0}
        for i in range(nr_of_samples):
            try:
                noisy_eos = return_noisy_eos(data, noise_sigma)
                for k, v in noisy_eos['fitted_parameters'].items():
                    dev = 100*(v-eos['fitted_parameters'][k])/eos['fitted_parameters'][k]
                    deviations[key][k].append(dev)
            except:
                deviations[key]['failed_runs']+=1
        stats[key] = {'mean_V0': np.mean(abs(np.asarray(deviations[key]['equilibrium_volume']))),
                      'mean_B0': np.mean(abs(np.asarray(deviations[key]['bulk_modulus']))),
                      'mean_B1': np.mean(abs(np.asarray(deviations[key]['bulk_modulus_derivative']))),
                      'mean_E0': np.mean(abs(np.asarray(deviations[key]['total_energy'])))}
    return stats, deviations

def plot_V0B0B1(stats, data, title):
    B = []
    B0_V0 = []
    B1_B0 = []
    B1_V0 = []
    for k, v in stats.items():
        B.append(data["BM_fit_data"][k]["bulk_modulus_ev_ang3"]*eV_over_ang3_to_GPa)
        B0_V0.append(v['mean_B0']/v['mean_V0'])
        B1_V0.append(v['mean_B1']/v['mean_V0'])
        B1_B0.append(v['mean_B1']/v['mean_B0'])
    
    c_B0V0, p_B0V0 = pearsonr(B, B0_V0)
    c_B1V0, p_B1V0 = pearsonr(B, B1_V0)
    
    plt.subplot(1,2,1)
    plt.semilogy(B, B0_V0, 'ro')
    plt.xlabel('bulk modulus [GPA]')
    plt.ylabel(r'absDevB0/absDevV0')
    plt.title(f'corr={np.round(c_B0V0,3)}, p-value={np.round(p_B0V0,3)}')
    plt.subplot(1,2,2)
    plt.semilogy(B, B1_V0, 'bx')
    plt.xlabel('bulk modulus [GPA]')
    plt.ylabel(r'absDevB1/absDevV0')
    plt.title(f'corr={np.round(c_B1V0,3)}, p-value={np.round(p_B1V0,3)}')
    fig = plt.gcf()
    fig.tight_layout(w_pad=1.0)    
    fig.suptitle(title, fontsize=14)
    fig.subplots_adjust(top=0.85)
    plt.savefig(f'{title}_statistics.png',
                dpi=300,
                orientation='landscape',
                bbox_inches='tight')
    return

if __name__ == '__main__':
    with open('/fs/home/wolloch/aiida/acwf-verification-scripts/3-analyze/outputs/results-oxides-verification-PBE-v1-fleur.json') as fhandle:
        fleur_oxides = json.load(fhandle)
    with open('/fs/home/wolloch/aiida/acwf-verification-scripts/3-analyze/outputs/results-unaries-verification-PBE-v1-fleur.json') as fhandle:
        fleur_unaries = json.load(fhandle)
            
    noise_sigma = 1E-4
    nr_of_samples = 100

    stats, deviations = get_statistics(fleur_unaries, noise_sigma, nr_of_samples)
    plot_V0B0B1(stats, fleur_unaries, f'fleur unaries with {noise_sigma*1000} meV noise')
    
    stats, deviations = get_statistics(fleur_oxides, noise_sigma, nr_of_samples)
    plot_V0B0B1(stats, fleur_oxides, f'fleur oxides with {noise_sigma*1000} meV noise')
    