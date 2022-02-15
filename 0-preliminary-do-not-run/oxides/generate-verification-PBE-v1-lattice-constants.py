"""
Uses results from WIEN2K and FLEUR to create new structures where
the volume is the avarage of the equilibrium volume for FLEUR and WIEN2K.
Also adds the Oxigen oxides to the set.
Produces "verification-PBE-v1_alat_and_vol.json"
"""

import json
import os
import sys

import numpy as np
import pylab as pl


#First part, using set2 results

with open('all_el_res_set2/results-wien2k-finite-T_delta-k_0.1_prec3.json') as fhandle:
    WIEN2K_data = json.load(fhandle)

with open('all_el_res_set2/results-fleur_testPrecise_16.json') as f2handle:
    FLEUR_data = json.load(f2handle)

all_systems = set(FLEUR_data['BM_fit_data'].keys())

collect={}

for syst in sorted(all_systems):
    wien_v = FLEUR_data["BM_fit_data"][syst]["min_volume"]
    fleur_v = WIEN2K_data["BM_fit_data"][syst]["min_volume"]
    av_volume = (wien_v+fleur_v)/2
    rel_diff = 100*abs(wien_v-fleur_v)/av_volume
    #if rel_diff>0.2:
    #    print(syst, rel_diff)
    element, configuration = syst.split('-')
    if configuration in ["X2O","XO2","XO"]:
        a0 = a0 = (av_volume*np.sqrt(2))**(1.0/3.0)
    else:
        a0 = av_volume**(1.0/3.0)

    collect[syst] = {
            "vol_rel_diff_wien2k_fleur": rel_diff,
            "average_vol_wien2k_fleur": av_volume,
            "new_latt_constant": a0
            }


#Second part, using set2-bis results.
#The set2-bis is a set containing the oxides for which the equilibrim volume
#was out of range and the rel error between FLEUR and WIEN2K was > 0.2%
#It also adds the oxigen oxides.

with open('all_el_res_set2bis/results-wien2k.set2-bis.dk-0.1invAng.json') as f4handle:
    new_WIEN2K_data = json.load(f4handle)

with open('all_el_res_set2bis/results-set2-bis-fleur_testPrecise_18_k_0.1.json') as f5handle:
    new_FLEUR_data = json.load(f5handle)

new_all_systems = set(new_FLEUR_data['BM_fit_data'].keys())

print("        Comparing set2 and se2-bis")
print("name    set2_av_vol  set2-bis_av_vol   set2_vol_err     set2-bis_vol_err ")
for syst in sorted(new_all_systems):
    wien_v = new_FLEUR_data["BM_fit_data"][syst]["min_volume"]
    fleur_v = new_WIEN2K_data["BM_fit_data"][syst]["min_volume"]
    av_volume = (wien_v+fleur_v)/2
    rel_diff = 100*abs(wien_v-fleur_v)/av_volume
    #if rel_diff>0.2:
    #    print(syst, rel_diff)
    element, configuration = syst.split('-')
    if configuration in ["X2O","XO2","XO"]:
        a0 = (av_volume*np.sqrt(2))**(1.0/3.0)
    else:
        a0 = av_volume**(1.0/3.0)

    try:
        collect[syst]
        print(syst, collect[syst]["average_vol_wien2k_fleur"], av_volume, collect[syst]["vol_rel_diff_wien2k_fleur"], rel_diff)
    except:
        pass

    collect[syst] = {
            "vol_rel_diff_wien2k_fleur": rel_diff,
            "average_vol_wien2k_fleur": av_volume,
            "new_latt_constant": a0
            }


with open("verification-PBE-v1_alat_and_vol.json", 'w') as f3handle:
    json.dump(collect, f3handle, indent=2, sort_keys=True)
