#!/usr/bin/env python
import json

# Set2 for unaries is obtained from the first results of Wien2K and the fits of
# Fleur, taking the average. Just 3 cases are manually fixed according to 
# private messages.

def alat_from_volume(volume, config):
    if config in ['BCC']:
        return (volume * 2) ** (1/3)
    if config in ['FCC', 'Diamond']:
        return (volume * 4) ** (1/3)
    if config in ['SC']:
        return volume ** (1/3)
    raise ValueError(f"Unknown config {config}")

# Remember to move here the final files from FLEUR and WIEN2K!
with open("results-unaries-set2-fleur-k_0.06.json") as fhandle:
    data_wien2k = json.load(fhandle)

with open("results-unaries-set2-wien2k-dk_0.06.json") as fhandle:
    data_fleur = json.load(fhandle)

VALID_CONFIGS = ['BCC', 'FCC', 'SC', 'Diamond']

alats_fleur = {k: {} for k in VALID_CONFIGS}
for key in data_fleur["BM_fit_data"]:
    volume = data_fleur["BM_fit_data"][key]['min_volume']
    element = key.partition('-')[0]
    config = key.rpartition('/')[2]
    alats_fleur[config][element] = alat_from_volume(volume, config)

alats_wien2k = {k: {} for k in VALID_CONFIGS}
for key in data_wien2k["BM_fit_data"]:
    volume = data_wien2k["BM_fit_data"][key]['min_volume']
    element = key.partition('-')[0]
    config = key.rpartition('/')[2]
    alats_wien2k[config][element] = alat_from_volume(volume, config)

# Check that the configs and systems are the same - so then I can loop on only one of them
for config in VALID_CONFIGS:
    assert set(alats_fleur[config].keys()) == set(alats_wien2k[config].keys())

# This will be used as reference
with open("lattice_parameters_unaries_set2.json") as fhandle:
    alats_set2 = json.load(fhandle)

alats_set3 = {k: {} for k in VALID_CONFIGS}
for config in alats_set3.keys():
    for element in alats_fleur[config]:
        wien2k = alats_wien2k[config][element]
        fleur = alats_fleur[config][element]
        alats_set3[config][element] = (wien2k + fleur)/2

        rel_err = (fleur - wien2k) / wien2k
        if abs(rel_err) > 0.003:
            print(f"WARNING! Difference of {100*rel_err:.2f}% (in alat, {3*100*rel_err:.2f}% for V0) for {element}-X/{config}: {wien2k:.3f} [WIEN2K] vs {fleur:.3f} [FLEUR] ")
        
        rel_err_wrt_previous_set = (alats_set3[config][element] - alats_set2[config][element]) / alats_set2[config][element]
        if abs(rel_err_wrt_previous_set) > 0.003:
            print(f"WARNING! Difference WRT THE PREVIOUS SET2 of {100*rel_err_wrt_previous_set:.2f}% (in alat, {3*100*rel_err_wrt_previous_set:.2f}% for V0) for {element}-X/{config}: {alats_set2[config][element]:.3f} [SET2] vs {alats_set3[config][element]:.3f} [tentative SET3] ")
        
# I don't dump it, I just want to check the errors
print("All check done. If no other output was produced, then the set2 was already good and there is no strong need to generate a set 3!")


# As a reference: The only three warnings from the data available on March 15th 2022 are:
#
##WARNING! Difference WRT THE PREVIOUS SET2 of 0.30% (in alat, 0.90% for V0) for Yb-X/SC: 3.371 [SET2] vs 3.381 [tentative SET3] 
##WARNING! Difference WRT THE PREVIOUS SET2 of 0.40% (in alat, 1.19% for V0) for Cm-X/Diamond: 5.346 [SET2] vs 5.367 [tentative SET3] 
##WARNING! Difference WRT THE PREVIOUS SET2 of -0.38% (in alat, -1.14% for V0) for Th-X/Diamond: 7.180 [SET2] vs 7.152 [tentative SET3] 
