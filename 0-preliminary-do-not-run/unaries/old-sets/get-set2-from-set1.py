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

with open("lattice_parameters_unaries_set1-WIEN2K.json") as fhandle:
    alats_wien2k_old = json.load(fhandle)

# alats_wien2k fixes (communicated by Email)
# He-sc (V0=21.536 Ang)
# Dy-bcc (V0=32.2 ang**3)
alats_wien2k_old['bcc']['Dy'] = alat_from_volume(32.2, 'BCC')
alats_wien2k_old['sc']['He'] = alat_from_volume(21.536, 'SC')

# Use "new" labels
alats_wien2k = {
    'BCC': alats_wien2k_old['bcc'],
    'FCC': alats_wien2k_old['fcc'],
    'Diamond': alats_wien2k_old['dia'],
    'SC': alats_wien2k_old['sc'],
}

with open("results-unaries-set1-fleur_testPrecise_19-k_0.1.json") as fhandle:
    data_fleur = json.load(fhandle)

# data_fleur fixes (fit with parabola since Birch-Murnaghan failed (minimum out of range)
# I only need the minimum, non need to set the rest
data_fleur["BM_fit_data"]["He-X/SC"] = {
      "min_volume": 21.233,
    }

alats_fleur = {'BCC': {}, 'FCC': {}, 'SC': {}, 'Diamond': {}}
for key in data_fleur["BM_fit_data"]:
    volume = data_fleur["BM_fit_data"][key]['min_volume']
    element = key.partition('-')[0]
    config = key.rpartition('/')[2]

    alats_fleur[config][element] = alat_from_volume(volume, config)

# Check that the configs and systems are the same - so then I can loop on only one of them
assert set(alats_fleur.keys()) == set(alats_wien2k.keys())
for config in alats_fleur.keys():
    assert set(alats_fleur[config].keys()) == set(alats_wien2k[config].keys())


alats_set2 = {'BCC': {}, 'FCC': {}, 'SC': {}, 'Diamond': {}}
for config in alats_set2.keys():
    for element in alats_fleur[config]:
        wien2k = alats_wien2k[config][element]
        fleur = alats_fleur[config][element]
        alats_set2[config][element] = (wien2k + fleur)/2

        rel_err = (fleur - wien2k) / wien2k
        if abs(rel_err) > 0.005:
            print(f"WARNING! Difference of {100*rel_err:.2f}% (in alat, {3*100*rel_err:.2f}% for V0) for {element}-X/{config}: {wien2k:.3f} [WIEN2K] vs {fleur:.3f} [FLEUR] ")
        

with open('lattice_parameters_unaries_set2.json', 'w') as fhandle:
    json.dump(alats_set2, fhandle, indent=2, sort_keys=True)
print("File 'lattice_parameters_unaries_set2.json' written.")
