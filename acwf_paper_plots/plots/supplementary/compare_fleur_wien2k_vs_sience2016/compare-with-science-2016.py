#!/usr/bin/env python
"""
This scripts reads the data of v3.1 of the "Science 2016 paper",
downloaded from https://molmod.ugent.be/deltacodesdft, and compares
with our new results, for the two AE codes.
The output is the % error (in %, so relative error multiplied by 100) for the discrepancy in the 3 Birch-Murnaghan parameters, for Wien2K and Fleur, old results (Science 2016) vs. new results (this work)
"""
import json
import os
import numpy as np
import pylab as pl
import acwf_paper_plots.quantities_for_comparison as qc


# The following list summarizes the subset of the 71 elements from the
# Science 2016 paper that have a FCC, BCC, SC or Diamond structure
# and thus overlap with the set of unaries of the current work.

# Note that we are removing Cr, Mn, Fe and Ni because they are treated
# using magnetization in the Science 2016 paper:
# > Include spin polarization for O, Cr and Mn (antiferromagnetic)
# > and Fe, Co, and Ni (ferromagnetic). 
# > The antiferromagnetism in the oxygen crystal should take place
# > between the O2 molecules.
# (O and Co are not in the list since they are not cubic structures).

# For reference, the string for these 4 elements is reported below, 
# commented:
# 
## Fe BCC
## Mn FCC
## Ni FCC
## Cr BCC

overlapping_elements_str = """Ag FCC
Al FCC
Ar FCC
Au FCC
Ca FCC
Cu FCC
Ir FCC
Kr FCC
Ne FCC
Pb FCC
Pd FCC
Pt FCC
Rh FCC
Rn FCC
Sr FCC
Xe FCC
Ba BCC
Cs BCC
K BCC
Mo BCC
Nb BCC
Rb BCC
Ta BCC
V BCC
W BCC
Po SC
Ge Diamond
Si Diamond
Sn Diamond"""

overlapping_elements = {}
for l in overlapping_elements_str.splitlines():
    element, structure = l.split()
    overlapping_elements[element] = {'structure': structure}

DATA_FOLDER = "../../../code-data"
with open(os.path.join(DATA_FOLDER, "labels.json")) as fhandle:
    labels_data = json.load(fhandle)
FLEUR_LABEL = labels_data['all-electron-keys']["FLEUR"]
WIEN2k_LABEL = labels_data['all-electron-keys']["WIEN2k"]

with open(os.path.join(DATA_FOLDER, labels_data['methods-main'][FLEUR_LABEL]["unaries"])) as fhandle:             
    raw = json.load(fhandle)
    for element in overlapping_elements:
        structure = overlapping_elements[element]['structure']
        fit_data = raw["BM_fit_data"][f"{element}-X/{structure}"]
        formula_unit_atoms = 2 if structure == "Diamond" else 1
        overlapping_elements[element]['fleur'] = [
            fit_data['min_volume'] / formula_unit_atoms,
            fit_data['bulk_modulus_ev_ang3'] * 160.21766208, # in GPa
            fit_data['bulk_deriv']
        ]

with open(os.path.join(DATA_FOLDER, labels_data['methods-main'][WIEN2k_LABEL]["unaries"])) as fhandle:             
    raw = json.load(fhandle)
    for element in overlapping_elements:
        structure = overlapping_elements[element]['structure']
        fit_data = raw["BM_fit_data"][f"{element}-X/{structure}"]
        formula_unit_atoms = 2 if structure == "Diamond" else 1
        overlapping_elements[element]['wien2k'] = [
            fit_data['min_volume'] / formula_unit_atoms,
            fit_data['bulk_modulus_ev_ang3'] * 160.21766208, # in GPa
            fit_data['bulk_deriv']
        ]

#print(overlapping_elements)

# READ data from Science paper
fleur_science = {}
with open('FLEUR-history.txt') as f:
    for l in f:
        if not l.strip() or l.startswith("#"):
            continue
        element, V0, B0, B1 = l.split()
        if element in overlapping_elements:
            fleur_science[element] = [float(V0), float(B0), float(B1)]

wien2k_science = {}
with open('WIEN2k-history.txt') as f:
    for l in f:
        if not l.strip() or l.startswith("#"):
            continue
        element, V0, B0, B1 = l.split()
        if element in overlapping_elements:
            wien2k_science[element] = [float(V0), float(B0), float(B1)]

# 33 overlapping, minus 4 magnetic == 29
assert len(overlapping_elements) == 29
assert len(fleur_science) == 29
assert len(wien2k_science) == 29

save_delta_new = []
save_delta_old = []
save_conf = []
for element in sorted(overlapping_elements):
    structure = overlapping_elements[element]['structure']
    wien2k_old = wien2k_science[element]
    fleur_old = fleur_science[element]
    save_delta_old.append(qc.delta(wien2k_old[0], wien2k_old[1]/160.21766208, wien2k_old[2], fleur_old[0], fleur_old[1]/160.21766208, fleur_old[2], 1, 0, 0))
    wien2k_new = overlapping_elements[element]['wien2k']
    fleur_new = overlapping_elements[element]['fleur']
    save_delta_new.append(qc.delta(wien2k_new[0], wien2k_new[1]/160.21766208, wien2k_new[2], fleur_new[0], fleur_new[1]/160.21766208, fleur_new[2], 1, 0, 0))

    wien2k_percent = 100 * np.abs((np.array(wien2k_old) - np.array(wien2k_new)) / np.array(wien2k_new))
    fleur_percent = 100 * np.abs((np.array(fleur_old) - np.array(fleur_new)) / np.array(fleur_new))

    save_conf.append(f"{element}-{structure}")

print('Delta_old, Delta_new')
correct_new = []
correct_old = []
for conf, new, old in zip(save_conf, save_delta_new, save_delta_old):
    correct_new.append(new)
    correct_old.append(old)
    print(conf, old, new)

xlabel = r"$\Delta$ of Science 2016"
ylabel = r"$\Delta$ of this work"
filename = "delta-new-vs-science-2016.pdf"
fig = pl.figure(figsize=(5,5))
#print(xdata)
#print(ydata)
pl.plot(correct_old, correct_new, 'o')
for i in range(len(correct_new)):
    if save_conf[i].split('-')[0] == 'Cu':
        pl.annotate(save_conf[i].split('-')[0], (correct_old[i]+0.01, correct_new[i]-0.01), fontsize=8)
    elif save_conf[i].split('-')[0] == 'Cs':
        pl.annotate(save_conf[i].split('-')[0], (correct_old[i]-0.01, correct_new[i]-0.05), fontsize=8)
    elif save_conf[i].split('-')[0] == 'K':
        pl.annotate(save_conf[i].split('-')[0], (-0.07, -0.03), fontsize=8)
        pl.plot([correct_old[i], -0.04], [correct_new[i], -0.02], '--k', linewidth=0.6)
    elif save_conf[i].split('-')[0] == 'Xe':
        pl.annotate(save_conf[i].split('-')[0], (-0.049, -0.09), fontsize=8)
        pl.plot([correct_old[i], -0.04], [correct_new[i], -0.06], '--k', linewidth=0.6)
    elif save_conf[i].split('-')[0] == 'Rn':
        pl.annotate(save_conf[i].split('-')[0], (+0.018, -0.05), fontsize=8)
        pl.plot([correct_old[i], +0.039], [correct_new[i], -0.02], '--k', linewidth=0.6)
    elif save_conf[i].split('-')[0] == 'Kr':
        pl.annotate(save_conf[i].split('-')[0], (+0.06, -0.09), fontsize=8)
        pl.plot([correct_old[i], +0.069], [correct_new[i], -0.06], '--k', linewidth=0.6)
    elif save_conf[i].split('-')[0] == 'Pb':
        pl.annotate(save_conf[i].split('-')[0], (+0.08, -0.05), fontsize=8)
        pl.plot([correct_old[i], +0.1], [correct_new[i], -0.02], '--k', linewidth=0.6)
    elif save_conf[i].split('-')[0] == 'Rb':
        pl.annotate(save_conf[i].split('-')[0], (+0.12, -0.09), fontsize=8)
        pl.plot([correct_old[i], +0.15], [correct_new[i], -0.06], '--k', linewidth=0.6)
    elif save_conf[i].split('-')[0] == 'Ba':
        pl.annotate(save_conf[i].split('-')[0], (+0.18, -0.09), fontsize=8)
        pl.plot([correct_old[i], +0.186], [correct_new[i], -0.085], '--k', linewidth=0.6)
    elif save_conf[i].split('-')[0] == 'Po':
        pl.annotate(save_conf[i].split('-')[0], (correct_old[i]+0.02, correct_new[i]-0.02), fontsize=8)
    elif save_conf[i].split('-')[0] == 'Ge':
        pl.annotate(save_conf[i].split('-')[0], (correct_old[i]+0.005, correct_new[i]+0.015), fontsize=8)
    elif save_conf[i].split('-')[0] in['Sr','Ag']:
        pl.annotate(save_conf[i].split('-')[0], (correct_old[i]-0.05, correct_new[i]+0.01), fontsize=8)
    elif save_conf[i].split('-')[0] == 'Si':
        pl.annotate(save_conf[i].split('-')[0], (correct_old[i]-0.01, correct_new[i]+0.015), fontsize=8)
    elif save_conf[i].split('-')[0] == 'Ca':
        pl.annotate(save_conf[i].split('-')[0], (correct_old[i]-0.05, correct_new[i]+0.01), fontsize=8)
    elif save_conf[i].split('-')[0] == 'Sn':
        pl.annotate(save_conf[i].split('-')[0], (correct_old[i]-0.01, correct_new[i]+0.02), fontsize=8)
    elif save_conf[i].split('-')[0] in ['Rb', 'Ba']:
        pl.annotate(save_conf[i].split('-')[0], (correct_old[i]+0.5, correct_new[i]+0.5), fontsize=8)
    elif save_conf[i].split('-')[0] == 'Ar':
        pl.annotate(save_conf[i].split('-')[0], (-0.07, 0.03), fontsize=8)
        pl.plot([correct_old[i], -0.04], [correct_new[i], 0.035], '--k', linewidth=0.6)
    elif correct_new[i]< 0.01:
        pl.annotate(save_conf[i].split('-')[0], (correct_old[i]+0.01, correct_new[i]-0.03), fontsize=8)
    else:
        pl.annotate(save_conf[i].split('-')[0], (correct_old[i]+0.01, correct_new[i]+0.01), fontsize=8)
pl.xlim((-0.1,1.3))
pl.ylim((-0.1,1.3))
pl.plot([-0.1, 1.3], [-0.1, 1.3], '-k', linewidth=1)
pl.xlabel(xlabel)
pl.ylabel(ylabel)
#adjust_text(texts)

pl.savefig(filename)
