"""
This scripts reads the data of v3.1 of the "Science 2016 paper",
downloaded from https://molmod.ugent.be/deltacodesdft, and compares
with our new results, for the two AE codes.
The output is the % error (in %, so relative error multiplied by 100) for the discrepancy in the 3 Birch-Murnaghan parameters, for Wien2K and Fleur, old results (Science 2016) vs. new results (this work)
"""
import json
import numpy as np

overlapping_elements_str = """Ag FCC
Al FCC
Ar FCC
Au FCC
Ca FCC
Cu FCC
Ir FCC
Kr FCC
Mn FCC
Ne FCC
Ni FCC
Pb FCC
Pd FCC
Pt FCC
Rh FCC
Rn FCC
Sr FCC
Xe FCC
Ba BCC
Cs BCC
Fe BCC
K BCC
Mo BCC
Nb BCC
Rb BCC
Ta BCC
V BCC
W BCC
Cr BCC
Po SC
Ge Diamond
Si Diamond
Sn Diamond"""

overlapping_elements = {}
for l in overlapping_elements_str.splitlines():
    element, structure = l.split()
    overlapping_elements[element] = {'structure': structure}

with open('unaries-fleur.json') as f:
    raw = json.load(f)
    for element in overlapping_elements:
        structure = overlapping_elements[element]['structure']
        fit_data = raw["BM_fit_data"][f"{element}-X/{structure}"]
        formula_unit_atoms = 2 if structure == "Diamond" else 1
        overlapping_elements[element]['fleur'] = [
            fit_data['min_volume'] / formula_unit_atoms,
            fit_data['bulk_modulus_ev_ang3'] * 160.21766208, # in GPa
            fit_data['bulk_deriv']
        ]

with open('unaries-wien2k.json') as f:
    raw = json.load(f)
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

assert len(overlapping_elements) == 33
assert len(fleur_science) == 33
assert len(wien2k_science) == 33

#print(fleur_science)
#print(wien2k_science)

print(r"\begin{tabular}{ll|rrr|rrr}")
print(r"&& \multicolumn{3}{c|}{FLEUR} & \multicolumn{3}{c}{WIEN2k} \\")
print(r"Element & Structure & $V_0$ error [\%] & $B_0$ error [\%] & $B_1$ error [\%] & $V_0$ error [\%] & $B_0$ error [\%] & $B_1$ error [\%] \\ \hline")
for element in sorted(overlapping_elements):
    structure = overlapping_elements[element]['structure']
    wien2k_old = wien2k_science[element]
    fleur_old = fleur_science[element]
    wien2k_new = overlapping_elements[element]['wien2k']
    fleur_new = overlapping_elements[element]['fleur']

    wien2k_percent = 100 * np.abs((np.array(wien2k_old) - np.array(wien2k_new)) / np.array(wien2k_new))
    fleur_percent = 100 * np.abs((np.array(fleur_old) - np.array(fleur_new)) / np.array(fleur_new))

    print(fr"{element:2s} & {structure:10s} & {fleur_percent[0]:6.2f} & {fleur_percent[1]:5.2f} & {fleur_percent[2]:6.2f} & {wien2k_percent[0]:6.2f} & {wien2k_percent[1]:5.2f} & {wien2k_percent[2]:6.2f} \\ ")
print(r'\end{tabular}')