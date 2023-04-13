#!/usr/bin/env runaiida
import os
import numpy as np
from aiida.orm import load_node, load_group
import sys

try:
    system = sys.argv[1]
except IndexError:
    print("pass on the command line the chemical formula")
    sys.exit(1)

GROUP = load_group(f'manual/bands/{system}/sssp')

#Group contains both relax and bands
IDS = [n.outputs.bands.pk for n in GROUP.nodes if n.process_label == 'QuantumEspressoCommonBandsWorkChain'] 

def get_volume(cell):
    return abs(np.dot(np.cross(cell[0], cell[1]), cell[2]))

def get_count(number, ascii=False):
    if number == 1:
        return ""
    else:
        if ascii:
            return str(number)
        else:
            return f"$_{number}$"

volumes = []
compositions = set()
for pk in IDS:
    bands = load_node(pk)
    volume = get_volume(bands.creator.inputs.structure.cell)
    volumes.append((volume, pk))
    # the composition is something like {'O': 1, 'Cs': 2}
    composition = bands.creator.inputs.structure.get_composition()

    assert len(composition) == 2 # Oxygen + something else
    oxygen_count = composition['O']
    del composition['O']
    other_element = list(composition.keys())[0]
    other_count = composition[other_element]

    compositions.add((oxygen_count, other_element, other_count))
volumes.sort()

# Get composition and write it in latex
assert len(compositions) == 1
oxygen_count, other_element, other_count = list(compositions)[0]
print(f"{other_element}{other_count if other_count != 1 else ''}O{oxygen_count if oxygen_count != 1 else ''}")
latex_formula = f"{other_element}{get_count(other_count)}O{get_count(oxygen_count)}"
text_formula = f"{other_element}{get_count(other_count, ascii=True)}O{get_count(oxygen_count, ascii=True)}"

print() 

OUT_FOLDER = f'out/{text_formula}/dataout'
os.makedirs(OUT_FOLDER, exist_ok=False)
print('# PK                                   Degauss (eV)     Fermi energy (eV)  volume (ang^3)')
for idx, (volume, pk) in enumerate(volumes, start=1):
    bands = load_node(pk)
    fermi_energy = bands.creator.res.fermi_energy # from the bands, should be the same as the SCF?
    print(f"# {pk} {bands.creator.res.degauss} {fermi_energy} {volume}")
    bands.export(
        f'{OUT_FOLDER}/{idx}-{pk}.py',
        fileformat='mpl_withjson',
        y_origin = fermi_energy, 
        plot_zero_axis=True, 
        y_min_lim=-3., 
        y_max_lim=+3., 
        title=f"{latex_formula} ($V = {volume:.2f}\\;\\mathrm{{\\AA}}^3$)", 
        overwrite=True
    )
