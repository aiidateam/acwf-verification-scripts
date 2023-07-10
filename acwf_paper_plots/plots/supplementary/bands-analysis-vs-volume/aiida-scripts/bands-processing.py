#!/usr/bin/env runaiida
import os
import numpy as np
from aiida.orm import load_node, load_group
import sys

try:
    element, configuration = sys.argv[1], sys.argv[2]
except IndexError:
    print("pass on the command line the element and the configuration")
    sys.exit(1)

GROUP = load_group('commonwf-oxides/set1/workflows/quantum_espresso')

#Group contains both relax and bands
WORKCHAINS = [n for n in GROUP.nodes if n.extras['element'] == element and n.extras['configuration'] == configuration] 
#print(WORKCHAINS)


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
for workchain in WORKCHAINS:
    structures = [workchain.outputs.structures[str(i)] for i in range(7)]
    pw_calcs_pk = [workchain.outputs.total_energies[str(i)].creator.inputs.parameters.creator.pk for i in range(7)]
    volumes = [get_volume(s.cell) for s in structures]
    volumes = list(zip(volumes, pw_calcs_pk))

    #print(structures)
    #print(volumes)

    # the composition is something like {'O': 1, 'Cs': 2}
    composition = structures[0].get_composition()

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
print(f"# {other_element}{other_count if other_count != 1 else ''}O{oxygen_count if oxygen_count != 1 else ''}")
latex_formula = f"{other_element}{get_count(other_count)}O{get_count(oxygen_count)}"
text_formula = f"{other_element}{get_count(other_count, ascii=True)}O{get_count(oxygen_count, ascii=True)}"

print() 

for idx, (volume, pk) in enumerate(volumes, start=1):
    pw_calc = load_node(pk)
    fermi_energy = pw_calc.res.fermi_energy # from the bands, should be the same as the SCF?
    print(f"{volume} {pw_calc.res.energy_smearing} # {pk}")

