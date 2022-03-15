#!/usr/bin/env runaiida
import os
import json
import sys

import numpy as np
import ase.build
import ase.data
import ase.visualize

from aiida.orm import StructureData, Group, QueryBuilder

# Cells, in units of alat
SC_CELL = [1.0, 1.0, 1.0]
BCC_CELL = [[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]
FCC_CELL = [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
DIAMOND_CELL = FCC_CELL

Z_MAX = 96
SET_NAME = 'unaries-verification-PBE-v1'
AIIDA_GROUP_LABEL = f'acwf-verification/{SET_NAME}/structures'
GENERATE_XSF_FILES = False

def get_cubic_unaries(element, alats):
    Z = ase.data.atomic_numbers[element]
    # From 1 to 96 they are all indicated (Z>96 they are set to 0.2, probably a placeholder)
    covalent_radius = ase.data.covalent_radii[Z]
    structures = {}

    # First-neighbor distance: the distance (in units of alat,
    # the cubic lattice parameter) of the first-neighbor distance
    # However, I don't use it here since I will use the a_lat directly
    for crystal_structure_label, cell in [
            ('Diamond', DIAMOND_CELL), # Coordination 4
            ('SC', SC_CELL), # Coordination 6
            ('BCC', BCC_CELL), # Coordination 8
            ('FCC', FCC_CELL), # Coordination 12
        ]:
        alat = alats[crystal_structure_label][element]

        num_atoms = 1
        if crystal_structure_label == "Diamond":
            num_atoms = 2
        atoms = ase.Atoms(
            symbols=f'{element}{"2" if num_atoms == 2 else ""}',
            pbc=True,
            cell=np.array(cell) * alat
        )
        if crystal_structure_label == "Diamond":
            atoms.set_scaled_positions(
                [
                    [0., 0., 0.],
                    [0.25, 0.25, 0.25],
                ])
        else:    
            # Probably the default is to set the atoms positions
            # to zero, but I do it again to be sure
            atoms.set_scaled_positions([[0., 0., 0.]])

        # I don't use this, for some elements it raises an exception?
        #atoms = ase.build.bulk(element, crystal_structure, a=1.)
        structures[crystal_structure_label] = atoms
    return structures

def get_aiida_structures(alats):
    all_structures = {}
    for Z in range(1, 96+1):
        element = ase.data.chemical_symbols[Z]
        unaries = get_cubic_unaries(element, alats)

        for crystal_structure_label, unary in unaries.items():
            structure = StructureData(ase=unary)
            configuration = f'X/{crystal_structure_label}'
            structure.set_extra('element', element)
            structure.set_extra('Z', Z)
            structure.set_extra('configuration', configuration)
            all_structures[(element, configuration)] = structure
    return all_structures



if __name__ == "__main__":
    with open("lattice_parameters_unaries-verification-PBE-v1.json") as fhandle:
        alats_wien2k = json.load(fhandle)

    subfolder = f'xsfs-{SET_NAME}'

    # Create XSF files
    if GENERATE_XSF_FILES:
        if os.path.exists(subfolder):
            print(f"Remove the '{subfolder}' folder first.")
            sys.exit(1)
        os.makedirs(subfolder)
        for Z in range(1, Z_MAX+1):
            element = ase.data.chemical_symbols[Z]
            unaries = get_cubic_unaries(element, alats_wien2k)
            for crystal_structure_label, unary in unaries.items():
                # Careful when storing as CIF, ASE will use low precision for the cell angles,
                # problematic for BCC!
                unary.write(f'{subfolder}/{element}-{crystal_structure_label}.xsf')
        print(f"{Z_MAX} structures written to the '{subfolder}' folder.")

    if GENERATE_XSF_FILES:
        print("WARNING: NOT WORKING ON AIIDA GROUP, AS GENERATE_XSF_FILES IS TRUE")
        sys.exit(0)

    # Create the group with structures
    all_structures = get_aiida_structures(alats_wien2k)

    qb = QueryBuilder()
    qb.append(Group, filters={'label': AIIDA_GROUP_LABEL}, tag='group')
    qb.append(StructureData, with_group='group', project=['extras.element', 'extras.configuration'])
    # Create a set of tuples
    existing = [(_[0], _[1]) for _ in qb.all()]
    # Check missing structures
    missing = set(all_structures.keys()).difference(existing)
    print(f"# {len(all_structures)} structures that should be in the group '{AIIDA_GROUP_LABEL}'")
    print(f"# {len(existing)} structures already in the group")
    print(f"# {len(missing)} structures that I will add to the group")

    group, _ = Group.objects.get_or_create(label=AIIDA_GROUP_LABEL)
    structures_to_add = [all_structures[extras] for extras in missing]
    # I need to store all these structures first
    for structure in structures_to_add:
        structure.store()
    # Now I add them to the group
    group.add_nodes(structures_to_add)
    if structures_to_add:
        print(f"# Nodes added to the group; current number of elements in the group: {len(group.nodes)}")
