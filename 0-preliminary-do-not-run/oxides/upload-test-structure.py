"""
A simple script to upload a test structure and add it to a new test group, to allow for simple
testing with the same infrastructure. The file is loaded from a CIF from disk, and the 'set name'
is the basename of the file. You need to also pass the code name (to essentially jump step 1, and
directly go to step 2 and run with your code). You also need to pass the element name and the
configuration.
"""
import pymatgen
from ase.data import chemical_symbols, atomic_numbers

import sys
from aiida.plugins import DataFactory
from aiida import orm

set_name = sys.argv[1] # Should be also the basename of a CIF in this folder
code_name = sys.argv[2]
element_name = sys.argv[3]
configuration = sys.argv[4]

Z = atomic_numbers[element_name]

STRUCTURES_FULL_GROUP_LABEL = f'acwf-verification/{set_name}/structures/{code_name}'

Structure = DataFactory('structure')
group, created = orm.Group.objects.get_or_create(label=STRUCTURES_FULL_GROUP_LABEL)

if not created and len(group.nodes) > 0:
    print(f"Stopping, non-empty group '{STRUCTURES_FULL_GROUP_LABEL}' already exists")
    sys.exit(1)

filename = f'{set_name}.cif'
pmg_structure = pymatgen.Structure.from_file(filename)
structure = Structure(pymatgen=pmg_structure)
structure.set_extra_many({'element': element_name, 'Z': Z, 'configuration': configuration})
structure.store()
group.add_nodes([structure])
print(f"Done, group {STRUCTURES_FULL_GROUP_LABEL} created.")
