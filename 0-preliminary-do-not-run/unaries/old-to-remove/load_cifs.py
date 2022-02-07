import os

import pymatgen
from ase.data import chemical_symbols

from aiida.plugins import DataFactory
from aiida import orm

### ADAPT IF NEEDED - FOR NOW, TESTING DATA
# For now I use a dict; in the future maybe I should get it from the name directly?
CIF_FOLDER = 'deltacodesdft-primCIFs-test2'
STRUCTURES_FULL_GROUP_LABEL = 'commonwf-oxides/monoelemental-structures-test2'
#####

Structure = DataFactory('structure')
group, _ = orm.Group.objects.get_or_create(label=STRUCTURES_FULL_GROUP_LABEL)

query = orm.QueryBuilder()
query.append(Structure, tag='structure', project=['extras', 'id'])
query.append(orm.Group, tag='group', filters={'label': STRUCTURES_FULL_GROUP_LABEL}, with_node='structure')
all_structures = {(res[0]['element'], res[0]['configuration']): res[1] for res in query.all()}

configuration = 'X' # constant (even O2 is considered X, I don't discuss the supercell size here)

structures = []
for Z in range(1, 96+1):
    element_name = chemical_symbols[Z]
    structure_pk = all_structures.get((element_name, configuration))
    if structure_pk is not None:
        print(f"- Skipping import of {element_name}-{configuration} as node {structure_pk} exists already in the group")
        continue

    filename = f'{CIF_FOLDER}/{element_name}.cif'
    if not os.path.exists(filename):
        # Skip missing elements (not all are there)
        continue

    pmg_structure = pymatgen.Structure.from_file(filename)
    structure = Structure(pymatgen=pmg_structure)
    structure.set_extra_many({'element': element_name, 'Z': Z, 'configuration': configuration})
    structures.append(structure)

if structures:
    print(f'{len(structures)} structures loaded from CIF and ready to add to group, now storing them in group...')
    for structure in structures:
        structure.store()
    group.add_nodes(structures)
    print("Done.")
else:
    print(f'No structure to add to the group: job done.')
