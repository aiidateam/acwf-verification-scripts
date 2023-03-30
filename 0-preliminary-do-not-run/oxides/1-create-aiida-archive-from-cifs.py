import pymatgen
from ase.data import chemical_symbols

from aiida.plugins import DataFactory
from aiida import orm

SET_NAME = 'oxides-verification-PBE-v1'

STRUCTURES_FULL_GROUP_LABEL = f'acwf-verification/{SET_NAME}/structures'

Structure = DataFactory('structure')
group, _ = orm.Group.objects.get_or_create(label=STRUCTURES_FULL_GROUP_LABEL)

query = orm.QueryBuilder()
query.append(Structure, tag='structure', project=['extras', 'id'])
query.append(orm.Group, tag='group', filters={'label': STRUCTURES_FULL_GROUP_LABEL}, with_node='structure')
all_structures = {(res[0]['element'], res[0]['configuration']): res[1] for res in query.all()}

structures = []
for Z in range(1, 96+1):
    element_name = chemical_symbols[Z]
    #if Z == 8: # Oxygen is not there
    #    continue
    for configuration in ['X2O', 'X2O5', 'XO2', 'X2O3', 'XO', 'XO3']:
        structure_pk = all_structures.get((element_name, configuration))
        if structure_pk is not None:
            print(f"- Skipping import of {element_name}-{configuration} as node {structure_pk} exists already in the group")
            continue
        filename = f'cifs-{SET_NAME}/POSCAR_{element_name}_{configuration}.cif'
        pmg_structure = pymatgen.core.Structure.from_file(filename)
        structure = Structure(pymatgen=pmg_structure)

        if SET_NAME != "set1":
            # For set1, we did not 'standardize' the structure.
            # In order to make sure the script generates the same 'set1',
            # we add this logic - but for any other set, we should
            # apply standardization. This is important so the cell vectors
            # are appropriately oriented w.r.t. the Cartesian axes.
            # For some codes (at least Quantum ESPRESSO), this will allow
            # proper detection of all symmetries (otherwise, e.g. for X2O
            # oxides, QE might detect only 4 out of 48 symmetries)
            from aiida.tools.data.array.kpoints import get_kpoints_path
            res = get_kpoints_path(structure, method='seekpath')
            structure = res['primitive_structure']

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
