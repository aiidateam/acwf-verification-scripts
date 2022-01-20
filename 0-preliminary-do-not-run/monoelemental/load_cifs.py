import pymatgen
from ase.data import chemical_symbols

from aiida.plugins import DataFactory
from aiida import orm

### ADAPT IF NEEDED - FOR NOW, TESTING DATA
# For now I use a dict; in the future maybe I should get it from the name directly?
CIF_FOLDER = 'cifs-test'
CIFS_DICT = {
    'Au': 'Au_mp-81_primitive.cif',
    'Ba': 'Ba_mp-122_primitive.cif',
    'Si': 'Si_mp-149_primitive.cif'
}
O_CIF = 'O2_mp-12957_primitive.cif'

STRUCTURES_FULL_GROUP_LABEL = 'commonwf-oxides/monoelemental-structures-test'
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
    try:
        filename = f'{CIF_FOLDER}/{CIFS_DICT[element_name]}'
    except KeyError:
        # For now I skip elements that are not in the CIFS_DICT
        continue
    pmg_structure = pymatgen.Structure.from_file(filename)
    structure = Structure(pymatgen=pmg_structure)
    structure.set_extra_many({'element': element_name, 'Z': Z, 'configuration': configuration})
    structures.append(structure)

    ## I also create one oxygen per other element - they are all the same, but with different extras so
    ## you can compute each of them with the right parameters (e.g. cutoffs).
    ## Possibly to be adapted in the future, see discussion in the script to submit monoelemental cyrstals,
    ## in folder 2-submit
    pmg_structure = pymatgen.Structure.from_file(filename = f'{CIF_FOLDER}/{O_CIF}')
    structure = Structure(pymatgen=pmg_structure)
    # Oxygen has Z=8; for the configuration of oxygen, I indicate for each of them the *other* element
    structure.set_extra_many({'element': 'O', 'Z': 8, 'configuration': element_name})
    structures.append(structure)


if structures:
    print(f'{len(structures)} structures loaded from CIF and ready to add to group, now storing them in group...')
    for structure in structures:
        structure.store()
    group.add_nodes(structures)
    print("Done.")
else:
    print(f'No structure to add to the group: job done.')
