#!/usr/bin/env runaiida
import sys

from aiida.plugins import DataFactory
from aiida import orm


def get_plugin_name():
    import os
    file_name = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        os.pardir, 'plugin_name.txt'
    )
    try:
        with open(file_name) as fhandle:
            plugin_name = fhandle.read().strip()
            # Simple check e.g. to make sure there are no weird characters,
            # newlines, ... - one might still make a typo, but at least we
            # do a basic check
            assert plugin_name.isidentifier()
        return plugin_name
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "You need to define a file `../plugin_name.txt`, containing the "
            "name of your plugin (siesta, quantum_espresso, ...) in the format "
            "expected by the aiida-common-workflows project"
        ) from exc

PLUGIN_NAME = get_plugin_name()

if __name__ == "__main__":

    try:
        SET_NAME = sys.argv[1]
    except IndexError:
        print("Pass as parameter the set name, e.g. oxides-verification-PBE-v1 or unaries-verification-PBE-v1")
        sys.exit(1)

    STRUCTURES_FULL_GROUP_LABEL = f'acwf-verification/{SET_NAME}/structures'
    STRUCTURES_GROUP_LABEL = f'acwf-verification/{SET_NAME}/structures/{PLUGIN_NAME}'

    group = orm.Group.collection.get(label=STRUCTURES_FULL_GROUP_LABEL)
    subgroup, _ = orm.Group.collection.get_or_create(label=STRUCTURES_GROUP_LABEL)

    #####################################################################################
    ## PLUGIN-SPECIFIC PART: ADD THE ELIF FOR YOUR CODE
    if PLUGIN_NAME == 'quantum_espresso':
        query = orm.QueryBuilder()
        query.append(orm.Node, project="attributes.element", tag='pseudo')
        query.append(orm.Group, filters={'label': 'SSSP/1.1/PBE/precision'}, with_node='pseudo')
        valid_elements = query.all(flat=True)
    elif PLUGIN_NAME == 'cp2k':
        valid_elements = [
                'H', 'He',
                'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
                'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
                'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
                'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
    #            'Fr',  'Ra', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og', 'Eu',
                ]
    #elif PLUGIN_NAME == 'xxx':
    #    yyy
    elif PLUGIN_NAME == 'gpaw':
        from ase.data import atomic_numbers
        query = orm.QueryBuilder()
        query.append(orm.Node, project="attributes.element", tag='pseudo')
        query.append(orm.Group, filters={'label': 'SSSP/1.1/PBE/precision'}, with_node='pseudo')
        valid_elements = [a for a in atomic_numbers.keys() if atomic_numbers[a] <= 83]
        for elements in ['Dy', 'Ce', 'Er', 'Eu', 'Gd', 'Ho', 'La', 'Lu', 'Nd', 'Pm', 'Pr', 'Sm', 'Tb', 'Tc', 'Tm', 'Yb' ]:
            valid_elements.remove(elements)
    else:
        raise ValueError(f"Unknown plugin name `{PLUGIN_NAME}`!")
    #####################################################################################

    print(f"Number of valid elements: {len(valid_elements)}")

    query = orm.QueryBuilder()
    query.append(orm.Node, tag='structure', project=['*'], filters={
        'extras.element': {'in': valid_elements}
    })
    query.append(orm.Group, tag='group', filters={'label': STRUCTURES_FULL_GROUP_LABEL}, with_node='structure')
    valid_structures = query.all(flat=True)
    subgroup.add_nodes(valid_structures)

    print(f"Structures from full group added to group '{STRUCTURES_GROUP_LABEL}'")
    print(f"Current group size: {len(subgroup.nodes)}")

