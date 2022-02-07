from aiida.plugins import DataFactory
from aiida import orm
from aiida.engine import calcfunction

SET_NAME = 'monoelemental-structures-test2'

@calcfunction
def get_oxygen_variant(oxygen_node, configuration):
    new_oxygen = oxygen_node.clone()
    new_oxygen.set_extra('Z', 8)
    new_oxygen.set_extra('element', 'O')
    new_oxygen.set_extra('configuration', configuration)

    return {'modified_oxygen': new_oxygen}

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

STRUCTURES_FULL_GROUP_LABEL = f'commonwf-oxides/{SET_NAME}'
STRUCTURES_GROUP_LABEL = f'commonwf-oxides/{SET_NAME}/structures/{PLUGIN_NAME}'

group = orm.Group.objects.get(label=STRUCTURES_FULL_GROUP_LABEL)
subgroup, _ = orm.Group.objects.get_or_create(label=STRUCTURES_GROUP_LABEL)

#####################################################################################
## PLUGIN-SPECIFIC PART: ADD THE ELIF FOR YOUR CODE
if PLUGIN_NAME == 'quantum_espresso':
    # First, get valid elements
    query = orm.QueryBuilder()
    query.append(orm.Node, project="attributes.element", tag='pseudo')
    query.append(orm.Group, filters={'label': 'SSSP/1.1/PBE/precision'}, with_node='pseudo')
    valid_elements = set(query.all(flat=True))
    print(f"Number of valid elements: {len(valid_elements)}")

    # Remove oxygen that will be treated separately
    valid_elements.remove("O")

    # Then, get the structures and put them into a new list
    query = orm.QueryBuilder()
    query.append(orm.Node, tag='structure', project=['*'], filters={
        'extras.element': {'in': list(valid_elements)}
    })
    query.append(orm.Group, tag='group', filters={'label': STRUCTURES_FULL_GROUP_LABEL}, with_node='structure')
    valid_structures = list(query.all(flat=True))

    # Finally, treat oxygen
    # For QE, I want to compute one for every possible combination of the two cutoffs
    # This table is obtained with the script utils/get_ecuts_table_qe.py (for set 1)

    # I hardcode here a table; NOTE: THIS IS ONLY VALID FOR THE USED PSEUDO LIBRARY AND PROTOCOL!
    # This was done for SSSP. In principle, it would be better to generate this table
    # by calling through the actual functions of the protocol.
    # the file in 2-submit/utils/get_ecuts_table_qe.py can still be used as a further check, if
    # you already run the oxides, to also assert that all calculatons of a given oxide family
    # with the same element used the same cutoffs.
    QE_ECUTS_TABLE = {
        'Ag': (75, 600),
        'Al': (75, 600),
        'Ar': (120, 600),
        'As': (75, 600),
        'Au': (75, 600),
        'B': (75, 600),
        'Ba': (75, 600),
        'Be': (75, 600),
        'Bi': (75, 600),
        'Br': (90, 720),
        'C': (75, 600),
        'Ca': (75, 600),
        'Cd': (90, 720),
        'Ce': (75, 600),
        'Cl': (100, 800),
        'Co': (90, 1080),
        'Cr': (75, 600),
        'Cs': (75, 600),
        'Cu': (90, 600),
        'Dy': (75, 600),
        'Er': (75, 600),
        'Eu': (75, 600),
        'F': (90, 600),
        'Fe': (90, 1080),
        'Ga': (90, 720),
        'Gd': (75, 600),
        'Ge': (75, 600),
        'H': (80, 600),
        'He': (75, 600),
        'Hf': (75, 600),
        'Hg': (75, 600),
        'Ho': (75, 600),
        'I': (75, 600),
        'In': (75, 600),
        'Ir': (75, 600),
        'K': (75, 600),
        'Kr': (75, 600),
        'La': (75, 600),
        'Li': (75, 600),
        'Lu': (75, 600),
        'Mg': (75, 600),
        'Mn': (90, 1080),
        'Mo': (75, 600),
        'N': (80, 600),
        'Na': (100, 600),
        'Nb': (75, 600),
        'Nd': (75, 600),
        'Ne': (75, 600),
        'Ni': (75, 600),
        'Os': (75, 600),
        'P': (75, 600),
        'Pb': (75, 600),
        'Pd': (75, 600),
        'Pm': (75, 600),
        'Po': (80, 640),
        'Pr': (75, 600),
        'Pt': (100, 800),
        'Rb': (75, 600),
        'Re': (75, 600),
        'Rh': (75, 600),
        'Rn': (200, 1600),
        'Ru': (75, 600),
        'S': (75, 600),
        'Sb': (75, 600),
        'Sc': (90, 720),
        'Se': (75, 600),
        'Si': (75, 600),
        'Sm': (75, 600),
        'Sn': (75, 600),
        'Sr': (75, 600),
        'Ta': (75, 600),
        'Tb': (75, 600),
        'Tc': (75, 600),
        'Te': (75, 600),
        'Ti': (75, 600),
        'Tl': (75, 600),
        'Tm': (75, 600),
        'V': (75, 600),
        'W': (75, 600),
        'Xe': (80, 600),
        'Y': (75, 600),
        'Yb': (75, 600),
        'Zn': (90, 720),
        'Zr': (75, 600)
    }


    possible_cutoffs = set(QE_ECUTS_TABLE.values())
    print(f"{len(possible_cutoffs)} possible cutoffs possible")
    oxygen_configurations = [f"{ecutwfc}-{ecutrho}" for ecutwfc, ecutrho in possible_cutoffs]

    # Get the oxygen
    query = orm.QueryBuilder()
    query.append(orm.Node, tag='structure', project=['*'], filters={
        'extras.element': 'O'
    })
    query.append(orm.Group, tag='group', filters={'label': STRUCTURES_FULL_GROUP_LABEL}, with_node='structure')
    # There should be only one result, I get it
    oxygen_node = list(query.all(flat=True))
    assert len(oxygen_node) == 1
    oxygen_node = oxygen_node[0]

    for configuration in oxygen_configurations:
        valid_structures.append(get_oxygen_variant(oxygen_node, orm.Str(configuration))['modified_oxygen'])
    
#elif PLUGIN_NAME == 'xxx':
#    yyy
else:
    raise ValueError(f"Unknown plugin name `{PLUGIN_NAME}`!")
#####################################################################################

# Now I add all valid structures to the group.
subgroup.add_nodes(valid_structures)

print(f"Structures from full group added to group '{STRUCTURES_GROUP_LABEL}'")
print(f"Current group size: {len(subgroup.nodes)}")

