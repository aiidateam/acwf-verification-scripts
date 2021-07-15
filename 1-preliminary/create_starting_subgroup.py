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

STRUCTURES_FULL_GROUP_LABEL = 'commonwf-oxides/set1/structures'
STRUCTURES_GROUP_LABEL = f'commonwf-oxides/set1/structures/{PLUGIN_NAME}'

group = orm.Group.objects.get(label=STRUCTURES_FULL_GROUP_LABEL)
subgroup, _ = orm.Group.objects.get_or_create(label=STRUCTURES_GROUP_LABEL)

#####################################################################################
## Adapt here
query = orm.QueryBuilder()
query.append(orm.Node, project="attributes.element", tag='pseudo')
query.append(orm.Group, filters={'label': 'SSSP/1.1/PBE/precision'}, with_node='pseudo')
valid_elements = query.all(flat=True)
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

