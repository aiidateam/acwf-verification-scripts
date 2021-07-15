import os

from aiida import orm
from aiida.common import LinkType
from aiida.common import LinkType
from aiida.tools.graph.graph_traversers import traverse_graph

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

STRUCTURES_GROUP_LABEL = f'commonwf-oxides/set1/structures/{PLUGIN_NAME}'
WORKFLOWS_GROUP_LABEL = f'commonwf-oxides/set1/workflows/{PLUGIN_NAME}'

group_node_query = orm.QueryBuilder().append(
    orm.Group, filters={'label': WORKFLOWS_GROUP_LABEL}, tag='groups',
).append(orm.Node, project='*', with_group='groups', filters={
    'and': [
        {'attributes.exit_status': {"!==": 0}},
        {'attributes.process_state': {"!==": 'waiting'}},
    ]
})
group_node_query.distinct()
wf_nodes = group_node_query.all(flat=True)

os.makedirs('outputs', exist_ok=True)
fname = f"outputs/errors-{PLUGIN_NAME}.txt"
with open(fname, 'w') as fhandle:
    for node in wf_nodes:
        structure = node.inputs.structure
        fhandle.write(f"{structure.extras['element']} {structure.extras['configuration']} ({structure.pk}) -> {node.pk}: {node.process_state.value} ({node.exit_status})\n")
print(f"'{fname}' written.")
