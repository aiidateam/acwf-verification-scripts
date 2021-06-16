import os

from aiida import orm
from aiida.common import LinkType
from aiida.common import LinkType
from aiida.tools.graph.graph_traversers import traverse_graph

PLUGIN_NAME = 'quantum_espresso'

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
