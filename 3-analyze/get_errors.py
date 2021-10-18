#!/usr/bin/env runaiida
import os
import json

from aiida import orm
from aiida.cmdline.utils.common import get_workchain_report

SET_NAME = 'set2'

def get_plugin_name():
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

STRUCTURES_GROUP_LABEL = f'commonwf-oxides/{SET_NAME}/structures/{PLUGIN_NAME}'
WORKFLOWS_GROUP_LABEL = f'commonwf-oxides/{SET_NAME}/workflows/{PLUGIN_NAME}'

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

data = {}
for node in wf_nodes:
    structure = node.inputs.structure
    print(f"{structure.extras['element']} {structure.extras['configuration']} ({structure.pk}) -> {node.pk}: {node.process_state.value} ({node.exit_status})")
    data[f"{structure.extras['element']}-{structure.extras['configuration']}"] = {
        'structure': structure.uuid, 
        'eos_workflow': node.uuid, 
        'eos_workflow_report': get_workchain_report(node, levelname='REPORT'),
        'process_state': node.process_state.value,
        'exit_status': node.exit_status
    }

fname = f"outputs/errors-{PLUGIN_NAME}.json"
os.makedirs('outputs', exist_ok=True)
with open(fname, 'w') as fhandle:
    json.dump(data, fhandle, indent=2, sort_keys=True)
print(f"'{fname}' written.")
