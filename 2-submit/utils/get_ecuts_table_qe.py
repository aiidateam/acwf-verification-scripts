#!/usr/bin/env runaiida
from collections import defaultdict
import os
import json

from aiida import orm
from aiida.cmdline.utils.common import get_workchain_report

SET_NAME = 'set1'

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

STRUCTURES_GROUP_LABEL = f'acwf-verification/{SET_NAME}/structures/{PLUGIN_NAME}'
WORKFLOWS_GROUP_LABEL = f'acwf-verification/{SET_NAME}/workflows/{PLUGIN_NAME}'

query = orm.QueryBuilder()
query.append(orm.Group, filters={'label': WORKFLOWS_GROUP_LABEL}, tag='groups')
query.append(orm.Node, with_group='groups', tag='eoswf', project='*')
query.append(orm.StructureData, with_outgoing='eoswf', project='extras')
query.append(orm.Float, with_incoming='eoswf', edge_filters={'label': {'like': r"total\_energies\_\_%"}}, tag='energy')
query.append(orm.CalcFunctionNode, with_outgoing='energy', tag='energycf')
query.append(orm.Dict, with_outgoing='energycf', tag='outparams')
query.append(orm.CalcJobNode, with_outgoing='outparams', tag='pw')
query.append(orm.Dict, with_outgoing='pw', tag='pwinput', edge_filters={'label': 'parameters'}, project='*')
query.distinct()
results = query.all()


data = defaultdict(set)
for eos_wf, structure_extras, params in results:
    data[structure_extras['element']].add((params['SYSTEM']['ecutwfc'], params['SYSTEM']['ecutrho']))

final_data = {}
for key, value in data.items():
    assert len(value) == 1, f"Problem for {key}: {value}"
    final_data[key] = sorted(value)[0]
print(f"### GETTING FROM SET '{SET_NAME}' ###")
print(final_data)