
from aiida.plugins import DataFactory, WorkflowFactory
from aiida import orm
from aiida.engine import submit

from aiida_common_workflows.common import ElectronicType, RelaxType, SpinType
from aiida_common_workflows.plugins import get_entry_point_name_from_class
from aiida_common_workflows.plugins import load_workflow_entry_point

PLUGIN_NAME = 'quantum_espresso'

STRUCTURES_GROUP_LABEL = f'commonwf-oxides/set1/structures/{PLUGIN_NAME}'
WORKFLOWS_GROUP_LABEL = f'commonwf-oxides/set1/workflows/{PLUGIN_NAME}'
CODE_LABEL = 'qe-6.7-pw@daint-mc'

Structure = DataFactory('structure')

query = orm.QueryBuilder()
query.append(Structure, tag='structure', project=['extras', '*'])
query.append(orm.Group, tag='group', filters={'label': STRUCTURES_GROUP_LABEL}, with_node='structure')
all_structures = {(res[0]['element'], res[0]['configuration']): res[1] for res in query.all()}

query = orm.QueryBuilder()
query.append(Structure, tag='structure', project='extras')
query.append(orm.Group, tag='group', filters={'label': STRUCTURES_GROUP_LABEL}, with_node='structure')
query.append(orm.WorkChainNode, tag='wc', project='id', with_incoming='structure')
query.append(orm.Group, tag='wcgroup', filters={'label': WORKFLOWS_GROUP_LABEL}, with_node='wc')
already_run_wc = {(res[0]['element'], res[0]['configuration']): res[1] for res in query.all()}
systems_keys_to_run = sorted(set(all_structures).difference(already_run_wc))

print(f'{len(already_run_wc)} run out of {len(all_structures)}.')    
print()
structure = all_structures[('Ag', 'X2O')]
print(f'Structure PK: {structure.pk}')

# For reference, from the command line
# aiida-common-workflows launch eos quantum_espresso -p precise -S 428 -X qe-6.7-pw@daint-mc -r none -s collinear -m 1 -w 3600 -d --engine-options='{"relax": {"account": "mr0"}}'

sub_process_cls = load_workflow_entry_point('relax', 'quantum_espresso')
sub_process_cls_name = get_entry_point_name_from_class(sub_process_cls).name
generator = sub_process_cls.get_input_generator()

engine_types = generator.get_engine_types()
engines = {}
# There should be only one
for engine in engine_types:
    engines[engine] = {
        'code': CODE_LABEL,
         'options': {
            'resources': {
                'num_machines': 1
            },
            'account': 'mr0',
            'max_wallclock_seconds': 3600
        }
    }

inputs = {
    'structure': structure,
    'generator_inputs': {  # code-agnostic inputs for the relaxation
        'engines': engines,
        'protocol': 'precise',
        'relax_type': RelaxType.NONE,
        'electronic_type': ElectronicType.METAL,
        'spin_type': SpinType.NONE,
    },
    'sub_process_class': sub_process_cls_name,
    'sub_process' : {  # optional code-dependent overrides
        'base': {
            'pw': {
                'settings' : orm.Dict(dict= {
                    'cmdline': ['-nk', '18'],
                })
            }
        }
    }
}

cls = WorkflowFactory('common_workflows.eos')
submit(cls, **inputs)
