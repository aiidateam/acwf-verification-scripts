from aiida.plugins import DataFactory, WorkflowFactory
from aiida import orm
from aiida.engine import submit

from aiida_common_workflows.common import ElectronicType, RelaxType, SpinType
from aiida_common_workflows.plugins import get_entry_point_name_from_class
from aiida_common_workflows.plugins import load_workflow_entry_point

PLUGIN_NAME = 'quantum_espresso'
CODE_LABEL = 'qe-6.8-pw@eiger-mc'
SET_NAME = 'set2'

STRUCTURES_GROUP_LABEL = f'commonwf-oxides/{SET_NAME}/structures/{PLUGIN_NAME}'
WORKFLOWS_GROUP_LABEL = f'commonwf-oxides/{SET_NAME}/workflows/{PLUGIN_NAME}'

Structure = DataFactory('structure')

query = orm.QueryBuilder()
query.append(Structure, tag='structure', project=['extras', '*'])
query.append(orm.Group, tag='group', filters={'label': STRUCTURES_GROUP_LABEL}, with_node='structure')
all_structures = {(res[0]['element'], res[0]['configuration']): res[1] for res in query.all()}

structure = all_structures[('Si', 'X2O')]
print(f'Structure PK: {structure.pk}')

sub_process_cls = load_workflow_entry_point('relax', 'quantum_espresso')
sub_process_cls_name = get_entry_point_name_from_class(sub_process_cls).name
generator = sub_process_cls.get_input_generator()

engine_types = generator.spec().inputs['engines']
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
            # 'queue_name': 'debug',
            'max_wallclock_seconds': 1700 # A bit less than 30 minutes (so we fit in the debug queue=partition)
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
            ## In order to make this work, you have perform the change discussed
            ## at the bottom of the file in the eos.py file.
            ## Otherwise the whole namespace is replaced.
            'pw': {
                'settings': orm.Dict(dict={
                    'cmdline': ['-nk', '32'],
                }), 
                'parameters': orm.Dict(dict={
                    'SYSTEM': {
                        'ecutwfc': 200,
                        'ecutrho': 2000
                    }
                })
            }
        }
    }
}

cls = WorkflowFactory('common_workflows.eos')
node = submit(cls, **inputs)
print(f"Submitted workflow with PK = {node.pk} for {structure.get_formula()}")


"""
NOTE: in order to make the recursive merge work (while we wait for the new implementation
of overrides), and not have the sub_process replace entire nodes or dicts,
you need to replace, in eos.py, in the method get_sub_workchain_builder, the call to

  builder._update(...)

with the following code:

        #### BEGIN REPLACEMENT ####
        # Instead of using _update (that completely replaces the topmost namespace)
        # or _merge (that merges recursively, but still replaces the whole Dict node)
        # we now use a strategy to recursively merge like _merge, but then for a Dict
        # we return a copy, with .update() called for its content dict.

        def _recursive_dict_update(dictionary, key, value):
            import collections

            if isinstance(value, collections.abc.Mapping):
                for inner_key, inner_value in value.items():
                    _recursive_dict_update(dictionary[key], inner_key, inner_value)
            else:
                dictionary[key] = value

        def _recursive_merge(dictionary, key, value):
            import collections            
            import copy

            if isinstance(value, collections.abc.Mapping):
                for inner_key, inner_value in value.items():
                    _recursive_merge(dictionary[key], inner_key, inner_value)
            else:
                if key in dictionary and isinstance(value, orm.Dict) and isinstance(dictionary[key], orm.Dict):
                    # If I want to replace a Dict node with a Dict node, instead
                    # of doing it, I recursively udpate the values *inside* the dict.
                    # I also check if key is in dictionary: if the dict was not existing before, we just go
                    # in the else statement, i.e. we add it.
                    new_dict = copy.deepcopy(dictionary[key].get_dict())
                    updating_dict = value.get_dict()
                    for inner_key, inner_value in updating_dict.items():
                        _recursive_dict_update(new_dict, inner_key, inner_value)
                    dictionary[key] = orm.Dict(dict=new_dict)
                else:
                    dictionary[key] = value

        for key, value in self.inputs.get('sub_process', {}).items():
            _recursive_merge(builder, key, value)
        ### END REPLACEMENT
"""
