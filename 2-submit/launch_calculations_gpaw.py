import sys
import time

#from re import S
from aiida.plugins import DataFactory, WorkflowFactory
from aiida import orm
from aiida.engine import submit

from aiida_common_workflows.common import ElectronicType, RelaxType, SpinType
from aiida_common_workflows.plugins import get_entry_point_name_from_class
from aiida_common_workflows.plugins import load_workflow_entry_point
from aiida_submission_controller import FromGroupSubmissionController

DRY_RUN = False
MAX_CONCURRENT = 30
PLUGIN_NAME = 'gpaw'
CODE_LABEL = 'gpaw.21.6.0@dtu_xeon40'

class EosSubmissionController(FromGroupSubmissionController):
    """A SubmissionController for submitting EOS with GPAW common workflows."""
    def __init__(self, code_label, *args, **kwargs):
        """Pass also a code label, that should be a code associated to an `quantumespresso.pw` plugin."""
        super().__init__(*args, **kwargs)
        self._code = orm.load_code(code_label)
        self._process_class = WorkflowFactory('common_workflows.eos')

    def get_extra_unique_keys(self):
        """Return a tuple of the keys of the unique extras that will be used to uniquely identify your workchains.

        Here: the chemical symbol of the element, and the configuration (XO, XO2, X2O3, ...).
        """
        return ['element', 'configuration']

    def get_inputs_and_processclass_from_extras(self, extras_values):
        """Return inputs and process class for the submission of this specific process.

        I just submit an ArithmeticAdd calculation summing the two values stored in the extras:
        ``left_operand + right_operand``.
        """
        structure = self.get_parent_node_from_extras(extras_values)

        sub_process_cls = load_workflow_entry_point('relax', 'gpaw')
        sub_process_cls_name = get_entry_point_name_from_class(sub_process_cls).name
        generator = sub_process_cls.get_input_generator()

        engine_types = generator.get_engine_types()
        engines = {}
        # There should be only one
        for engine in engine_types:
            engines[engine] = {
                'code': self._code,
                'options': {
                    'resources': {
                        'num_machines': 1
                    },
                    'withmpi':True,
                    'max_wallclock_seconds': 1 * 60 * 60, 
                }
            }

        generator.get_builder(structure, engines)

        inputs = {
            'structure': structure,
            'generator_inputs': {  # code-agnostic inputs for the relaxation
                'engines': engines,
                'protocol': 'oxides_validation',
                'relax_type': RelaxType.NONE,
                'electronic_type': ElectronicType.METAL,
                'spin_type': SpinType.NONE,
            },
            'sub_process_class': sub_process_cls_name,
        }

        return inputs, self._process_class

if __name__ == "__main__":
    try:
        SET_NAME = sys.argv[1]
    except IndexError:
        print("Pass as parameter the set name, e.g. oxides-verification-PBE-v1 or unaries-verification-PBE-v1")
        sys.exit(1)

    STRUCTURES_GROUP_LABEL = f'acwf-verification/{SET_NAME}/structures/{PLUGIN_NAME}'
    WORKFLOWS_GROUP_LABEL = f'acwf-verification/{SET_NAME}/workflows/{PLUGIN_NAME}'

    controller = EosSubmissionController(
        parent_group_label=STRUCTURES_GROUP_LABEL,
        code_label=CODE_LABEL,
        group_label=WORKFLOWS_GROUP_LABEL,
        max_concurrent=MAX_CONCURRENT)
    
    print('Already run    :', controller.num_already_run)
    print('Max concurrent :', controller.max_concurrent)
    print('Available slots:', controller.num_available_slots)
    print('Active slots   :', controller.num_active_slots)
    print('Still to run   :', controller.num_to_run)
    print()

    run_processes = controller.submit_new_batch(dry_run=DRY_RUN)
    for run_process_extras, run_process in run_processes.items():
        if run_process is None:
            print(f'{run_process_extras} --> To be run')    
        else:
            print(f'{run_process_extras} --> PK = {run_process.pk}')

    print()
