#!/usr/bin/env runaiida
import time

#from re import S
from aiida.plugins import DataFactory, WorkflowFactory
from aiida import orm
from aiida.engine import submit

import numpy as np

from aiida_common_workflows.common import ElectronicType, RelaxType, SpinType
from aiida_common_workflows.plugins import get_entry_point_name_from_class
from aiida_common_workflows.plugins import load_workflow_entry_point
from aiida_submission_controller import FromGroupSubmissionController

DRY_RUN = False
MAX_CONCURRENT = 200
PLUGIN_NAME = 'quantum_espresso'
CODE_LABEL = 'qe-7.0-pw@daint-mc'
# NPOOLS = 32 # for eiger (128 CPUS)
NPOOLS = 9 # for daint (36 CPUs)
SET_NAME = 'test-Er-Diamond0'

class EosSubmissionController(FromGroupSubmissionController):
    """A SubmissionController for submitting EOS with Quantum ESPRESSO common workflows."""
    def __init__(self, code_label, smearing_type, degauss_meV, *args, **kwargs):
        """Pass also a code label, that should be a code associated to an `quantumespresso.pw` plugin."""
        super().__init__(*args, **kwargs)
        self._code = orm.load_code(code_label)
        self._process_class = WorkflowFactory('common_workflows.eos')
        self._smearing_type = smearing_type
        self._degauss_meV = degauss_meV

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
                    'max_wallclock_seconds': 3600
                }
            }

        dim_data = structure.get_dimensionality()
        assert dim_data['dim'] == 3
        reference_volume = dim_data['value']

        actual_volumes = np.linspace(45, 200, 40)
        scale_factors = actual_volumes / reference_volume

        inputs = {
            'structure': structure,
            'scale_factors': orm.List(list=scale_factors.tolist()),
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
                            'cmdline': ['-nk', str(NPOOLS)],
                        }),
                        # This works because we patched the eos.py code
                        'parameters': orm.Dict(dict = {
                            'SYSTEM': {
                                'smearing': self._smearing_type,
                                'degauss': self._degauss_meV / 13605.662285137 # convert to Ry
                            }
                        })
                    }
                }
            }
        }

        return inputs, self._process_class

if __name__ == "__main__":
    for smearing_type, degauss_meV_values in [
            ['fermi-dirac', [30., 50., 70., 90., 110.]],
            ['gaussian', [50., 70., 90., 110.]],
            ['cold', [50., 70., 90., 110.]],
        ]:
        for degauss_meV in degauss_meV_values:

            print("#"*72)
            print(f"#### SMEARING = {smearing_type} -> DEGAUSS = {degauss_meV}")
            print("#"*72)

            STRUCTURES_GROUP_LABEL = f'commonwf-oxides/{SET_NAME}/structures/{PLUGIN_NAME}'
            WORKFLOWS_GROUP_LABEL = f'commonwf-oxides/{SET_NAME}/workflows/{PLUGIN_NAME}/smearing-{smearing_type}/{degauss_meV}'

            controller = EosSubmissionController(
                parent_group_label=STRUCTURES_GROUP_LABEL,
                code_label=CODE_LABEL,
                group_label=WORKFLOWS_GROUP_LABEL,
                max_concurrent=MAX_CONCURRENT,
                smearing_type=smearing_type,
                degauss_meV=degauss_meV
            )
            
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
