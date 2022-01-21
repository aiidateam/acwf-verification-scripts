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
MAX_CONCURRENT = 200
PLUGIN_NAME = 'quantum_espresso'
CODE_LABEL = 'qe-6.8-pw@eiger-mc'
SET_NAME = 'monoelemental-structures-test'
SUFFIX = ''
ACCOUNT = 'mr0'
                    
#MODIFICATIONS
#CODE_LABEL = 'qe-6.8-pw@imxge'
#SUFFIX = '/imxge'
#ACCOUNT = None

#CODE_LABEL = 'qe-6.7-pw@thanos'
#SUFFIX = '/thanos'
#ACCOUNT = None
#END OF MODIFICATIONS


STRUCTURES_GROUP_LABEL = f'commonwf-oxides/{SET_NAME}' # I could add here /{PLUGIN_NAME} if I also filter by those working for the current code; for now I do them all.
WORKFLOWS_GROUP_LABEL = f'commonwf-oxides/{SET_NAME}/workflows/{PLUGIN_NAME}{SUFFIX}'

# This table is obtained with the script utils/get_ecuts_table_qe.py (for set 1)
# THIS IS ONLY VALID FOR THE USED PSEUDO LIBRARY AND PROTOCOL!
# TODO: replace this with a call to the protocol functions to get the cutoff instead,
# in the function below get_ecuts.
# utils/get_ecuts_table_qe.py can still be used as a further check, to assert that
# all calculatons of a given oxide family with the same element used the same cutoffs.
ECUTS_TABLE = {
    'Ag': (75.0, 600.0),
    'Al': (75.0, 600.0),
    'Ar': (120.0, 600.0),
    'As': (75.0, 600.0),
    'Au': (75.0, 600.0),
    'B': (75.0, 600.0),
    'Ba': (75.0, 600.0),
    'Be': (75.0, 600.0),
    'Bi': (75.0, 600.0),
    'Br': (90.0, 720.0),
    'C': (75.0, 600.0),
    'Ca': (75.0, 600.0),
    'Cd': (90.0, 720.0),
    'Ce': (75.0, 600.0),
    'Cl': (100.0, 800.0),
    'Co': (90.0, 1080.0),
    'Cr': (75.0, 600.0),
    'Cs': (75.0, 600.0),
    'Cu': (90.0, 600.0),
    'Dy': (75.0, 600.0),
    'Er': (75.0, 600.0),
    'Eu': (75.0, 600.0),
    'F': (90.0, 600.0),
    'Fe': (90.0, 1080.0),
    'Ga': (90.0, 720.0),
    'Gd': (75.0, 600.0),
    'Ge': (75.0, 600.0),
    'H': (80.0, 600.0),
    'He': (75.0, 600.0),
    'Hf': (75.0, 600.0),
    'Hg': (75.0, 600.0),
    'Ho': (75.0, 600.0),
    'I': (75.0, 600.0),
    'In': (75.0, 600.0),
    'Ir': (75.0, 600.0),
    'K': (75.0, 600.0),
    'Kr': (75.0, 600.0),
    'La': (75.0, 600.0),
    'Li': (75.0, 600.0),
    'Lu': (75.0, 600.0),
    'Mg': (75.0, 600.0),
    'Mn': (90.0, 1080.0),
    'Mo': (75.0, 600.0),
    'N': (80.0, 600.0),
    'Na': (100.0, 600.0),
    'Nb': (75.0, 600.0),
    'Nd': (75.0, 600.0),
    'Ne': (75.0, 600.0),
    'Ni': (75.0, 600.0),
    'Os': (75.0, 600.0),
    'P': (75.0, 600.0),
    'Pb': (75.0, 600.0),
    'Pd': (75.0, 600.0),
    'Pm': (75.0, 600.0),
    'Po': (80.0, 640.0),
    'Pr': (75.0, 600.0),
    'Pt': (100.0, 800.0),
    'Rb': (75.0, 600.0),
    'Re': (75.0, 600.0),
    'Rh': (75.0, 600.0),
    'Rn': (200.0, 1600.0),
    'Ru': (75.0, 600.0),
    'S': (75.0, 600.0),
    'Sb': (75.0, 600.0),
    'Sc': (90.0, 720.0),
    'Se': (75.0, 600.0),
    'Si': (75.0, 600.0),
    'Sm': (75.0, 600.0),
    'Sn': (75.0, 600.0),
    'Sr': (75.0, 600.0),
    'Ta': (75.0, 600.0),
    'Tb': (75.0, 600.0),
    'Tc': (75.0, 600.0),
    'Te': (75.0, 600.0),
    'Ti': (75.0, 600.0),
    'Tl': (75.0, 600.0),
    'Tm': (75.0, 600.0),
    'V': (75.0, 600.0),
    'W': (75.0, 600.0),
    'Xe': (80.0, 600.0),
    'Y': (75.0, 600.0),
    'Yb': (75.0, 600.0),
    'Zn': (90.0, 720.0),
    'Zr': (75.0, 600.0)
 }

def get_ecuts(element, configuration):
    """Return the two cutoffs for monoelemental solids.
    
    :return: a tuple of (ecutwfc, ecutrho).

    :note: for simplicity, I create one oxygen per other element, and for each of them use the cutoffs.
        In reality, already from the data above, we see that we have only few possible combinations of cutoffs, e.g.
        only 10 in the example above: 
           (75,600), (80,600), (80,640), (90,600), (90,720), (90,1080), (100,600), (100,800), (120,600), (200,1600)
        Therefore one can think to a different approach, e.g. using the `configuration` to have instead the two cutoffs
        separated by a dash, and parsing those to return the cutoffs, like:

          if element == "O":
              ecutwfc_str, ecutrho_str = configuration.split('-')
              ecutwfc = float(ecutwfc_str)
              ecutrho = float(ecutrho_str)
        
        Then also the plotting routines need to be updated.

    For elements other than oxygen, you should return the cutoffs used for element + oxygen.
    For oxygen, you should check the configuration to know what the other element is, and you should return
    the cutoffs used for other_element + oxygen.

    You can check the cutoffs looking at the ones already used, or you call through to the correct
    methods to get the cutoff via the protocol.
    """
    if element == "O":
        # The configuration, in the case of monoelemental oxygen, represents the other element, e.g. 'Na'
        # indicates oxygen with sodium.
        # Therefore, I use that as the key.
        ecutwfc, ecutrho = ECUTS_TABLE[configuration]
    else:
        assert configuration == "X", f"The configuration for monoelemental '{element}' is not 'X'"
        ecutwfc, ecutrho = ECUTS_TABLE[element]

    return ecutwfc, ecutrho

class EosSubmissionController(FromGroupSubmissionController):
    """A SubmissionController for submitting EOS with Quantum ESPRESSO common workflows."""
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

        sub_process_cls = load_workflow_entry_point('relax', 'quantum_espresso')
        sub_process_cls_name = get_entry_point_name_from_class(sub_process_cls).name
        generator = sub_process_cls.get_input_generator()

        ecutwfc, ecutrho = get_ecuts(*extras_values)

        # Here I use the old syntax since I didn't update the code of the common workflows, yet
        engine_types = generator.spec().inputs['engines']
        engines = {}
        # There should be only one
        for engine in engine_types:
            engine_dict = {
                'code': CODE_LABEL,
                'options': {
                    'resources': {
                        'num_machines': 1
                    },
                    'max_wallclock_seconds': 3600
                }
            }
            if ACCOUNT is not None:
                engine_dict['options']['account'] = ACCOUNT
            engines[engine] = engine_dict

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
                ### IMPORTANT! TO MAKE THIS WORK, YOU HAVE TO MAKE THE CHANGE TO eos.py
                ### AS DESCRIBED AT THE BOTTOM OF THE FILE first-tests/launch_example_one_calc_only_template.py
                'base': {
                    'pw': {
                        'settings' : orm.Dict(dict={
                            'cmdline': ['-nk', '32'],
                        }), 
                        'parameters': orm.Dict(dict={
                            'SYSTEM': {
                                'ecutwfc': ecutwfc,
                                'ecutrho': ecutrho
                            }
                        })
                    }
                }
            }
        }

        return inputs, self._process_class

if __name__ == "__main__":
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
