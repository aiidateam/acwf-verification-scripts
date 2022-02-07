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
SET_NAME = 'monoelemental-structures-test2'
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

STRUCTURES_GROUP_LABEL = f'commonwf-oxides/{SET_NAME}/structures/{PLUGIN_NAME}'
WORKFLOWS_GROUP_LABEL = f'commonwf-oxides/{SET_NAME}/workflows/{PLUGIN_NAME}{SUFFIX}'


# THIS MUST BE THE SAME AS IN THE 1-preliminary/create_starting_subgroup_monoelemental.py file!
QE_ECUTS_TABLE = {
    'Ag': (75, 600),
    'Al': (75, 600),
    'Ar': (120, 600),
    'As': (75, 600),
    'Au': (75, 600),
    'B': (75, 600),
    'Ba': (75, 600),
    'Be': (75, 600),
    'Bi': (75, 600),
    'Br': (90, 720),
    'C': (75, 600),
    'Ca': (75, 600),
    'Cd': (90, 720),
    'Ce': (75, 600),
    'Cl': (100, 800),
    'Co': (90, 1080),
    'Cr': (75, 600),
    'Cs': (75, 600),
    'Cu': (90, 600),
    'Dy': (75, 600),
    'Er': (75, 600),
    'Eu': (75, 600),
    'F': (90, 600),
    'Fe': (90, 1080),
    'Ga': (90, 720),
    'Gd': (75, 600),
    'Ge': (75, 600),
    'H': (80, 600),
    'He': (75, 600),
    'Hf': (75, 600),
    'Hg': (75, 600),
    'Ho': (75, 600),
    'I': (75, 600),
    'In': (75, 600),
    'Ir': (75, 600),
    'K': (75, 600),
    'Kr': (75, 600),
    'La': (75, 600),
    'Li': (75, 600),
    'Lu': (75, 600),
    'Mg': (75, 600),
    'Mn': (90, 1080),
    'Mo': (75, 600),
    'N': (80, 600),
    'Na': (100, 600),
    'Nb': (75, 600),
    'Nd': (75, 600),
    'Ne': (75, 600),
    'Ni': (75, 600),
    'Os': (75, 600),
    'P': (75, 600),
    'Pb': (75, 600),
    'Pd': (75, 600),
    'Pm': (75, 600),
    'Po': (80, 640),
    'Pr': (75, 600),
    'Pt': (100, 800),
    'Rb': (75, 600),
    'Re': (75, 600),
    'Rh': (75, 600),
    'Rn': (200, 1600),
    'Ru': (75, 600),
    'S': (75, 600),
    'Sb': (75, 600),
    'Sc': (90, 720),
    'Se': (75, 600),
    'Si': (75, 600),
    'Sm': (75, 600),
    'Sn': (75, 600),
    'Sr': (75, 600),
    'Ta': (75, 600),
    'Tb': (75, 600),
    'Tc': (75, 600),
    'Te': (75, 600),
    'Ti': (75, 600),
    'Tl': (75, 600),
    'Tm': (75, 600),
    'V': (75, 600),
    'W': (75, 600),
    'Xe': (80, 600),
    'Y': (75, 600),
    'Yb': (75, 600),
    'Zn': (90, 720),
    'Zr': (75, 600)
}

def get_ecuts(element, configuration):
    """Return the two cutoffs for monoelemental solids.
    
    :return: a tuple of (ecutwfc, ecutrho).

    :note: for QE, for oxygen, the configuration is of the form XXX-YYY where XXX and YYY are
        ecutwfc and ecutrho, respectively.
        There are converted to floats. For the other elements, the configuration is always
        the string 'X'.

    For elements other than oxygen, you should return the cutoffs used for element + oxygen.
    For oxygen, you should check the configuration to know what the other element is, and you should return
    the cutoffs used for other_element + oxygen.

    You can check the cutoffs looking at the ones already used, or you call through to the correct
    methods to get the cutoff via the protocol.
    """
    if element == "O":
        # The configuration, in the case of monoelemental oxygen, 
        # represents the two cutoffs
        ecutwfc_str, ecutrho_str = configuration.split('-')
        ecutwfc = float(ecutwfc_str)
        ecutrho = float(ecutrho_str)
    else:
        assert configuration == "X", f"The configuration for monoelemental '{element}' is not 'X'"
        ecutwfc, ecutrho = QE_ECUTS_TABLE[element]

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
