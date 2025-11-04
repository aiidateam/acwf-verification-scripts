#!/usr/bin/env runaiida
import json

# from re import S
from aiida.plugins import WorkflowFactory

#!/usr/bin/env runaiida
from aiida import orm

from aiida_common_workflows.common import ElectronicType, RelaxType, SpinType
from aiida_common_workflows.plugins import get_entry_point_name_from_class
from aiida_common_workflows.plugins import load_workflow_entry_point
from aiida_submission_controller import FromGroupSubmissionController

DRY_RUN = False
MAX_CONCURRENT = 1
PLUGIN_NAME = "quantum_espresso"
CODE_LABEL = "qe-7.3-gf-pw@thor"


class EnergyMagnetizationSubmissionController(FromGroupSubmissionController):
    """A SubmissionController for submitting EOS with Quantum ESPRESSO common workflows."""

    def __init__(
        self, code_label, configuration_total_magnetizations_mapping, *args, **kwargs
    ):
        """Pass also a code label, that should be a code associated to an `quantumespresso.pw` plugin.

        Also pass a mapping from configuration, e.g., `Ni-X/FCC` to a list of total magnetizations to be used
        in the energy vs magnetization calculations.
        """
        super().__init__(*args, **kwargs)
        self._code = orm.load_code(code_label)
        self._process_class = WorkflowFactory("common_workflows.em")
        self._configuration_total_magnetizations_mapping = (
            configuration_total_magnetizations_mapping
        )

    def get_extra_unique_keys(self):
        """Return a tuple of the keys of the unique extras that will be used to uniquely identify your workchains.

        Here: the chemical symbol of the element, and the configuration (XO, XO2, X2O3, ...).
        """
        return ["element", "configuration"]

    def get_inputs_and_processclass_from_extras(self, extras_values):
        """Return inputs and process class for the submission of this specific process.

        I just submit an ArithmeticAdd calculation summing the two values stored in the extras:
        ``left_operand + right_operand``.
        """
        structure = self.get_parent_node_from_extras(extras_values)

        sub_process_cls = load_workflow_entry_point("relax", "quantum_espresso")
        sub_process_cls_name = get_entry_point_name_from_class(sub_process_cls).name
        generator = sub_process_cls.get_input_generator()

        key = f"{extras_values[0]}-{extras_values[1]}"
        total_magnetizations = self._configuration_total_magnetizations_mapping[key][
            "fixed_total_magnetizations_mu_B/cell"
        ]

        engine_types = generator.spec().inputs["engines"]
        engines = {}
        # There should be only one
        for engine in engine_types:
            engines[engine] = {
                "code": CODE_LABEL,
                "options": {
                    "resources": {
                        "num_machines": 1,
                        "num_mpiprocs_per_machine": 48,
                    },
                    "queue_name": "short",
                    "max_wallclock_seconds": 1200,
                },
            }

        inputs = {
            "structure": structure,
            "fixed_total_magnetizations": total_magnetizations,
            "generator_inputs": {  # code-agnostic inputs for the relaxation
                "engines": engines,
                "protocol": "verification-PBE-v1",
                "relax_type": RelaxType.NONE,
                "electronic_type": ElectronicType.METAL,
                "spin_type": SpinType.COLLINEAR,
            },
            "sub_process_class": sub_process_cls_name,
            "sub_process": {  # optional code-dependent overrides
                "base": {
                    "pw": {
                        "settings": orm.Dict(
                            {
                                "cmdline": ["-nk", "8"],
                            }
                        )
                    }
                }
            },
        }

        return inputs, self._process_class


if __name__ == "__main__":

    SET_NAME = "unaries-verification-PBE-magnetic-d-block-OFClBr-v1"
    STRUCTURES_GROUP_LABEL = f"acwf-verification/{SET_NAME}/structures"
    WORKFLOWS_GROUP_LABEL = f"acwf-verification/{SET_NAME}/workflows"

    with open(
        "../1-preliminary/magnetic-verification/fixed_total_magnetiations-unaries-PBE-magnetic-d-block-OFClBr-v1.json",
        "r",
    ) as f:
        configuration_total_magnetizations_mapping = json.load(f)

    controller = EnergyMagnetizationSubmissionController(
        parent_group_label=STRUCTURES_GROUP_LABEL,
        code_label=CODE_LABEL,
        configuration_total_magnetizations_mapping=configuration_total_magnetizations_mapping,
        group_label=WORKFLOWS_GROUP_LABEL,
        max_concurrent=MAX_CONCURRENT,
    )

    print("Already run    :", controller.num_already_run)
    print("Max concurrent :", controller.max_concurrent)
    print("Available slots:", controller.num_available_slots)
    print("Active slots   :", controller.num_active_slots)
    print("Still to run   :", controller.num_to_run)
    print()

    run_processes = controller.submit_new_batch(dry_run=DRY_RUN)
    for run_process_extras, run_process in run_processes.items():
        if run_process is None:
            print(f"{run_process_extras} --> To be run")
        else:
            print(f"{run_process_extras} --> PK = {run_process.pk}")

    print()
