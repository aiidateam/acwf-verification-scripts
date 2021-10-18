#!/usr/bin/env runaiida

# The version of the script will be placed in the json file containing the results.
# We should change this number anytime this script or `eos_utils.eosfit_31_adapted` is modified.
__version__ = "0.0.2"

SET_NAME = 'set2'

from collections import Counter
import json
import os

import numpy as np
import tqdm
from eos_utils.eosfit_31_adapted import BM, echarge

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

from aiida import orm
from aiida.common import LinkType


if __name__ == "__main__":
    # Get all nodes in the output group (EOS workflows)
    group_node_query = orm.QueryBuilder().append(
        orm.Group, filters={'label': WORKFLOWS_GROUP_LABEL}, tag='groups',
    ).append(orm.Node, project='*', with_group='groups')
    group_node_query.distinct()
    wf_nodes = group_node_query.all(flat=True)

    states = []
    data_to_print = {}
    warning_lines = []

    uuid_mapping = {}
    all_missing_outputs = []
    completely_off = []
    failed_wfs = []
    all_eos_data = {}
    all_BM_fit_data = {}

    # Initialize the progress bar as a variable so we can dynamically set its description
    progress_bar = tqdm.tqdm(wf_nodes)
    for node in progress_bar:
        structure = node.inputs.structure
        element = structure.extras['element']
        configuration = structure.extras['configuration']

        uuid_mapping[f'{element}-{configuration}'] = {
            'structure': structure.uuid,
            'eos_workflow': node.uuid
        }

        # Set the progress bar description; using :16s to minimize length change of the description
        description = f"{element} {configuration} ({structure.pk})"
        progress_bar.set_description(f"{description:16s}")
        progress_bar.refresh()

        # Get the state (possibly adding the exit status if it's finished) and add to a list
        state = node.process_state.value
        if state == 'finished':
            state = f'{state}.{node.exit_status}'
        states.append(state)

        # Initialize to None if the outputs are not there
        eos_data = None
        BM_fit_data = None

        # For successfully finished workflows, fit the EOS
        if node.process_state.value == 'finished' and node.exit_status == 0:
            # Check if the required outputs are there
            outputs = node.get_outgoing(link_type=LinkType.RETURN).nested()
            missing_outputs = tuple(output for output in ('structures', 'total_energies') if output not in outputs)
            if missing_outputs:
                all_missing_outputs.append({'element': element, 'configuration': configuration, 'missing': missing_outputs})
                warning_lines.append(f"  WARNING! MISSING OUTPUTS: {missing_outputs}")
            else:
                # Extract volumes and energies for this system
                volumes = []
                energies = []
                for index, sub_structure in sorted(outputs['structures'].items()):
                    volumes.append(sub_structure.get_cell_volume())
                    energies.append(outputs['total_energies'][index].value)
                energies = [e for _, e in sorted(zip(volumes, energies))]
                volumes = sorted(volumes)
                # List as I need to JSON-serialize it
                eos_data = (np.array([volumes, energies]).T).tolist()

                # Check if the central point was completely off (i.e. the minimum of the energies is
                # on the very left or very right of the volume range)
                min_loc = np.array(energies).argmin()
                if min_loc == 0:
                    # Side is whether the minimum occurs on the left side (small volumes) or right side (large volumes)
                    completely_off.append({'element': element, 'configuration': configuration, 'side': 'left'})
                elif min_loc == len(energies) - 1:
                    completely_off.append({'element': element, 'configuration': configuration, 'side': 'right'})   

                try:
                    # I need to pass a numpy array
                    min_volume, E0, bulk_modulus_internal, bulk_deriv, residuals = BM(np.array(eos_data))
                    bulk_modulus_GPa = bulk_modulus_internal * echarge * 1.0e21
                    #1 eV/Angstrom3 = 160.21766208 GPa
                    bulk_modulus_ev_ang3 = bulk_modulus_GPa / 160.21766208
                    data_to_print[(structure.extras['element'], structure.extras['configuration'])] = (
                        min_volume, E0, bulk_modulus_GPa, bulk_deriv)
                    BM_fit_data = {
                        'min_volume': min_volume,
                        'E0': E0,
                        'bulk_modulus_ev_ang3': bulk_modulus_ev_ang3,
                        'bulk_deriv': bulk_deriv,
                        'residuals': residuals[0]
                    }
                    if residuals[0] > 1.e-6:
                        warning_lines.append(f"WARNING! High fit residuals: {residuals[0]} for {structure.extras['element']} {structure.extras['configuration']}")
                except ValueError:
                    # If we cannot find a minimum
                    # Note that BM_fit_data was already set to None at the top
                    warning_lines.append(f"WARNING! Unable to fit for {structure.extras['element']} {structure.extras['configuration']}")

        elif (node.process_state.value == 'finished' and node.exit_status != 0) or (node.process_state.value == 'excepted'):
            failed_wfs.append({
                'element': element,
                'configuration': configuration,
                'process_state': node.process_state.value,
                'exit_status': node.exit_status,
            })
        all_eos_data[f'{element}-{configuration}'] = eos_data
        all_BM_fit_data[f'{element}-{configuration}'] = BM_fit_data

    data = {
        'script_version': __version__,
        'set_name': SET_NAME,
        # Mapping from strings like "He-X2O" to a dictionary with the UUIDs of the structure and the EOS workflow
        'uuid_mapping': uuid_mapping,
        # A list of dictionaries with information on the workchains that did not finish with a 0 exit code
        'failed_wfs': failed_wfs,
        # A list of dictionaries that indicate for which elements and configurations there are missing outputs,
        # if any (for the workchains that finished with 0 exit code)
        'missing_outputs': all_missing_outputs,
        # A list of dictionaries that indicate which elements and configurations have been computed completely
        # off-centre (meaning that the minimum of all computed energies is on either of the two edges, i.e. for
        # the smallest or largest volume)
        'completely_off': completely_off,
        # Dictionary with the EOS data (volumes and energies datapoints). The keys are the same as the `uuid_mapping`.
        # Values can be None.
        'eos_data': all_eos_data,
        # Birch-Murnaghan fit data. See above for the keys. Can be None.
        'BM_fit_data': all_BM_fit_data
    }

    # Print some statistics on the results
    warning_lines.append("")
    warning_lines.append("Counter of states: " + str(Counter(states)))
    good_cnt = len([eos_data for eos_data in data['eos_data'] if eos_data is not None])
    warning_lines.append("")
    warning_lines.append(f"Minimum completely off for {len(completely_off)}/{good_cnt}")
    warning_lines.append("Completely off systems (symbol indicates if the minimum is on the very left or right):")
    for system in data['completely_off']:
        warning_lines.append(
            f"- {system['element']} {system['configuration']} "
            f"({'<' if system['side'] == 'left' else '>'})"
        )

    fname = f"outputs/results-warnings-{PLUGIN_NAME}.txt"
    with open(fname, 'w') as fhandle:
        for line in warning_lines:
            fhandle.write(f"{line}\n")
            print(line)
    print(f"Warning log written to: '{fname}'.")

    # Output results to file
    os.makedirs('outputs', exist_ok=True)
    fname = f"outputs/results-{PLUGIN_NAME}.json"
    with open(fname, 'w') as fhandle:
        json.dump(data, fhandle, indent=2, sort_keys=True)
    print(f"Output results written to: '{fname}'.")
