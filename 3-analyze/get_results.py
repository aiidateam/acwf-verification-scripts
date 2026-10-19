#!/usr/bin/env runaiida

# The version of the script will be placed in the json file containing the results.
# We should change this number anytime this script or `eos_utils.eosfit_31_adapted` is modified.
import sys
import json
import os

import tqdm
import numpy as np

from collections import Counter
from eos_utils.eosfit_31_adapted import BM, echarge

from aiida import orm
from aiida.common import LinkType, NotExistentAttributeError
from aiida_common_workflows.workflows.relax.workchain import CommonRelaxWorkChain


__version__ = "0.0.4"

def extract_from_failed(node):
    """
    For EoS workchain with exit status different than zero, try to extract info on the completed volumes.

    This can be done only if the calculated volumes have a "relaxed_structure" output.
    In fact the EoS workchain calles <Code>CommonRelaxWorkChains and they have a common
    output interface, not an input one!
    """
    ens=[]
    vols=[]
    streses=[]
    num_atoms = None
    num_attempt_vols = 0
    for i in node.base.links.get_outgoing(link_type=LinkType.CALL_WORK).all():
        num_attempt_vols = num_attempt_vols + 1
        if i.node.is_finished_ok:
            if num_atoms is None:
                try:
                    num_atoms = len(i.node.outputs.relaxed_structure.sites)
                except NotExistentAttributeError:
                    # there is no output structure, can not retrieve the volume!
                    break
            ens.append(i.node.outputs.total_energy.value)
            vol = i.node.outputs.relaxed_structure.get_cell_volume()
            vols.append(vol)
            try:
                stress = i.node.outputs.stress.get_array('stress').tolist()
            except AttributeError:
                stress = None
            stresses.append(stress)

    return vols,ens,streses,num_atoms,num_attempt_vols


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

if __name__ == "__main__":
    try:
        SET_NAME = sys.argv[1]
    except IndexError:
        print("Pass as parameter the set name, e.g. oxides-verification-PBE-v1 or unaries-verification-PBE-v1")
        sys.exit(1)

    STRUCTURES_GROUP_LABEL = f'acwf-verification/{SET_NAME}/structures/{PLUGIN_NAME}'
    WORKFLOWS_GROUP_LABEL = f'acwf-verification/{SET_NAME}/workflows/{PLUGIN_NAME}'

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
    all_missing_outputs = {}
    completely_off = []
    failed_wfs = []
    all_eos_data = {}
    all_stress_data = {}
    all_BM_fit_data = {}
    num_atoms_in_sim_cell = {}

    # Initialize the progress bar as a variable so we can dynamically set its description
    progress_bar = tqdm.tqdm(wf_nodes)
    for node in progress_bar:
        structure = node.inputs.structure
        element = structure.base.extras.all['element']
        configuration = structure.base.extras.all['configuration']

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
        stress_data = None
        BM_fit_data = None
        num_atoms = None

        volumes = []
        energies = []
        stresses = []

        # For successfully finished workflows, collect the data from outputs
        if node.process_state.value == 'finished' and node.exit_status == 0:
            # Extract volumes and energies for this system
            outputs = node.base.links.get_outgoing(link_type=LinkType.RETURN).nested()
            for index, sub_structure in sorted(outputs['structures'].items()):
                if num_atoms is None:
                    num_atoms = len(sub_structure.sites)
                else:
                    assert num_atoms == len(sub_structure.sites), (
                        f"Number of atoms changes between structures for {element} {configuration}!"
                    )
                volumes.append(sub_structure.get_cell_volume())
                energy_node = outputs['total_energies'][index]
                energies.append(energy_node.value)
                parent_workflows_links = energy_node.base.links.get_incoming(link_type=LinkType.RETURN).all()
                parent_workflows = [
                    triple.node for triple in parent_workflows_links
                    if issubclass(triple.node.process_class, CommonRelaxWorkChain)]
                assert len(parent_workflows) == 1, "Error retrieving the parent Relax workflow!"
                parent_workflow = parent_workflows[0]
                try:
                    stress = parent_workflow.outputs.stress.get_array('stress').tolist()
                except AttributeError:
                    stress = None
                stresses.append(stress)
        # For failed workflows, check if some volumes concluded succesfully, if more than 80% of vol are ok, go on with fit
        elif (node.process_state.value == 'finished' and node.exit_status != 0) or (node.process_state.value == 'excepted'):
            volumes,energies,stresses,num_atoms,num_attempt_vols = extract_from_failed(node)
            if len(volumes)/float(num_attempt_vols) < 0.8:
                # Not enough volumes, list the material as failed
                failed_wfs.append({
                    'element': element,
                    'configuration': configuration,
                    'process_state': node.process_state.value,
                    'exit_status': node.exit_status,
                })
                # Return the info collected so far, eos_data, stress_data, BM_fit_data are still None
                all_eos_data[f'{element}-{configuration}'] = eos_data
                num_atoms_in_sim_cell[f'{element}-{configuration}'] = num_atoms
                all_stress_data[f'{element}-{configuration}'] = stress_data
                all_BM_fit_data[f'{element}-{configuration}'] = BM_fit_data
                # Exit loop = no fit attempted
                continue
            all_missing_outputs[f'{element}-{configuration}'] = num_attempt_vols-len(volumes)
            warning_lines.append(f"  WARNING! MISSING OUTPUTS: {num_attempt_vols-len(volumes)}")
        # We return all None also for the case not covered by the two others if. For instance, fall here the materials still running.
        else:
            #print(element,configuration,"probably still running?")
            all_eos_data[f'{element}-{configuration}'] = eos_data
            num_atoms_in_sim_cell[f'{element}-{configuration}'] = num_atoms
            all_stress_data[f'{element}-{configuration}'] = stress_data
            all_BM_fit_data[f'{element}-{configuration}'] = BM_fit_data
            # Exit loop = no fit attempted
            continue
        
        energies = [e for _, e in sorted(zip(volumes, energies))]
        volumes = sorted(volumes)
        # List as I need to JSON-serialize it
        eos_data = (np.array([volumes, energies]).T).tolist()
        stress_data = list(zip(volumes, stresses))

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
            data_to_print[
                (
                    structure.base.extras.all['element'], structure.base.extras.all['configuration']
                    )
                ] = (
                    min_volume, E0, bulk_modulus_GPa, bulk_deriv
                    )
            BM_fit_data = {
                'min_volume': min_volume,
                'E0': E0,
                'bulk_modulus_ev_ang3': bulk_modulus_ev_ang3,
                'bulk_deriv': bulk_deriv,
                'residuals': residuals[0]
            }
            if residuals[0] > 1.e-3:
                warning_lines.append(
                    f"WARNING! High fit residuals: {residuals[0]} for {structure.base.extras.all['element']} "
                    f"{structure.base.extras.all['configuration']}"
                    )
        except ValueError:
            # If we cannot find a minimum
            # Note that BM_fit_data was already set to None at the top
            warning_lines.append(
                f"WARNING! Unable to fit for {structure.base.extras.all['element']} {structure.base.extras.all['configuration']}"
                )

        all_eos_data[f'{element}-{configuration}'] = eos_data
        num_atoms_in_sim_cell[f'{element}-{configuration}'] = num_atoms
        all_stress_data[f'{element}-{configuration}'] = stress_data
        all_BM_fit_data[f'{element}-{configuration}'] = BM_fit_data

    data = {
        'script_version': __version__,
        'set_name': SET_NAME,
        # Mapping from strings like "He-X2O" to a dictionary with the UUIDs of the structure and the EOS workflow
        'uuid_mapping': uuid_mapping,
        # A list of dictionaries with information on the workchains that did not finish with a 0 exit code
        'failed_wfs': failed_wfs,
        # A dictionary that indicate for which elements and configurations there are missing outputs,
        # (only for the workchains that still had enough volumes to be considered for a fit)
        'missing_outputs': all_missing_outputs,
        # A list of dictionaries that indicate which elements and configurations have been computed completely
        # off-centre (meaning that the minimum of all computed energies is on either of the two edges, i.e. for
        # the smallest or largest volume)
        'completely_off': completely_off,
        # Dictionary with the EOS data (volumes and energies datapoints). The keys are the same as the `uuid_mapping`.
        # Values can be None.
        'eos_data': all_eos_data,
        'stress_data': all_stress_data,
        # Birch-Murnaghan fit data. See above for the keys. Can be None.
        'BM_fit_data': all_BM_fit_data,
        'num_atoms_in_sim_cell': num_atoms_in_sim_cell
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

    fname = f"outputs/warnings-{SET_NAME}-{PLUGIN_NAME}.txt"
    with open(fname, 'w') as fhandle:
        for line in warning_lines:
            fhandle.write(f"{line}\n")
            print(line)
    print(f"Warning log written to: '{fname}'.")

    # Output results to file
    os.makedirs('outputs', exist_ok=True)
    fname = f"outputs/results-{SET_NAME}-{PLUGIN_NAME}.json"
    with open(fname, 'w') as fhandle:
        json.dump(data, fhandle, indent=2, sort_keys=True)
    print(f"Output results written to: '{fname}'.")
