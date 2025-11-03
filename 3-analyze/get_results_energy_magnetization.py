#!/usr/bin/env runaiida

import sys
import json
import os

import tqdm

from aiida import orm
from aiida.common import LinkType


__version__ = "0.0.4"

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
        WORKFLOWS_GROUP_LABEL = f"acwf-verification/{SET_NAME}/workflows"
    except IndexError:
        print("Pass as parameter the set name and the workflows group label.")
        sys.exit(1)


    # Get all nodes in the output group (EOS workflows)
    group_node_query = orm.QueryBuilder().append(
        orm.Group, filters={'label': WORKFLOWS_GROUP_LABEL}, tag='groups',
    ).append(orm.Node, project='*', with_group='groups')
    group_node_query.distinct()
    wf_nodes = group_node_query.all(flat=True)

    uuid_mapping = {}
    failed_wfs = []
    all_total_magnetizations_energies = {}
    all_fermi_energies = {}

    progress_bar = tqdm.tqdm(wf_nodes)
    for node in progress_bar:
        structure = node.inputs.structure
        element = structure.extras['element']
        configuration = structure.extras['configuration']

        uuid_mapping[f'{element}-{configuration}'] = {
            'structure': structure.uuid,
            'eos_workflow': node.uuid
        }

        description = f"{element} {configuration} ({structure.pk})"
        progress_bar.set_description(f"{description:16s}")
        progress_bar.refresh()

        total_magnetizations_energies = []
        fermi_energies = []
        num_atoms = len(node.inputs.structure.sites)

        # For successfully finished workflows, collect the data from outputs
        if node.is_finished_ok:
            outputs = node.get_outgoing(link_type=LinkType.RETURN).nested()
            for index, tot_magnetization in sorted(outputs['total_magnetizations'].items()):
                energy_node = outputs['total_energies'][index]
                total_magnetizations_energies.append((tot_magnetization.value, energy_node.value))
                fermi_energy_up = outputs['fermi_energies_up'][index]
                fermi_energy_down = outputs['fermi_energies_down'][index]
                fermi_energies.append((fermi_energy_down.value, fermi_energy_up.value))

            all_total_magnetizations_energies[f'{element}-{configuration}'] = total_magnetizations_energies
            all_fermi_energies[f'{element}-{configuration}'] = fermi_energies
        else:
            print(f"Workflow {node.pk} (configuration: {element}-{configuration}) did not finish successfully (exit status: {node.exit_status}).")
            failed_wfs.append({
                    'element': element,
                    'configuration': configuration,
                    'process_state': node.process_state.value,
                    'exit_status': node.exit_status,
                })

    data = {
        'total_magnetization_energy_data': all_total_magnetizations_energies,
        'fermi_energies': all_fermi_energies,
        'script_version': __version__,
        'uuid_mapping': uuid_mapping,
        'failed_wfs': failed_wfs,
        'set_name': SET_NAME,
    }

    os.makedirs('outputs', exist_ok=True)
    fname = f"outputs/results-{SET_NAME}-{PLUGIN_NAME}.json"
    with open(fname, 'w') as fhandle:
        json.dump(data, fhandle, indent=2, sort_keys=True)
    print(f"Output results written to: '{fname}'.")
