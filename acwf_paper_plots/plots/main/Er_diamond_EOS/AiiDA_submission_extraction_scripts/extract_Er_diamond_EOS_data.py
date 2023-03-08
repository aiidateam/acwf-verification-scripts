#!/usr/bin/env runaiida

# The version of the script will be placed in the json file containing the results.
# We should change this number anytime this script or `eos_utils.eosfit_31_adapted` is modified.
import numpy as np
import json

from aiida import orm
from aiida.common import LinkType
from aiida_common_workflows.workflows.relax.workchain import CommonRelaxWorkChain
from scipy.interpolate import interp1d

SET_NAME = "test-Er-Diamond0"
PLUGIN_NAME = "quantum_espresso"
ONLY_FULLY_COMPLETE = False


smearing_name_mapping = {
    'fermi-dirac': 'FD',
    'cold': 'Cold',
    'gaussian': 'Gaussian'
}


if __name__ == "__main__":
    STRUCTURES_GROUP_LABEL = f'commonwf-oxides/{SET_NAME}/structures/{PLUGIN_NAME}'
    WORKFLOWS_GROUP_LABEL_PREFIX = f'commonwf-oxides/{SET_NAME}/workflows/{PLUGIN_NAME}/smearing-'

    group_labels = orm.QueryBuilder().append(
        orm.Group, filters={'label': {'like': f"{WORKFLOWS_GROUP_LABEL_PREFIX}%"}}, project='label').all(flat=True)

    def sort_labels(group_name):
        """Return A tuple with the degauss as a float, to have a proper sorting."""
        smearing_string = group_name[len(WORKFLOWS_GROUP_LABEL_PREFIX):]
        try:
            smearing_type, degauss_meV = smearing_string.split('/')
        except ValueError:
            return ('zz', 1e99)
        degauss_meV = float(degauss_meV)
        return(smearing_type, degauss_meV)

    data = []
    warning_lines = []
    for label in sorted(group_labels, key=sort_labels):
        smearing_string = label[len(WORKFLOWS_GROUP_LABEL_PREFIX):]
        try:
            smearing_type, degauss_meV = smearing_string.split('/')
        except ValueError:
            # E.g. if there are too many slashes, e.g. for /OLD systems
            continue
        degauss_meV = float(degauss_meV)

        if smearing_type == 'fermi-dirac' and abs(degauss_meV - 110.) < 1.e-6:
            print(f"SKIPPING MANUALLY {smearing_type} @ {degauss_meV} meV")
            continue
        print(smearing_type, degauss_meV)

        # Get all nodes in the output group (EOS workflows)
        group_node_query = orm.QueryBuilder().append(
            orm.Group, filters={'label': label}, tag='groups',
        ).append(orm.Node, project='*', with_group='groups')
        group_node_query.distinct()
        wf_nodes = group_node_query.all(flat=True)
        assert len(wf_nodes) == 1
        # This is the EOS workflow node
        node = wf_nodes[0]

        structure = node.inputs.structure
        element = structure.extras['element']
        configuration = structure.extras['configuration']

        # Get the state (possibly adding the exit status if it's finished) and add to a list
        state = node.process_state.value
        if state == 'finished':
            state = f'{state}.{node.exit_status}'
        print(f"  `-> {state}`")

        assert len(node.inputs.structure.sites) == 2, "Error, expecting 2 atoms per simulation cell!"
        num_atoms_per_cell = 2

        # Initialize to None if the outputs are not there
        eos_data = None
        stress_data = None
        BM_fit_data = None
        num_atoms = None

        min_value = None

        label_annotation = ""
        # For successfully finished workflows, fit the EOS
        if node.process_state.value == 'finished':
            if node.exit_status == 0:
                # Check if the required outputs are there
                outputs = node.get_outgoing(link_type=LinkType.RETURN).nested()
                missing_outputs = tuple(output for output in ('structures', 'total_energies') if output not in outputs)
                if missing_outputs:
                    #all_missing_outputs.append({'element': element, 'configuration': configuration, 'missing': missing_outputs})
                    warning_lines.append(f"  WARNING! MISSING OUTPUTS for workflow {node.pk}: {missing_outputs}")
                 # Extract volumes and energies for this system
                volumes = []
                energies = []
                stresses = []
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
                    parent_workflows_links = energy_node.get_incoming(link_type=LinkType.RETURN).all()
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
            else:
                #label_annotation = " (partial)"
                label_annotation = ""
                volumes = []
                energies = []
                stresses = []
                called_subprocesses = node.get_outgoing(link_type=LinkType.CALL_WORK)
                for called in called_subprocesses:
                    called = called.node
                    if issubclass(called.process_class, CommonRelaxWorkChain):
                        try:
                            volume = called.outputs.structure.get_cell_volume()
                        except AttributeError:
                            volume = called.inputs.structure.get_cell_volume()
                        try:
                            energy = called.outputs.total_energy.value
                        except AttributeError:
                            print(f"Skipping Workflow PK = {called.pk} for {smearing_type} ({degauss_meV} meV), V={volume:.3f}ang^3: no energy")
                            # Skip this calc
                            continue
                        stress = called.outputs.stress.get_array('stress').tolist()
                        volumes.append(volume)
                        energies.append(energy)
                        stresses.append(stress)


            energies = [e for _, e in sorted(zip(volumes, energies))]
            volumes = sorted(volumes)
            # List as I need to JSON-serialize it
            eos_data = (np.array([volumes, energies]).T).tolist()
            #stress_data = list(zip(volumes, stresses))
            print("  ", len(eos_data))
            if ONLY_FULLY_COMPLETE and label_annotation:
                continue
        data.append({
            'volumes': volumes,
            'energies': energies,
            'smearing_type':smearing_type,
            'degauss_meV': degauss_meV
        })

    if warning_lines:
        print("WARNINGS:")
        for line in warning_lines:
            print(line)

    with open('Er-diamond-EOS-data.json', 'w') as fhandle:
        json.dump(data, fhandle)