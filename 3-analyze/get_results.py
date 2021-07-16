from collections import Counter
from eos_utils.eosfit_31_adapted import BM, echarge
import os

import numpy as np
import pylab as pl
import tqdm

PLOT = True

def get_plugin_name():
    import os
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

STRUCTURES_GROUP_LABEL = f'commonwf-oxides/set1/structures/{PLUGIN_NAME}'
WORKFLOWS_GROUP_LABEL = f'commonwf-oxides/set1/workflows/{PLUGIN_NAME}'

from aiida import orm
from aiida.common import LinkType

from aiida.common import LinkType

def validate_element(element_name):
    from ase.data import atomic_numbers
    assert element_name in atomic_numbers.keys(), f"Invalid element name {element_name}"

def get_conf_nice(configuration_string):
    ret_pieces = []
    for char in configuration_string:
        if char in "0123456789":
            ret_pieces.append(f"$_{char}$")
        else:
            ret_pieces.append(char)
    return "".join(ret_pieces)

def get_cottenier_data():
    all_data = {}

    # volume_factor is there because we compute per unit-cell, while the data
    # from Stefaan is per atom
    for configuration, volume_factor in [('X2O3', 10), ('X2O5', 14), ('X2O', 3), ('XO2', 3), ('XO3', 4), ('XO', 2)]:
        with open(os.path.join(os.path.dirname(__file__), 'cottenier_data_set1', f'{configuration}_temps.dat')) as fhandle:
            lines = fhandle.readlines()

        for line in lines:
            pieces = line.split()
            # I only take the first columns, more columns are typically comments
            element = pieces[0]
            validate_element(element)
            V0 = float(pieces[1]) * volume_factor
            B0_GPa = float(pieces[2])
            B0prime = float(pieces[3])

            all_data[(element, configuration)] = (V0, B0_GPa, B0prime)

    return all_data

def birch_murnaghan(V,E0,V0,B0,B01):
    r = (V0/V)**(2./3.)
    return (E0 +
            9./16. * B0 * V0 * (
            (r-1.)**3 * B01 + 
            (r-1.)**2 * (6. - 4.* r)))
   

if __name__ == "__main__":
    if PLOT:
        PLOT_FOLDER = f'outputs/plots-{PLUGIN_NAME}'
        os.makedirs(PLOT_FOLDER, exist_ok=True)

    group_node_query = orm.QueryBuilder().append(
        orm.Group, filters={'label': WORKFLOWS_GROUP_LABEL}, tag='groups',
    ).append(orm.Node, project='*', with_group='groups')
    group_node_query.distinct()
    wf_nodes = group_node_query.all(flat=True)

    states = []
    completely_off = []
    good_cnt = 0
    bad_pks = []
    cottenier_data = get_cottenier_data()
    data_to_print = {}
    output_lines = []

    progress_bar = tqdm.tqdm(wf_nodes)
    for node in progress_bar:
        structure = node.inputs.structure
        description = f"{structure.extras['element']} {structure.extras['configuration']} ({structure.pk})"
        # Using 16s to minimize length change of the description
        progress_bar.set_description(f"{description:16s}")
        progress_bar.refresh()

        state = node.process_state.value
        if state == 'finished':
            state = f'{state}.{node.exit_status}'
        states.append(state)
        if node.process_state.value == 'finished' and node.exit_status == 0:
            energies = {}
            volumes = {}
            pks = {}
            outputs = node.get_outgoing(link_type=LinkType.RETURN).nested()
            missing_outputs = tuple(output for output in ('structures', 'total_energies') if output not in outputs)
            if missing_outputs:
                output_lines.append(f"  WARNING! MISSING OUTPUTS: {missing_outputs}")
            else:
                volumes = []
                energies = []
                for index, sub_structure in sorted(outputs['structures'].items()):
                    volumes.append(sub_structure.get_cell_volume())
                    energies.append(outputs['total_energies'][index].value)
                energies = [e for _, e in sorted(zip(volumes, energies))]
                volumes = sorted(volumes)
                min_loc = np.array(energies).argmin()
                if min_loc == 0 or min_loc == len(energies) - 1:
                    completely_off.append((structure.extras['element'], structure.extras['configuration']))
                good_cnt += 1
                    
                data = np.array([volumes, energies]).T
                try:
                    min_volume, E0, bulk_modulus_internal, bulk_deriv, residuals = BM(data)
                    bulk_modulus_GPa = bulk_modulus_internal * echarge * 1.0e21
                    #1 eV/Angstrom3 = 160.21766208 GPa
                    bulk_modulus_ev_ang3 = bulk_modulus_GPa / 160.21766208
                    data_to_print[(structure.extras['element'], structure.extras['configuration'])] = (
                        min_volume, E0, bulk_modulus_GPa, bulk_deriv)
                    if residuals[0] > 1.e-6:
                        output_lines.append(f"WARNING! High fit residuals: {residuals[0]} for {structure.extras['element']} {structure.extras['configuration']}")
                    fit_ok = True
                except ValueError:
                    # If we cannot find a minimum
                    output_lines.append(f"WARNING! Unable to fit for {structure.extras['element']} {structure.extras['configuration']}")
                    fit_ok = False

                cottenier_data_this_material = cottenier_data[(structure.extras['element'], structure.extras['configuration'])]

                dense_volumes = np.linspace(min(volumes), max(volumes), 100)
                if fit_ok:
                    eos_fit_energy = birch_murnaghan(dense_volumes, E0, min_volume, bulk_modulus_ev_ang3, bulk_deriv)
                    # Note: Use the same E0 for the WIEN2K, this is how the Delta is supposed to be computed
                    eos_fit_energy_WIEN2K = birch_murnaghan(
                        dense_volumes, E0, 
                        cottenier_data_this_material[0], # min_volume
                        cottenier_data_this_material[1] / 160.21766208, # bulk_modulus_ev_ang3
                        cottenier_data_this_material[2] # bulk_deriv
                    )

                if PLOT:
                    fig = pl.figure()
                    if fit_ok:
                        pl.plot(volumes, energies - E0, 'ob', label=f'{PLUGIN_NAME}')
                    else:
                        # No -E0 since it's not computed; moreover I don't plot even Wien2K
                        # since I wouldn't know how to align it
                        pl.plot(volumes, energies, 'ob', label='QE')
                    
                    if fit_ok:
                        pl.plot(dense_volumes, eos_fit_energy_WIEN2K - E0, '-r', label='WIEN2K data')
                        pl.plot(dense_volumes, eos_fit_energy - E0, '-b', label=f'{PLUGIN_NAME} fit')
                        pl.fill_between(dense_volumes, eos_fit_energy - E0, eos_fit_energy_WIEN2K - E0, alpha=0.5, color='red')
                    
                    pl.legend(loc='upper center')
                    pl.xlabel("Cell volume ($\\AA^2$)")
                    if fit_ok:
                        pl.ylabel("$E_{tot} - E_0$ (eV)")
                    else:
                        pl.ylabel("$E_{tot} - E_0$ (eV)")
                    conf_nice = get_conf_nice(structure.extras['configuration'])
                    if fit_ok:
                        pl.title(f"{structure.extras['element']} ({conf_nice})")
                    else:
                        pl.title(f"{structure.extras['element']} ({conf_nice}) [FIT FAILED, PLOTTING ABSOLUTE ENERGIES]")
                    pl.savefig(f"{PLOT_FOLDER}/{structure.extras['element']}-{structure.extras['configuration']}.png")
                    pl.close(fig)

        elif (node.process_state.value == 'finished' and node.exit_status != 0) or (node.process_state.value == 'excepted'):
            bad_pks.append(node.pk)

    os.makedirs('outputs', exist_ok=True)
    fitted_params_fname = f"outputs/fitted_parameters-{PLUGIN_NAME}.dat"
    with open(fitted_params_fname, "w") as f:
        f.write(f'# Element Configuration V0 E0 B0 B0prime\n')
        for element in sorted(data_to_print):
            data = data_to_print[element]
            f.write(f'{element[0]} {element[1]} {data[0]} {data[1]} {data[2]} {data[3]}\n')
    print(f"Fitted parameters written to: '{fitted_params_fname}'")
    print(f"Plots written to: '{PLOT_FOLDER}'")
    
    output_lines.append("")
    output_lines.append("Counter of states: " + str(Counter(states)))

    output_lines.append("")
    output_lines.append(f"Minimum completely off for {len(completely_off)}/{good_cnt}")
    output_lines.append("Completely off systems:")
    for system in completely_off:
        output_lines.append(f"- {system[0]} {system[1]}")
    
    os.makedirs('outputs', exist_ok=True)
    fname = f"outputs/results-{PLUGIN_NAME}.txt"
    with open(fname, 'w') as fhandle:
        for line in output_lines:
            fhandle.write(f"{line}\n")
    print(f"'{fname}' written.")
