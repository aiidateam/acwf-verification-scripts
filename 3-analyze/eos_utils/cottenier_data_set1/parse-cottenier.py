#!/usr/bin/env python
import json
import os

def validate_element(element_name):
    from ase.data import atomic_numbers
    assert element_name in atomic_numbers.keys(), f"Invalid element name {element_name}"


def get_cottenier_data():
    all_data = {}

    # volume_factor is there because we compute per unit-cell, while the data
    # from Stefaan is per atom
    for configuration, volume_factor in [('X2O3', 10), ('X2O5', 14), ('X2O', 3), ('XO2', 3), ('XO3', 4), ('XO', 2)]:
        with open(os.path.join(os.path.dirname(__file__),  f'{configuration}_temps.dat')) as fhandle:
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

if __name__ == "__main__":
    cottenier_data = get_cottenier_data()
    data = {
        'BM_fit_data': {}
    }
    for (element, configuration), fit_data in cottenier_data.items():
        data['BM_fit_data'][f"{element}-{configuration}"] = {
            'min_volume': fit_data[0],
            'E0': None,
            'bulk_modulus_ev_ang3': fit_data[1] / 160.21766208, # convert from GPa to eV/angstrom^3,
            'bulk_deriv': fit_data[2],
            'residuals': None
        }
    with open('results-cottenier-wien2k.json', 'w') as fhandle:
        json.dump(data, fhandle, indent=2, sort_keys=True)
