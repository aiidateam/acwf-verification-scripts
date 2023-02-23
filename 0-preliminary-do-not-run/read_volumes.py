import pymatgen
import ase.io
from ase.data import chemical_symbols


def get_num_atoms_in_formula_unit(configuration):
    """Return the number ot atoms per formula unit.

    For unaries, it's the number in the primitive cell: 1 for SC, BCC, FCC; 2 for diamond.
    """
    mapping = {
        'XO': 2,
        'XO2': 3,
        'X2O': 3,
        'X2O3': 5,
        'XO3': 4,
        'X2O5': 7,
        'X/Diamond': 2,
        'X/SC': 1,
        'X/FCC': 1,
        'X/BCC': 1
    }

    try:
        return mapping[configuration]
    except KeyError:
        raise ValueError(f"Unknown number of atoms in formula unit for configuration '{configuration}'")


def get_table():
    print('  & FCC & BCC & SC & Diamond & X2O & X2O5 & XO2 & X2O3 & XO & XO3 \\\\')
    print('\hline')
    SET_NAME = 'oxides-verification-PBE-v1'
    for Z in range(1, 96+1):
        dict_vols = {}
        element_name = chemical_symbols[Z]
        for configuration in ['X2O', 'X2O5', 'XO2', 'X2O3', 'XO', 'XO3']:
            filename = f'oxides/cifs-{SET_NAME}/POSCAR_{element_name}_{configuration}.cif'
            pmg_structure = pymatgen.Structure.from_file(filename)
            number_units = len(pmg_structure.sites)/get_num_atoms_in_formula_unit(configuration)
            volume = pmg_structure.volume/number_units
            dict_vols[configuration] = volume
        for configuration in ['FCC', 'BCC', 'SC', 'Diamond']:
            filename = f'unaries/xsfs-unaries-verification-PBE-v1/{element_name}-{configuration}.xsf'
            stru_ase = ase.io.read(filename)
            number_units = len(stru_ase)/get_num_atoms_in_formula_unit(f'X/{configuration}')
            volume = stru_ase.get_volume()/number_units
            dict_vols[configuration] = volume
        string = f'{element_name}'
        for i in ['FCC', 'BCC', 'SC', 'Diamond', 'X2O', 'X2O5', 'XO2', 'X2O3', 'XO', 'XO3']:
            string = string + f' & {round(dict_vols[i],5)}'
        string = string + ' \\\\'
        print(string)


if __name__ == "__main__":
    get_table()
