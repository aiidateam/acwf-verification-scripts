import json
from ase.data import chemical_symbols
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
from pymatgen.core import Structure


Z_max = 103

UNARIES_OXIDES = 'oxides'

CONFIGURATIONS = {
    'unaries': ['SC', 'BCC', 'FCC', 'Diamond'],
    'oxides': ['X2O', 'X2O5', 'XO2', 'X2O3', 'XO', 'XO3']
}

PRIM_CONV_LATTICE_CONST_SCALING = {
    'unaries': {
        'SC': 1,
        'BCC': 2 / np.sqrt(3),
        'FCC': np.sqrt(2),
        'Diamond': np.sqrt(2)
    },
    'oxides': {
        'X2O': np.sqrt(2),
        'X2O5': 1,
        'XO2': np.sqrt(2),
        'X2O3': 1,
        'XO': np.sqrt(2),
        'XO3': 1
    }
}



def plot_lattice_constants(lattice_constants):
    nrows = len(CONFIGURATIONS[UNARIES_OXIDES])
    fig, ax = plt.subplots(nrows, 1, figsize=(12, 4 * nrows))
    plt.subplots_adjust(hspace=0.8)

    for i, crystal_lattice in enumerate(CONFIGURATIONS[UNARIES_OXIDES]):
        for _, xc in enumerate(['PBE', 'PBEsol', 'LDA']):
            x = np.arange(1, Z_max + 1)
            y = np.zeros(len(x))
            for j in range(1, Z_max + 1):
                el = chemical_symbols[j]
                y[j-1] = lattice_constants[xc][crystal_lattice].get(el, np.NaN)
            ax[i].plot(x, y, label=xc, marker='o', markersize=3, linewidth=1)
            
        ax[i].set_title(f'{crystal_lattice}')
        ax[i].set_ylabel('Lattice constant (Å)')
        ax[i].set_xlim(1, Z_max)
        ax[i].minorticks_on()
        ax[i].legend(ncols=3)

        ax[i].xaxis.set_major_locator(MultipleLocator(10))
        ax[i].minorticks_on()
        ax[i].xaxis.set_minor_locator(MultipleLocator(1))
        ax[i].grid(which='major', axis='x', color='#ccc', linestyle='-')
        ax[i].grid(which='minor', axis='x', color='#eee', linestyle='-')
        ax[i].grid(which='major', axis='y', color='#ccc', linestyle='-')
        ax[i].grid(which='minor', axis='y', color='#eee', linestyle='dotted')

        ticks_x = np.arange(1, Z_max+1)
        sec = ax[i].secondary_xaxis(location='top')
        sec.set_xticks(ticks_x)
        sec.set_xticklabels([chemical_symbols[Z] for Z in ticks_x], fontsize=7)
        sec.tick_params(rotation=90)
        plt.savefig(f'lattice_constants_{UNARIES_OXIDES}.png', bbox_inches='tight')
        
def check_lattice_constants(lattice_constants):
    functional_structures_dir_mapping = {
        'PBE': f'xsfs-{UNARIES_OXIDES}-verification-PBE-missin-actinides-v1',
        'PBEsol': f'xsfs-{UNARIES_OXIDES}-verification-PBEsol-v1',
        'LDA': f'xsfs-{UNARIES_OXIDES}-verification-LDA-v1'
    }

    for xc in ['PBE', 'PBEsol', 'LDA']:
        for config in CONFIGURATIONS[UNARIES_OXIDES]:
            for Z in range(1, Z_max + 1):
                el = chemical_symbols[Z]
                try:
                    s = Structure.from_file(
                        f'{functional_structures_dir_mapping[xc]}/{el}-{config}.xsf'
                    )
                except FileNotFoundError:
                    continue
                
                assert np.isclose(s.lattice.abc[0], s.lattice.abc[1]), f"Non-cubic lattice found for {el}-{config}"
                assert np.isclose(s.lattice.abc[1], s.lattice.abc[2]), f"Non-cubic lattice found for {el}-{config}"
                assert np.isclose(
                    s.lattice.abc[0] * PRIM_CONV_LATTICE_CONST_SCALING[UNARIES_OXIDES][config], 
                    lattice_constants[xc][config].get(el)
                    ), f"Non-matching lattice constant found for {config} {el} in {xc} calculation, expected {s.lattice.abc[0] * PRIM_CONV_LATTICE_CONST_SCALING[UNARIES_OXIDES][config]} but got {lattice_constants[xc][config].get(el)}"

if __name__ == "__main__":
    prefix = 'lattice_parameters_' if UNARIES_OXIDES == 'unaries' else 'conv_lattice_constants_'
    
    with open(f'{prefix}{UNARIES_OXIDES}-verification-LDA-v1.json', 'r') as f:
        lda_alat = json.load(f)
        
    with open(f'{prefix}{UNARIES_OXIDES}-verification-PBE-missing-actinides-v1.json', 'r') as f:
        pbe_alat = json.load(f)

    with open(f'{prefix}{UNARIES_OXIDES}-verification-PBEsol-v1.json', 'r') as f:
        pbesol_alat = json.load(f)
        
    with open(f'{prefix}{UNARIES_OXIDES}-verification-PBE-v1.json', 'r') as f:
        pbe_1_96 = json.load(f)

        for k in pbe_1_96:
            # This is the old one that was calculated by accident
            pbe_1_96[k].pop('Cm')

        for k in pbe_1_96:
            pbe_alat[k].update(pbe_1_96[k])
        
    lattice_constants = {
        'LDA': lda_alat,
        'PBE': pbe_alat,
        'PBEsol': pbesol_alat
    }

    plot_lattice_constants(lattice_constants)
    
    check_lattice_constants(lattice_constants)