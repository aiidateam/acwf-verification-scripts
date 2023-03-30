"""
Script to convert thec cifs in cifs-oxides-verification-PBE-v1 into xsf. We use pymatgen.
"""
import os
from pymatgen.io.cif import CifParser

os.mkdir('xsfs-oxides-verification-PBE-v1')

for i in os.listdir('cifs-oxides-verification-PBE-v1'):
    cif = CifParser(f'cifs-oxides-verification-PBE-v1/{i}')
    structure = cif.get_structures()[0]
    name = i.replace('POSCAR_', '')
    name = name.replace('_','-')
    name = name.replace('.cif','.xsf')
    structure.to('xsf', name)
    os.replace(name, 'xsfs-unaries-verification-PBE-v1/'+name)

