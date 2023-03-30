"""
Script to convert thec cifs in cifs-oxides-verification-PBE-v1 into xsf. We use pymatgen.
"""
import os
import pymatgen, pymatgen.core

subfolder = 'xsfs-oxides-verification-PBE-v1'
os.mkdir(subfolder)

for i in os.listdir('cifs-oxides-verification-PBE-v1'):
    structure = pymatgen.core.Structure.from_file(f'cifs-oxides-verification-PBE-v1/{i}')
    name = i.replace('POSCAR_', '')
    name = name.replace('_','-')
    name = name.replace('.cif','.xsf')
    structure.to(filename=os.path.join(subfolder, name))

