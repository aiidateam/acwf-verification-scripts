#!/usr/bin/env python
"""A script to create the 'set2' of structures.

This is obtained starting from 'set1', but changing all lattice constants to
improved optimized values.
"""
import os
import shutil
import sys

if os.path.exists("cifs-set2"):
    shutil.rmtree("cifs-set2")
os.mkdir("cifs-set2")

lattice_constants = {}

with open ("set2-lattice-constant.txt") as fileopen:
    for index, line in enumerate(fileopen.readlines()):
        name_to_modify, latt_constant = line.split(";")
        lattice_constants[name_to_modify] = float(latt_constant)

for name_to_modify, latt_constant in lattice_constants.items():
    cif_basename = f'POSCAR_{name_to_modify.replace("-","_")}.cif'
    with open("cifs-set1/" + cif_basename) as fhandle:
        lines = fhandle.readlines()
    with open("cifs-set2/" + cif_basename, "w") as file_to_write:        
        save_key = []
        save_value = []
        for index, line in enumerate(lines):
            if "_cell_length" in line:
                save_key.append(line.split()[0])
                save_value.append(float(line.split()[1]))
                file_to_write.write(line.split()[0] + "               " + str(latt_constant) + "\n")
            else:
                file_to_write.write(line)

        if sorted(save_key) != ["_cell_length_a", "_cell_length_b", "_cell_length_c"]:
            raise ValueError(f"{cif_basename}: Missing one lattice constant in original file!")

        if not all(abs(this_value - save_value[0]) < 1.e-12 for this_value in save_value):
            raise ValueError(f"{cif_basename} Lattice vectors have different length in original file!")

        # Check if lattice is over 5% different
        rel_difference = (latt_constant-save_value[0])/latt_constant
        if abs(rel_difference) > 0.05:
            print((f"WARNING: {cif_basename}: Very different lattice constant ({rel_difference * 100:.1f}%).\n     Original: {save_value[0]}, new: {latt_constant}"), file=sys.stderr)

