import json
import os
import sys
import shutil
import numpy as np
import pylab as pl

if os.path.exists("cifs-verification-PBE-v1"):
    shutil.rmtree("cifs-verification-PBE-v1")
os.mkdir("cifs-verification-PBE-v1")


with open('verification-PBE-v1_alat_and_vol.json') as fhandle:
    NEW_alat = json.load(fhandle)

for name_to_modify, data in NEW_alat.items():
    #print(data["new_latt_constant"])
    element, configuration = name_to_modify.split('-')
    
    #Oxygen oxides, taking structure from F
    if element == "O":
        print(name_to_modify)
        cif_basename = f'POSCAR_F_{configuration}.cif'
        cif_newname = f'POSCAR_O_{configuration}.cif'
        with open("cifs-set2/" + cif_basename) as fhandle:
            lines = fhandle.readlines()
        with open("cifs-verification-PBE-v1/" + cif_newname, "w") as file_to_write:
            #save_key = []
            #save_value = []
            for index, line in enumerate(lines):
                if "F" in line:
                    new_line = line.replace("F","O")
                    file_to_write.write(new_line)
                elif "_cell_length" in line:
                    #save_key.append(line.split()[0])
                    #save_value.append(float(line.split()[1]))
                    file_to_write.write(line.split()[0] + "               " + str(data["new_latt_constant"]) + "\n")
                else:
                    file_to_write.write(line)

    else:
        cif_basename = f'POSCAR_{name_to_modify.replace("-","_")}.cif'
        with open("cifs-set2/" + cif_basename) as fhandle:
            lines = fhandle.readlines()
        with open("cifs-verification-PBE-v1/" + cif_basename, "w") as file_to_write:
            save_key = []
            save_value = []
            for index, line in enumerate(lines):
                if "_cell_length" in line:
                    save_key.append(line.split()[0])
                    save_value.append(float(line.split()[1]))
                    file_to_write.write(line.split()[0] + "               " + str(data["new_latt_constant"]) + "\n")
                else:
                    file_to_write.write(line)

        if sorted(save_key) != ["_cell_length_a", "_cell_length_b", "_cell_length_c"]:
            raise ValueError(f"{cif_basename}: Missing one lattice constant in original file!")

        if not all(abs(this_value - save_value[0]) < 1.e-12 for this_value in save_value):
            raise ValueError(f"{cif_basename} Lattice vectors have different length in original file!")

        # Check if lattice is over 5% different
        diff = data["new_latt_constant"]-save_value[0]
        av = data["new_latt_constant"]+save_value[0]/2
        rel_difference = diff/av
        if abs(rel_difference) > 0.01:
            print((f"WARNING: {cif_basename}: Very different lattice constant ({rel_difference * 100:.1f}%).\n     Original: {save_value[0]}, new: {data['new_latt_constant']}"), file=sys.stderr)

