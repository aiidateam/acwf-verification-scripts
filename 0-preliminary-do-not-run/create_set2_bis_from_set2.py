#!/usr/bin/env python
"""A script to create the 'set2-bis' of structures.

This an intermediate set to generate new cental volumes for 7 structures:
 Ar-XO2   51.35417444630869  4.17218524915 
 He-X2O   91.44084234877081  5.0569082255012 
 He-XO    31.34727051532348  3.5391984073003 
 Kr-XO2   47.54670488377530  4.0664160182532
 Ne-X2O   95.93779745747102  5.1384828791316 
 Ne-XO    47.08525282205443  4.0532180515519
 Rb-XO3  113.24793199207410  4.8381213822293
and to add the oxigen oxides:
 O-X2O3   95.726216868625   4.5745
 O-X2O5  116.057167806952   4.8778
 O-X2O    27.095308         3.3713437
 O-XO2    27.127714         3.3726872
 O-XO3    44.701078149000   3.5490
 O-XO     15.910169         2.8231238
"""
import os
import shutil
import sys

if os.path.exists("cifs-set2-bis"):
    shutil.rmtree("cifs-set2-bis")
os.mkdir("cifs-set2-bis")

lattice_constants = {}

new_a0 = {
        "Ar-XO2":  4.17218524915, 
        "He-X2O":  5.0569082255012,
        "He-XO" :  3.5391984073003,
        "Kr-XO2":  4.0664160182532,
        "Ne-X2O":  5.1384828791316,
        "Ne-XO":   4.0532180515519,
        "Rb-XO3":  4.8381213822293
        }

oxigen_oxides = {
        "X2O3":  4.5745,
        "X2O5":  4.8778,
        "X2O":   3.3713437,
        "XO2":   3.3726872,
        "XO3":   3.5490,
        "XO":    2.8231238
    }


for name_to_modify, latt_constant in new_a0.items():
    cif_basename = f'POSCAR_{name_to_modify.replace("-","_")}.cif'
    with open("cifs-set2/" + cif_basename) as fhandle:
        lines = fhandle.readlines()
    with open("cifs-set2-bis/" + cif_basename, "w") as file_to_write:        
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
        print(name_to_modify, save_value[0], latt_constant)
        rel_difference = (latt_constant-save_value[0])/latt_constant
        if abs(rel_difference) > 0.05:
            print((f"WARNING: {cif_basename}: Very different lattice constant ({rel_difference * 100:.1f}%).\n     Original: {save_value[0]}, new: {latt_constant}"), file=sys.stderr)

for conf, latt_constant in oxigen_oxides.items():
    cif_basename = f'POSCAR_F_{conf}.cif'
    cif_newname = f'POSCAR_O_{conf}.cif'
    with open("cifs-set2/" + cif_basename) as fhandle:
        lines = fhandle.readlines()
    with open("cifs-set2-bis/" + cif_newname, "w") as file_to_write:
        save_key = []
        save_value = []
        for index, line in enumerate(lines):
            if "F" in line:
                new_line = line.replace("F","O")
                file_to_write.write(new_line)
            elif "_cell_length" in line:
                save_key.append(line.split()[0])
                save_value.append(float(line.split()[1]))
                file_to_write.write(line.split()[0] + "               " + str(latt_constant) + "\n")
            else:
                file_to_write.write(line)

