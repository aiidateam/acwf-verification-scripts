import os
import shutil
import sys

if os.path.exists("cifs-set2"):
    shutil.rmtree("cifs-set2")
os.mkdir("cifs-set2")

fileopen = open("set2-lattice-constant.txt","r")
lines = fileopen.readlines()
for index, line in enumerate(lines):
    latt_constant = line.split(";")[1]
    name_to_modify = line.split(";")[0]
    name = "POSCAR_"+name_to_modify.replace("-","_")+".cif"
    str_file = open("cifs-set1/"+name,"r")
    file_to_write = open("cifs-set2/"+name,"w")
    other_lines = str_file.readlines()
    save_key = []
    save_value = []
    for index2, line2 in enumerate(other_lines):
        if "_cell_length" in line2:
            save_key.append(line2.split()[0])
            save_value.append(line2.split()[1])
            file_to_write.write(line2.split()[0] + "               " + latt_constant)
        else:
            file_to_write.write(line2)
    if len(save_key) != 3:
        print((f"{name} No three lattice constant in original file!"), file=sys.stderr)
    if "_cell_length_a" not in save_key or "_cell_length_b" not in save_key or "_cell_length_c" not in save_key:
        print((f"{name} Missing one lattice constant in original file!"), file=sys.stderr)
    if not all(element == save_value[0] for element in save_value):
        print((f"{name} Structure not cubic in original file!"), file=sys.stderr)
    if abs(float(latt_constant)-float(save_value[0]))/float(latt_constant)*100 > 5:
        print((f"{name} Very different lattice constant. Original: {save_value[0]}, new: {latt_constant.strip()}"), file=sys.stderr)
    file_to_write.close()
    str_file.close()



