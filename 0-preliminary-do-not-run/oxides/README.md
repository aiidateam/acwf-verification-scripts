# Preliminary helper scripts

These scripts shouldn't be run by all users: they are supposed to be run
only once, when the structures need to be imported.
(So: go to folder 1 if you are not sure why you are here :-) )

**Note**: before running, in each of the three scripts, change the variable `SET_NAME` to the correct name of the set of structures you want to process

- `1-create-aiida-archive-from-cifs.py`: read the folder with the CIF files, convert to AiiDA StructureData and put them in a properly named group

- `2-get-uuids-csv.py`: create a CSV containing, for each material, the UUID of the corresponding AiiDA structure (StructureData)

- `3-export_structure_group.sh`: this will generate the AiiDA archive file containing the structures and the group


**Note2**: the original structure for the oxides are in cif format, we later on created a script to transform them is xsf
