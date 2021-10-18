# Import the structures

Import the full set of structures for your set (the older ones are for historial reasons, use
only the most recent set, currently `set2`).

E.g. use `verdi archive import commonwf-oxides_set2_structures.aiida`.

The nodes are imported in a group `commonwf-oxides/set2/structures`.
Each node has extras including `Z` (atomic number), `element` (chemical element)
and `configuration` (a string identifying the anonymous formula like `XO`, `X2O3`, ...).

These can be checked e.g. with:
```python
for node in Group.get(label='commonwf-oxides/set2/structures').nodes:
    print(node.extras)
```

For reference, a CSV file with the element name, the configuration (`XO`, `XO2`, `X2O3`, ...) and
the UUID of the structure is also provided, for convenience, in the file `set2_structures_uuids.csv`.

# Create the actual (sub)group of structures that you want to run
It is highly probable that your code might not be able to deal with *all* systems (e.g.
you might not have the pseudos for elements with high Z; or you might first want to run
for a very small set).

As a first step, in the main directory of the repository (i.e., one level above the folder
named `1-preliminary`), create a `plugin_name.txt` file, containing on a single line the
exact plugin name (referred as `PLUGIN_NAME` in the following) of your plugin (e.g., `quantum_espresso` for Quantum ESPRESSO).

Then, adapt the `create_starting_subgroup.py`, addign the logic to set the `valid_elements` list for your specific plugin (add an `elif` statement for your plugin where marked by the comments). Feel free to push back your changes if you want (in that case, make sure that all imports that only your code might need are within the `if` statement, so other people can still run the script).

Then, run the script to create the appropriate subgroup `commonwf-oxides/set2/structures/PLUGIN_NAME`.
