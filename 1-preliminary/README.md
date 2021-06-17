# Import the structures

Import the full set of structures for your set (for now there is only 1 set).

E.g. use `verdi archive import commonwf-oxides_set1_structures.aiida`.

The nodes are imported in a group `commonwf-oxides/set1/structures`.
Each node has extras including `Z` (atomic number), `element` (chemical element)
and `configuration` (a string identifying the anonymous formula like `XO`, `X2O3`, ...).

These can be checked e.g. with:
```python
for node in Group.get(label='commonwf-oxides/set1/structures').nodes:
    print(node.extras)
```

For reference, a CSV file with the element name, the configuration (`XO`, `XO2`, `X2O3`, ...) and
the UUID of the structure is also provided, for convenience, in the file `set1_structures_uuids.csv`.

# Create the actual (sub)group of structures that you want to run
It is highly probable that your code might not be able to deal with *all* systems (e.g.
you might not have the pseudos for elements with high Z; or you might first want to run
for a very small set).

Then edit the `create_starting_subgroup.py` replacing the `PLUGIN_NAME = 'quantum_espresso'`,
and changing the logic to set the `valid_elements` list.

Then, run the script to create the appropriate subgroup `commonwf-oxides/set1/structures/PLUGIN_NAME`.
