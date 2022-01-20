out_lines = []
for node in Group.get(label='commonwf-oxides/monoelemental-structures-test').nodes:
    out_lines.append(f"{node.extras['element']},{node.extras['configuration']},{node.uuid}")
out_lines.sort()

fname = "monoelemental_structures_test_uuids.csv"

with open(fname, "w") as fhandle:
    fhandle.write("#element,configuration,UUID\n")
    for line in out_lines:
        fhandle.write(f"{line}\n")

print(f"File '{fname}' written.")

    
