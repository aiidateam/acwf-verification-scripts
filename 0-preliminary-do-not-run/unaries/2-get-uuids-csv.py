#!/usr/bin/env runaiida
SET_NAME = 'unaries-verification-PBE-v1'

out_lines = []
for node in Group.get(label=f'acwf-verification/{SET_NAME}/structures').nodes:
    out_lines.append(f"{node.extras['element']},{node.extras['configuration']},{node.uuid}")
out_lines.sort()

fname = f"{SET_NAME}_structures_uuids.csv"

with open(fname, "w") as fhandle:
    fhandle.write("#element,configuration,UUID\n")
    for line in out_lines:
        fhandle.write(f"{line}\n")

print(f"File '{fname}' written.")


