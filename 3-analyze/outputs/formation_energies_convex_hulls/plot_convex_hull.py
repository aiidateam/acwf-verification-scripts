#!/usr/bin/env python
import json
import os
import pylab as pl

def get_plugin_name():
    file_name = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        os.pardir, os.pardir, os.pardir, 'plugin_name.txt'
    )
    try:
        with open(file_name) as fhandle:
            plugin_name = fhandle.read().strip()
            # Simple check e.g. to make sure there are no weird characters,
            # newlines, ... - one might still make a typo, but at least we
            # do a basic check
            assert plugin_name.isidentifier()
        return plugin_name
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "You need to define a file `../../plugin_name.txt`, containing the "
            "name of your plugin (siesta, quantum_espresso, ...) in the format "
            "expected by the aiida-common-workflows project"
        ) from exc

PLUGIN_NAME = get_plugin_name()


def draw_segment(x1, y1, x2, y2, is_hull):
    pl.plot([x1, x2], [y1, y2], "-b" if is_hull else "--g")


if __name__ == "__main__":
    # Plotting now
    OUT_FOLDER = 'convex-hulls'
    fname = f'{OUT_FOLDER}/convex-hull-data-{PLUGIN_NAME}.json'
    with open(fname) as fhandle:
        json_data = json.load(fhandle)

    # Check data version
    assert json_data['data_version'] == "0.0.1"
    all_data_flat = json_data['all_data_flat']
    all_convex_hull = json_data['all_convex_hull']

    for element in sorted(all_data_flat):
        data = all_data_flat[element]
        convex_hull = all_convex_hull[element]

        fig = pl.figure()

        for percentage_oxide, formation_energy, formula_ascii in data:
            pl.plot([percentage_oxide], [formation_energy], 'ob' if formula_ascii in convex_hull['on_hull'] else 'or')
            formula = "".join(f"$_{char}$" if char in "0123456789" else char for char in formula_ascii)
            pl.text(percentage_oxide, formation_energy, formula)

        # Plot convex hull now
        convex_hull_idx = sorted(convex_hull['on_hull'].values())
        for left_idx, right_idx in zip(convex_hull_idx[:-1], convex_hull_idx[1:]):
            draw_segment(
                data[left_idx][0], data[left_idx][1],
                data[right_idx][0], data[right_idx][1],
                is_hull=True
            )

        pl.xlim([0, 1])
        pl.xlabel('Oxygen fraction')
        pl.ylabel('eV/atom')
        pl.title(f"Convex hull ({element}-O) - {PLUGIN_NAME}\nDISCLAIMER: no correction applied!")

        fname = f'{OUT_FOLDER}/convex-hull-{PLUGIN_NAME}-{element}-O.pdf'
        pl.savefig(fname)
        print(f"File '{fname}' written.")
        pl.close(fig)
