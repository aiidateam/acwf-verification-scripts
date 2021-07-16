#!/usr/bin/env python

def get_plugin_name():
    import os
    file_name = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        os.pardir, os.pardir, 'plugin_name.txt'
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


## ADAPT TO YOUR PLUGIN NAME
PLUGIN_NAME = get_plugin_name()

## ADAPT TO THE CORRECT LIST OF ELEMENTS
element_list = [
    'Ag', 'Al', 'Ar', 'As', 'Au', 'B', 'Ba', 'Be', 'Bi', 'Br', 'C', 'Ca', 'Cd', 
    'Ce', 'Cl', 'Co', 'Cr', 'Cs', 'Cu', 'Dy', 'Er', 'Eu', 'F', 'Fe', 'Ga', 'Gd',
    'Ge', 'H', 'He', 'Hf', 'Hg', 'Ho', 'I', 'In', 'Ir', 'K', 'Kr', 'La', 'Li',
    'Lu', 'Mg', 'Mn', 'Mo', 'N', 'Na', 'Nb', 'Nd', 'Ne', 'Ni', 'Os', 'P', 'Pb',
    'Pd', 'Pm', 'Po', 'Pr', 'Pt', 'Rb', 'Re', 'Rh', 'Rn', 'Ru', 'S', 'Sb', 'Sc',
    'Se', 'Si', 'Sm', 'Sn', 'Sr', 'Ta', 'Tb', 'Tc', 'Te', 'Ti', 'Tl', 'Tm', 'V',
    'W', 'Xe', 'Y', 'Yb', 'Zn', 'Zr']

with open('tex-template/results.tex', 'w') as fhandle:
    fhandle.write(r"""\documentclass{article}
    \usepackage[top=0.5cm, bottom=1.5cm, left=1.5cm, right=0.5cm, landscape]{geometry}

    \usepackage{graphicx}

    \begin{document}

    """)

    for element in element_list:
        for configuration in ['XO', 'XO2', 'XO3', 'X2O', 'X2O3', 'X2O5']:
            fhandle.write("\\IfFileExists{../../outputs/plots-%s/%s-%s.png}"  % (PLUGIN_NAME, element, configuration))
            fhandle.write("{\\includegraphics[width=0.15\\linewidth]{../../outputs/plots-%s/%s-%s}}" % (PLUGIN_NAME, element, configuration))
            fhandle.write("{\\includegraphics[width=0.15\\linewidth]{missing}}\n")
        fhandle.write("\n")

    fhandle.write(r"\end{document}")
    fhandle.write("\n")
