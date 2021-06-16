#!/usr/bin/env python

## ADAPT TO YOUR CODE
CODE_NAME = "quantum_espresso"

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
            fhandle.write("\\IfFileExists{../../outputs/plots-%s/%s-%s.png}"  % (CODE_NAME, element, configuration))
            fhandle.write("{\\includegraphics[width=0.15\\linewidth]{../../outputs/plots-%s/%s-%s}}" % (CODE_NAME, element, configuration))
            fhandle.write("{\\includegraphics[width=0.15\\linewidth]{missing}}\n")
        fhandle.write("\n")

    fhandle.write(r"\end{document}")
    fhandle.write("\n")
