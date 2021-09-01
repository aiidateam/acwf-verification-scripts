#!/usr/bin/env python
import os
import sys

def get_plugin_name():
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

if __name__ == "__main__":
    print(f"Creating LaTeX file for the PNGs of plugin '{PLUGIN_NAME}',")
    try:
        arg = sys.argv[1]
        if arg == "many":
            folder_string = "many"
            print(f"using plots where the plugin is compared with many others")
        else:
            folder_string = f"{PLUGIN_NAME}-vs-{sys.argv[1]}"
            print(f"using plots where the plugin is compared with '{sys.argv[1]}'.")
    except IndexError:
        folder_string = PLUGIN_NAME
        print(f"using plots where the plugin is not compared with any other plugin.")
        print(f"NOTE: If you want to use the plot comparing with another plugin, pass the other")
        print(f"      plugin name as a command-line parameter.")

    if not os.path.exists(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)), os.pardir,
             "outputs", "plots-%s" % (folder_string))):
        print("ERROR! No folder ../outputs/plots-%s found." % (folder_string))
        print("       Did you run the `../outputs/generate_plots.py` script?")
        sys.exit(1)


    with open('tex-template/results.tex', 'w') as fhandle:
        fhandle.write(r"""\documentclass{article}
        \usepackage[top=0.5cm, bottom=1.5cm, left=1.5cm, right=0.5cm, landscape]{geometry}

        \usepackage{graphicx}

        \begin{document}

        """)

        for element in element_list:
            for configuration in ['XO', 'XO2', 'XO3', 'X2O', 'X2O3', 'X2O5']:
                fhandle.write("\\IfFileExists{../../outputs/plots-%s/%s-%s.png}"  % (folder_string, element, configuration))
                fhandle.write("{\\includegraphics[width=0.15\\linewidth]{../../outputs/plots-%s/%s-%s}}" % (folder_string, element, configuration))
                fhandle.write("{\\includegraphics[width=0.15\\linewidth]{missing}}\n")
            fhandle.write("\n")

        fhandle.write(r"\end{document}")
        fhandle.write("\n")

    print()
    print("-"*72)
    print("LaTeX file generated. Now you can:")
    print("  1. enter the `tex-template` subfolder;")
    print("  2. run `pdflatex results.tex`.")
    print("  3. inspect the `results.pdf` file.")
    print("(Don't forget to rename the file to include the plugin name(s) if you are going to share it).")
