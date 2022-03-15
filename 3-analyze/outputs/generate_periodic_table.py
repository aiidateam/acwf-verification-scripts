#!/usr/bin/env python
import json
import os
import sys

import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
import tqdm

import quantities_for_comparison as qc

SHOW_IN_BROWSER=False
DEFAULT_PREFACTOR = 100
DEFAULT_wb0 = 1.0/8.0
DEFAULT_wb1 = 1.0/64.0
EXPECTED_SCRIPT_VERSION = "0.0.3"


from bokeh.models import (
    ColumnDataSource,
    LinearColorMapper,
    LogColorMapper,
    ColorBar,
    BasicTicker,
)
from bokeh.plotting import figure, output_file
from bokeh.io import show as show_, export_png
from bokeh.sampledata.periodic_table import elements
from bokeh.transform import dodge
from matplotlib.colors import Normalize, LogNorm, to_hex
from matplotlib.cm import (
    plasma,
    inferno,
    magma,
    viridis,
    cividis,
    turbo,
    ScalarMappable,
)
from pandas import options
from typing import List
import warnings
from bokeh.io import export_svg


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
            #assert plugin_name.isidentifier()
        return plugin_name
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "You need to define a file `../../plugin_name.txt`, containing the "
            "name of your plugin (siesta, quantum_espresso, ...) in the format "
            "expected by the aiida-common-workflows project"
        ) from exc


def abs_V0_rel_diff(*args, **kwargs):
    return abs(qc.V0_rel_diff(*args, **kwargs))
def abs_B0_rel_diff(*args, **kwargs):
    return abs(qc.B0_rel_diff(*args, **kwargs))
def abs_B1_rel_diff(*args, **kwargs):
    return abs(qc.B1_rel_diff(*args, **kwargs))

quantity_for_comparison_map = {
    "delta_per_formula_unit": qc.delta,
    "B0_rel_diff": qc.B0_rel_diff,
    "V0_rel_diff": qc.V0_rel_diff,
    "B1_rel_diff": qc.B1_rel_diff,
    "abs_V0_rel_diff": abs_V0_rel_diff,
    "abs_B0_rel_diff": abs_B0_rel_diff,
    "abs_B1_rel_diff": abs_B1_rel_diff,            
    "rel_errors_vec_length": qc.rel_errors_vec_length,
    "epsilon": qc.epsilon
}


PLUGIN_NAME = get_plugin_name()

if __name__ == "__main__":
    try:
        SET_NAME = sys.argv[1]
    except IndexError:
        print(f"The first argument must be the set name.")
        sys.exit(1)

    try:
        QUANTITY = sys.argv[2]
    except IndexError:
        print(f"The second argument must be the quantity to use for comparison. Choose among {quantity_for_comparison_map.keys()}")
        sys.exit(1)

    if QUANTITY not in quantity_for_comparison_map.keys():
        print(f"The second argument must be the quantity to use for comparison. Choose among {quantity_for_comparison_map.keys()}")
        sys.exit(1)

    try:
        other_code = sys.argv[3]
    except IndexError:
        print(f"The third argument must be the code to compare with")
        sys.exit(1)

    try:
        with open(f'results-{SET_NAME}-{PLUGIN_NAME}.json') as fhandle:
            reference_plugin_data = json.load(fhandle)
    except OSError:
        print(f"No data found for your plugin '{PLUGIN_NAME}' (set '{SET_NAME}'). Did you run `./get_results.py` first?")
        sys.exit(1)

    try:
        with open(f'results-{SET_NAME}-{other_code}.json') as fhandle:
            compare_plugin_data = json.load(fhandle)
    except OSError:
        print(f"No data found for code '{other_code}' (set '{SET_NAME}'). A file results-{SET_NAME}-{other_code}.json is expected")
        sys.exit(1)

    if not reference_plugin_data['script_version'] == EXPECTED_SCRIPT_VERSION:
        raise ValueError(
            f"This script only works with data generated at version {EXPECTED_SCRIPT_VERSION}. "
            f"Please re-run ./get_results.py to update the data format for {PLUGIN_NAME}!"
            )

    print(f"Using data for plugin '{PLUGIN_NAME}' (set '{SET_NAME}') compared with {other_code}.")

    #name_file = f'histo-{QUANTITY}-{SET_NAME}-{PLUGIN_NAME}.pdf'

    all_systems = set(reference_plugin_data['eos_data'].keys())
    all_systems = set(reference_plugin_data['BM_fit_data'].keys())
    #all_systems.update(compare_plugin_data['BM_fit_data'].keys())

    collect = {
        "X/Diamond" : {"elements": [], "values": []},
        "X/FCC" : {"elements": [], "values": []},
        "X/BCC" : {"elements": [], "values": []},
        "X/SC" : {"elements": [], "values": []},
        "X2O3" : {"elements": [], "values": []},
        "X2O5" : {"elements": [], "values": []},
        "XO2" : {"elements": [], "values": []},
        "XO3" : {"elements": [], "values": []},
        "XO" : {"elements": [], "values": []},
        "X2O" : {"elements": [], "values": []}
        }

    #print(f"comparing with {all_args[index]}")
    progress_bar = tqdm.tqdm(sorted(all_systems))
    for element_and_configuration in progress_bar:
        progress_bar.set_description(f"{element_and_configuration:12s}")
        progress_bar.refresh()

        element, configuration = element_and_configuration.split('-')
        # Get the data for the reference plugin
        ref_BM_fit_data = reference_plugin_data['BM_fit_data'][f'{element}-{configuration}']
    
        if ref_BM_fit_data is None:
            continue
    
        scaling_factor_ref = qc.get_volume_scaling_to_formula_unit(
                reference_plugin_data['num_atoms_in_sim_cell'][f'{element}-{configuration}'],
                element, configuration
            )

        V0=ref_BM_fit_data['min_volume']/scaling_factor_ref
        B0=ref_BM_fit_data['bulk_modulus_ev_ang3']
        B01=ref_BM_fit_data['bulk_deriv']

        # Get the data for the compare_with plugin, if specified (and if the EOS worked for the 
        # reference plugin, otherwise we don't know which E0 to use)
        try:
            compare_BM_fit_data = compare_plugin_data['BM_fit_data'][f'{element}-{configuration}']
            if compare_BM_fit_data is None:
                # No fitting data in the plugin to compare with.
                # Raise this exception that is catched one line below, so
                # it will set `compare_eos_fit_energy` to None.
                raise KeyError                    
        except KeyError:
            # Set to None if fit data is missing (if we are here, the EOS points
            # are there, so it means that the fit failed). I will still plot the
            # points
            continue

        scaling_factor_comp = qc.get_volume_scaling_to_formula_unit(
                compare_plugin_data['num_atoms_in_sim_cell'][f'{element}-{configuration}'],
                element, configuration
            )

        CV0=compare_BM_fit_data['min_volume']/scaling_factor_comp
        CB0=compare_BM_fit_data['bulk_modulus_ev_ang3']
        CB01=compare_BM_fit_data['bulk_deriv']

        quant = quantity_for_comparison_map[QUANTITY](V0,B0,B01,CV0,CB0,CB01,DEFAULT_PREFACTOR,DEFAULT_wb0,DEFAULT_wb1)

        collect[configuration]["values"].append(quant)
        collect[configuration]["elements"].append(element)

    # Way to decide whether is a unaries or an oxides set is a bit fragile.
    if collect["X/Diamond"]["values"]:
        unaries = True
        list_confs = ["X/Diamond","X/FCC","X/BCC","X/SC"]
    else:
        unaries = False
        list_confs = ["X2O3","X2O5","X2O","XO2","XO3","XO"]

    width = 1050
    cmap = "plasma"
    alpha = 0.65
    extended = True
    log_scale = False
    cbar_height = None
    cbar_standoff = 12
    cbar_fontsize = 14
    blank_color = "#c4c4c4"
    under_value = None
    under_color = "#140F0E"
    over_value = None
    over_color = "#140F0E"
    special_elements = None
    special_color = "#6F3023"

    options.mode.chained_assignment = None

    # Assign color palette based on input argument
    if cmap == "plasma":
        cmap = plasma
        bokeh_palette = "Plasma256"
    elif cmap == "magma":
        cmap = magma
        bokeh_palette = "Magma256"
    else:
        raise ValueError("Unknown color map")

    # Define number of and groups
    period_label = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    group_range = [x for x in range(1, 19)]

    #We "fake" that La-Yb has period 9 and group from 5 to 18, also that
    #Th-Lr has period 10 and group from 5 to 18. It is just to place them in
    #the correct point in the chart.
    count=0
    for i in range(56, 70):
        elements.period[i] = 9
        elements.group[i] = count + 4
        count += 1

    count = 0
    for i in range(88, 102):
        elements.period[i] = 10
        elements.group[i] = count + 4
        count += 1

    per = [int(i) for i in elements["period"]]
    grou = [int(i) for i in elements["group"]]

    min_data = min([min(collect[i]["values"]) for i in list_confs])
    max_data = max([max(collect[i]["values"]) for i in list_confs])

    color_list={}

    for conf in list_confs:
        data_elements = collect[conf]["elements"]
        data = collect[conf]["values"]

        if len(data) != len(data_elements):
            raise ValueError("Unequal number of atomic elements and data points")

        # Define matplotlib and bokeh color map
        if log_scale:
            for datum in data:
                if datum < 0:
                    raise ValueError(
                        f"Entry for element {datum} is negative but log-scale is selected"
                    )
            color_mapper = LogColorMapper(
                palette=bokeh_palette, low=min_data, high=max_data
            )
            norm = LogNorm(vmin=min_data, vmax=max_data)
        else:
            color_mapper = LinearColorMapper(
                palette=bokeh_palette, low=min_data, high=max_data
            )
            norm = Normalize(vmin=min_data, vmax=max_data)
            
        color_scale = ScalarMappable(norm=norm, cmap=cmap).to_rgba(data, alpha=None)

        # Set blank color
        color_list[conf] = [blank_color] * len(elements)


        # Compare elements in dataset with elements in periodic table
        for i, data_element in enumerate(data_elements):
            element_entry = elements.symbol[
                elements.symbol.str.lower() == data_element.lower()
            ]
            if element_entry.empty == False:
                element_index = element_entry.index[0]
            else:
                warnings.warn("Invalid chemical symbol: " + data_element)
            if color_list[conf][element_index] != blank_color:
                warnings.warn("Multiple entries for element " + data_element)
            elif under_value is not None and data[i] <= under_value:
                color_list[conf][element_index] = under_color
            elif over_value is not None and data[i] >= over_value:
                color_list[conf][element_index] = over_color
            else:
                color_list[conf][element_index] = to_hex(color_scale[i])

    if unaries:
        # Define figure properties for visualizing data
        source = ColumnDataSource(
            data=dict(
                group=grou,
                period=per,
                top=[i-0.45 for i in per],
                bottom=[i+0.45 for i in per],
                left=[i-0.45 for i in grou],
                right=[i+0.45 for i in grou],
                sym=elements["symbol"],
                atomic_number=elements["atomic number"],
                type_color_dia=color_list["X/Diamond"],
                type_color_sc=color_list["X/SC"],
                type_color_bcc=color_list["X/BCC"],
                type_color_fcc=color_list["X/FCC"],
            )
        )

        # Plot the periodic table
        p = figure(x_range=[0,19], y_range=[11,0], tools="save")
        p.plot_width = width
        p.outline_line_color = None
        p.background_fill_color = None
        p.border_fill_color = None
        p.toolbar_location = "above"
        p.quad(left="left", right="group", top="period", bottom="bottom", source=source, alpha=alpha, color="type_color_dia")
        p.quad(left="left", right="group", top="top", bottom="period", source=source, alpha=alpha, color="type_color_sc")
        p.quad(left="group", right="right", top="period", bottom="bottom", source=source, alpha=alpha, color="type_color_bcc")
        p.quad(left="group", right="right", top="top", bottom="period", source=source, alpha=alpha, color="type_color_fcc")
        p.axis.visible = False

        #The reference block
        p.quad(left=5,right=6.5,bottom=0,top=1.5,fill_color="white")
        p.quad(left=5,right=6.5,bottom=1.5,top=3,fill_color="white")
        p.quad(left=6.5,right=8,bottom=0,top=1.5,fill_color="white")
        p.quad(left=6.5,right=8,bottom=1.5,top=3,fill_color="white")
        xx=[5.75,5.75,7.25,7.25]
        yy=[0.75,2.25,0.75,2.25]
        text=["SC","Diamond","FCC","BCC"]
        sou = ColumnDataSource(dict(x=xx, y=yy, text=text))
        p.text(
            x="x",
            y="y",
            text="text",
            source=sou,
            text_font_style="bold",
            text_font_size="17pt",
            text_align= "center",
            text_baseline= "middle",
            angle=-45,
            angle_units="grad"
            )

    else:
        # Define figure properties for visualizing data
        source = ColumnDataSource(
            data=dict(
                group=grou,
                period=per,
                top=[i-0.45 for i in per],
                bottom=[i+0.45 for i in per],
                left=[i-0.45 for i in grou],
                right=[i+0.45 for i in grou],
                midup=[i-0.15 for i in per],
                middown=[i+0.15 for i in per],
                sym=elements["symbol"],
                atomic_number=elements["atomic number"],
                type_color_X2O3=color_list["X2O3"],
                type_color_X2O5=color_list["X2O5"],
                type_color_X2O=color_list["X2O"],
                type_color_XO2=color_list["XO2"],
                type_color_XO3=color_list["XO3"],
                type_color_XO=color_list["XO"],
            )
        )

        # Plot the periodic table
        p = figure(x_range=[0,19], y_range=[11,0], tools="save")
        p.plot_width = width
        p.outline_line_color = None
        p.background_fill_color = None
        p.border_fill_color = None
        p.toolbar_location = "above"
        p.quad(left="left", right="group", top="top", bottom="midup", source=source, alpha=alpha, color="type_color_X2O3")
        p.quad(left="left", right="group", top="midup", bottom="middown", source=source, alpha=alpha, color="type_color_X2O")
        p.quad(left="left", right="group", top="middown", bottom="bottom", source=source, alpha=alpha, color="type_color_XO3")
        p.quad(left="group", right="right", top="top", bottom="midup", source=source, alpha=alpha, color="type_color_X2O5")
        p.quad(left="group", right="right", top="midup", bottom="middown", source=source, alpha=alpha, color="type_color_XO2")
        p.quad(left="group", right="right", top="middown", bottom="bottom", source=source, alpha=alpha, color="type_color_XO")
        p.axis.visible = False

        #The reference block
        p.quad(left=5,right=6.5,bottom=0,top=1,fill_color="white")
        p.quad(left=5,right=6.5,bottom=1,top=2,fill_color="white")
        p.quad(left=5,right=6.5,bottom=2,top=3,fill_color="white")
        p.quad(left=6.5,right=8,bottom=0,top=1,fill_color="white")
        p.quad(left=6.5,right=8,bottom=1,top=2,fill_color="white")
        p.quad(left=6.5,right=8,bottom=2,top=3,fill_color="white")
        xx=[5.75,5.75,5.75,7.25,7.25,7.25]
        yy=[0.5,1.5,2.5,0.5,1.5,2.5]
        text=["X2O3","X2O","XO3","X2O5","XO2","XO"]
        sou = ColumnDataSource(dict(x=xx, y=yy, text=text))
        p.text(
            x="x",
            y="y",
            text="text",
            source=sou,
            text_font_style="bold",
            text_font_size="17pt",
            text_align= "center",
            text_baseline= "middle",
            )


    #Add element name
    text_props = {
        "source": source,
        "angle": 0,
        "color": "black",
        "text_align": "center",
        "text_baseline": "middle",
    }
    #x = dodge("group", -0.4, range=p.x_range)
    #y = dodge("period", 0.3, range=p.y_range)
    p.text(
        x="group",
        y="period",
        text="sym",
        text_font_style="bold",
        text_font_size="16pt",
        **text_props,
    )
    #p.text(x=x, y=y, text="atomic_number", text_font_size="11pt", **text_props)

    color_bar = ColorBar(
        color_mapper=color_mapper,
        ticker=BasicTicker(desired_num_ticks=10),
        border_line_color=None,
        label_standoff=cbar_standoff,
        location=(0, 0),
        orientation="vertical",
        scale_alpha=alpha,
        major_label_text_font_size=f"{cbar_fontsize}pt",
    )

    if cbar_height is not None:
        color_bar.height = cbar_height

    p.add_layout(color_bar, "right")
    p.grid.grid_line_color = None

     # Open in a browser
    if SHOW_IN_BROWSER:
        output_file("periodic-table-plot.html")
        show_(p)
    else:
        try:
            export_png(p, filename="periodic-table-plot.png")
        except RuntimeError as exc:
            msg = str(exc)
            msg = f"""

ERROR GENERATING THE IMAGE!
The original error message was:
{msg}

Please check the following:
- Bokeh instructions here: https://docs.bokeh.org/en/latest/docs/user_guide/export.html#additional-dependencies
- That you installed the requirements.txt file, and in particular that you installed
  `pip install selenium chromedriver-binary`
  (to use with Chrome)
- that you have a recent version of Chrome
- that you downloaded from https://chromedriver.chromium.org/ and put in your PATH
  the chromedriver executable for the *SAME* version of Chrome that you have
  (note that Chrome typically self-updates, so check even if this script
  used to work, check that now Chrome is not more recent than the chromedriver you had installed;
  in this case udpate it).
"""
            raise RuntimeError(msg)

