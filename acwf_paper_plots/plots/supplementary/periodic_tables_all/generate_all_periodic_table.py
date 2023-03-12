#!/usr/bin/env python
import json
import os
import sys

import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
import tqdm

import acwf_paper_plots.quantities_for_comparison as qc

SHOW_IN_BROWSER=False
DEFAULT_wb0 = 1.0/20.0
DEFAULT_wb1 = 1.0/400.0
# Default prefactor if not indicated: 1.
PREFACTOR_DICT = {'nu': 100.}
EXPECTED_SCRIPT_VERSION = ["0.0.3","0.0.4"]
EXCELLENT_AGREEMENT_THRESHOLD = {
    'nu': 0.1, 'epsilon': 0.07,
    'delta_per_formula_unit': 0., # I put zero, it's not used in this script anyway
    'delta_per_formula_unit_over_b0': 0. # I put zero, it's not used in this script anyway

    }
PRINT_NON_EXCELLENT = False

## IN ORDER TO PLOT ALL
# Whether to use
USE_AE_AVERAGE_AS_REFERENCE = True
# The following line is ony used if USE_AE_AVERAGE_AS_REFERENCE is False
REFERENCE_CODE_LABEL = None
SET_MAX_SCALE_DICT = {}
ONLY_CODES = None

# ## IN ORDER TO PLOT ONLY THE AE COMPARISON
# # Whether to use
# USE_AE_AVERAGE_AS_REFERENCE = False
# # The following line is ony used if USE_AE_AVERAGE_AS_REFERENCE is False
# REFERENCE_CODE_LABEL = "FLEUR"
# SET_MAX_SCALE_DICT = {'nu': 0.350000000001, 'epsilon': 0.2}
# ONLY_CODES = ['WIEN2k']

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


def abs_V0_rel_diff(*args, **kwargs):
    return abs(qc.V0_rel_diff(*args, **kwargs))
def abs_B0_rel_diff(*args, **kwargs):
    return abs(qc.B0_rel_diff(*args, **kwargs))
def abs_B1_rel_diff(*args, **kwargs):
    return abs(qc.B1_rel_diff(*args, **kwargs))

quantity_for_comparison_map = {
    "delta_per_formula_unit": qc.delta,
    "delta_per_formula_unit_over_b0": qc.delta_over_b0,
    "B0_rel_diff": qc.B0_rel_diff,
    "V0_rel_diff": qc.V0_rel_diff,
    "B1_rel_diff": qc.B1_rel_diff,
    "abs_V0_rel_diff": abs_V0_rel_diff,
    "abs_B0_rel_diff": abs_B0_rel_diff,
    "abs_B1_rel_diff": abs_B1_rel_diff,            
    f"nu": qc.nu,
    "epsilon": qc.epsilon
}


def plot_periodic_table(SET_NAME, QUANTITY):
    SET_MAX_SCALE = SET_MAX_SCALE_DICT.get(QUANTITY, None)
    prefactor = PREFACTOR_DICT.get(QUANTITY, 1.)

    DATA_FOLDER = "../../../code-data"
    with open(os.path.join(DATA_FOLDER, "labels.json")) as fhandle:
        labels_data = json.load(fhandle)
    
    if USE_AE_AVERAGE_AS_REFERENCE:
        reference_data_files = labels_data['references']['all-electron average']
        reference_short_label = "ae"
    else:
        reference_data_files = labels_data['methods-main'][REFERENCE_CODE_LABEL]
        reference_short_label = labels_data['methods-main'][REFERENCE_CODE_LABEL]['short_label']
    try:
        with open(os.path.join(DATA_FOLDER, reference_data_files[SET_NAME])) as fhandle:
            compare_plugin_data = json.load(fhandle)
            if not compare_plugin_data['script_version'] in EXPECTED_SCRIPT_VERSION:
                raise ValueError(
                    f"This script only works with data generated at version {EXPECTED_SCRIPT_VERSION}. "
                    f"Please re-run ./get_results.py to update the data format for the all-electron dataset!"
                    )
                #sys.exit(1)
    except OSError:
        print(f"No data found for the all-electron dataset (set '{SET_NAME}'), it is the reference and must be present")
        sys.exit(1)


    code_results = {}
    short_labels = {}
    for code_label in labels_data['methods-main']:
        if ONLY_CODES is not None and code_label not in ONLY_CODES:
            continue
        short_labels[code_label] = labels_data['methods-main'][code_label]['short_label']
        with open(os.path.join(DATA_FOLDER, labels_data['methods-main'][code_label][SET_NAME])) as fhandle:
            code_results[code_label] = json.load(fhandle)
            if not code_results[code_label]['script_version'] in EXPECTED_SCRIPT_VERSION:
                raise ValueError(
                    f"This script only works with data generated at version {EXPECTED_SCRIPT_VERSION}. "
                    f"Please re-run ./get_results.py to update the data format for {code_label}! Skipping it"
                    )
                #code_results.pop(label)

    for plugin, plugin_data in code_results.items():

        print(f"Using data for method '{plugin}' (set '{SET_NAME}') compared with {reference_short_label}.")


        all_systems = set(plugin_data['eos_data'].keys())
        all_systems = set(plugin_data['BM_fit_data'].keys())
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

        progress_bar = tqdm.tqdm(sorted(all_systems))
        for element_and_configuration in progress_bar:
            progress_bar.set_description(f"{element_and_configuration:12s}")
            progress_bar.refresh()

            element, configuration = element_and_configuration.split('-')
            # Get the data for the reference plugin
            ref_BM_fit_data = plugin_data['BM_fit_data'][f'{element}-{configuration}']
        
            if ref_BM_fit_data is None:
                continue
        
            scaling_factor_ref = qc.get_volume_scaling_to_formula_unit(
                    plugin_data['num_atoms_in_sim_cell'][f'{element}-{configuration}'],
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

            quant = quantity_for_comparison_map[QUANTITY](V0,B0,B01,CV0,CB0,CB01,prefactor,DEFAULT_wb0,DEFAULT_wb1)

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
        alpha = 0.7
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
        elif cmap == "viridis":
            cmap = viridis
            bokeh_palette = "Viridis256"
        elif cmap == "inferno":
            cmap = inferno
            bokeh_palette = "Inferno256"
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

        # I put a large (for min) and small (default for empty sets)
        min_data = min([min(collect[i]["values"], default=100) for i in list_confs])
        max_data = max([max(collect[i]["values"], default=0) for i in list_confs])

        # If it's reallly all empty
        if min_data > max_data:
            min_data = max_data

        color_list={}

        if SET_MAX_SCALE:
            high = SET_MAX_SCALE
        else:
            high = max_data

        # EXPORT JSON: a dictionary with key = element+config, value = measure
        data_to_export = {}
        for conf in list_confs:
            
            data_to_export.update(dict(zip(
                (f'{element}-{conf}' for element in collect[conf]["elements"]),
                collect[conf]["values"])))
        with open(f"{QUANTITY}-{SET_NAME}-{short_labels[plugin].replace(' ', '_')}-vs-{reference_short_label.replace(' ', '_')}.json", 'w') as fhandle:
            json.dump(data_to_export, fhandle)

        non_excellent = []
        tot_count = 0
        for conf in list_confs:
            data_elements = collect[conf]["elements"]
            tot_count += len(data_elements)
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
                    palette=bokeh_palette, low=min_data, high=high
                )
                norm = LogNorm(vmin=min_data, vmax=high)
            else:
                low = 0.
                color_mapper = LinearColorMapper(
                    palette=bokeh_palette, low=low, high=high
                )
                norm = Normalize(vmin=low, vmax=high)
                
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
                if data[i] > high:
                    print(f"** WARNING! {data_element}-{conf} has value {data[i]} > max of colorbar ({high})")
                if data[i] > EXCELLENT_AGREEMENT_THRESHOLD[QUANTITY]:
                    non_excellent.append(f"{data_element}({conf})")
        if PRINT_NON_EXCELLENT:
            print(f">>> Non excellent agreement ({QUANTITY} >= {EXCELLENT_AGREEMENT_THRESHOLD[QUANTITY]}) for {len(non_excellent)}/{tot_count} systems: {','.join(non_excellent)}")


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
            p.toolbar.logo = None
            p.toolbar.tools = []
            p.toolbar_location = None
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
            p.toolbar.logo = None
            p.toolbar.tools = []
            p.toolbar_location = None
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

        reference_label = 'all-electron average' if USE_AE_AVERAGE_AS_REFERENCE else REFERENCE_CODE_LABEL
        p.title = f"{QUANTITY} for {plugin} vs. {reference_label}"
        p.title.text_font_size = '16pt'

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
                export_png(p, filename=f"periodic-table-{SET_NAME}-{short_labels[plugin].replace(' ', '_')}-vs-{reference_short_label.replace(' ', '_')}-{QUANTITY}.png")
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



if __name__ == "__main__":
    
    for SET_NAME in ['unaries', 'oxides']:
        for QUANTITY in [
            'epsilon', 'nu', 'delta_per_formula_unit', 
            'delta_per_formula_unit_over_b0']:
            plot_periodic_table(SET_NAME, QUANTITY)