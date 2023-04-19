#!/usr/bin/env python
import json
import os
import sys

import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
import tqdm

import acwf_paper_plots.quantities_for_comparison as qc

# As found in the paper, nu and eps can be roughly related via just a multiplication: nu=NU_EPS_FACTOR*eps
# Use this to set a consistent maximum colorbar value
NU_EPS_FACTOR=1.65

SHOW_IN_BROWSER=False
DEFAULT_wb0 = 1.0/20.0
DEFAULT_wb1 = 1.0/400.0
# Default prefactor if not indicated: 1.
PREFACTOR_DICT = {'nu': 100.}
EXPECTED_SCRIPT_VERSION = ["0.0.3","0.0.4"]
# NOTE! in the code, I call the function e.g. 'delta_per_formula_unit', but in reality I then already divide by
# the number of atoms in the formula unit, so the numbers I get are per atom.
# Therefore, the UNICODE name has 'per atom' since it is shown in the final plot
UNICODE_QUANTITY = {'nu': 'ν', 'epsilon': 'ε', 'delta_per_formula_unit': 'Δ per atom', 'delta_per_formula_unit_over_b0': 'Δ/B₀ per atom'}
EXCELLENT_AGREEMENT_THRESHOLD = {
    'nu': 0.10, 'epsilon': 0.06,
    'delta_per_formula_unit': 0., # I put zero, it's not used in this script anyway
    'delta_per_formula_unit_over_b0': 0. # I put zero, it's not used in this script anyway
    }
GOOD_AGREEMENT_THRESHOLD = {
    'nu': 0.33, 'epsilon': 0.20,
    'delta_per_formula_unit': 0., # I put zero, it's not used in this script anyway
    'delta_per_formula_unit_over_b0': 0. # I put zero, it's not used in this script anyway
    }
OUTLIER_THRESHOLD = {
    'nu': 1.0 * NU_EPS_FACTOR, 'epsilon': 1.0,
    'delta_per_formula_unit': 0., # I put zero, it's not used in this script anyway
    'delta_per_formula_unit_over_b0': 0. # I put zero, it's not used in this script anyway
    }
PRINT_NON_EXCELLENT = False

## --------------------------------------------------
## "Constants" that might need to be changed, depeding on what Figure is generated

# Whether to use
USE_AE_AVERAGE_AS_REFERENCE = True
# The following line is ony used if USE_AE_AVERAGE_AS_REFERENCE is False
REFERENCE_CODE_LABEL = "FLEUR@LAPW+LO"
SKIP_PLOT_FOR_QUANTITIES = ['delta_per_formula_unit', 'delta_per_formula_unit_over_b0']
LABELS_KEY = 'methods-main'
ONLY_CODES = None #["CASTEP@PW|C19MK2", "Quantum ESPRESSO@PW|SSSP-prec-v1.3"] #["ABINIT@PW|PseudoDojo-v0.5", "BigDFT@DW|HGH-K(Valence)"]

CBAR_MAX_DICT = {}

CMAP_TYPE = "quality"
SET_MAX_SCALE_DICT = {"nu": 1.0*NU_EPS_FACTOR, "epsilon":1.0}
OUTLIER_COLOR = "#bf0000" # darker red
#CBAR_MAX_DICT = {"nu": 0.4*NU_EPS_FACTOR, "epsilon":0.4}
HIGHLIGHT = {}

EXPORT_JSON=False

PRINT_LATEX_CODE=False

SET_NAMES = ['unaries', 'oxides']
QUANTITIES = ['epsilon', 'nu', 'delta_per_formula_unit', 'delta_per_formula_unit_over_b0']

## ------------------------------------------------------------------------------------------------
## Override the default variables based on the input argument

if len(sys.argv) == 2:
    if sys.argv[1] == "MAIN":
        # FIGURE 2 IN MAIN TEXT
        USE_AE_AVERAGE_AS_REFERENCE = False
        REFERENCE_CODE_LABEL = "FLEUR@LAPW+LO"
        LABELS_KEY = 'methods-main'
        ONLY_CODES = ["WIEN2k@(L)APW+lo+LO"]
        CBAR_MAX_DICT = {"nu": 0.4*NU_EPS_FACTOR, "epsilon":0.4}

    if sys.argv[1] == "SI-all-tables":
        # Section S14
        USE_AE_AVERAGE_AS_REFERENCE = True
        LABELS_KEY = 'methods-main'
        ONLY_CODES = None
        EXPORT_JSON=True
        PRINT_LATEX_CODE=True

    if sys.argv[1] == "SI-29-vs-960-highlight":
        # Figure S39
        USE_AE_AVERAGE_AS_REFERENCE = True
        LABELS_KEY = 'methods-main'
        ONLY_CODES = ["CASTEP@PW|C19MK2", "Quantum ESPRESSO@PW|SSSP-prec-v1.3"]
        QUANTITIES=["epsilon"]
        HIGHLIGHT = {
            "unaries": {
                "epsilon" : {
                    "CASTEP@PW|C19MK2": {
                        "X/SC": ["Po"],
                        "X/FCC": ["Ne", "Al", "Ar", "Ca", "Cu", "Kr", "Sr", "Rh", "Pd", "Ag", "Xe", "Ir", "Pt", "Au", "Pb", "Rn"],
                        "X/BCC": ["K", "V", "Rb", "Nb", "Mo", "Cs", "Ba", "Ta", "W"],
                        "X/Diamond": ["Si", "Ge", "Sn"],
                        },
                    "Quantum ESPRESSO@PW|SSSP-prec-v1.3": {
                        "X/SC": ["Po"],
                        "X/FCC": ["Ne", "Al", "Ar", "Ca", "Cu", "Kr", "Sr", "Rh", "Pd", "Ag", "Xe", "Ir", "Pt", "Au", "Pb", "Rn"],
                        "X/BCC": ["K", "V", "Rb", "Nb", "Mo", "Cs", "Ba", "Ta", "W"],
                        "X/Diamond": ["Si", "Ge", "Sn"],
                        },
                }
            }
        }

    if sys.argv[1] == "SI-VASP-1":
        # S27
        USE_AE_AVERAGE_AS_REFERENCE = True
        LABELS_KEY = 'methods-supplementary'
        ONLY_CODES = ["VASP@PW|PBErec-PAW54*|defCutoff", "VASP@PW|PBErec-PAW54*|800Cutoff"]
        QUANTITIES=["epsilon"]
        SET_NAMES = ['unaries']

    if sys.argv[1] == "SI-VASP-2":
        # S27
        USE_AE_AVERAGE_AS_REFERENCE = True
        LABELS_KEY = 'methods-main'
        ONLY_CODES = ["VASP@PW|GW-PAW54*"]
        QUANTITIES=["epsilon"]
        SET_NAMES = ['unaries']


    if sys.argv[1] == "SI-PSEUDODOJO-SECTION-1":
        # Section S16
        USE_AE_AVERAGE_AS_REFERENCE = False
        REFERENCE_CODE_LABEL = "ABINIT@PW|PseudoDojo-v0.4"
        LABELS_KEY = 'methods-supplementary'
        ONLY_CODES = ["CASTEP@PW|PseudoDojo-v0.4-trim"]
        QUANTITIES=["epsilon"]

    if sys.argv[1] == "SI-PSEUDODOJO-SECTION-2":
        # Section S16
        USE_AE_AVERAGE_AS_REFERENCE = False
        REFERENCE_CODE_LABEL = "ABINIT@PW|PseudoDojo-v0.4"
        LABELS_KEY = 'methods-supplementary'
        ONLY_CODES = ["Quantum ESPRESSO@PW|PseudoDojo-v0.4-trim"]
        QUANTITIES=["epsilon"]

    if sys.argv[1] == "SI-PSEUDODOJO-SECTION-3":
        # Section S16
        USE_AE_AVERAGE_AS_REFERENCE = False
        REFERENCE_CODE_LABEL = "CASTEP@PW|PseudoDojo-v0.4-trim"
        LABELS_KEY = 'methods-supplementary'
        ONLY_CODES = ["Quantum ESPRESSO@PW|PseudoDojo-v0.4-trim"]
        QUANTITIES=["epsilon"]

    if sys.argv[1] == "SI-PSEUDODOJO-SECTION-4":
        # Section S16
        USE_AE_AVERAGE_AS_REFERENCE = False
        REFERENCE_CODE_LABEL = "SIRIUS/CP2K@PW|SSSP-prec-v1.2"
        LABELS_KEY = 'methods-supplementary'
        ONLY_CODES = ["Quantum ESPRESSO@PW|SSSP-prec-v1.2"]
        QUANTITIES=["epsilon"]

## ------------------------------------------------------------------------------------------------

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
from bokeh.colors import RGB
from matplotlib.colors import Normalize, LogNorm, to_hex, LinearSegmentedColormap
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

def make_quality_matching_cmap(quantity):
    """
    Custom colormap matching the excellent/good/bad thresholds
    """
    exc_thresh = EXCELLENT_AGREEMENT_THRESHOLD[quantity]
    good_thresh = GOOD_AGREEMENT_THRESHOLD[quantity]
    outl_thresh = OUTLIER_THRESHOLD[quantity]

    colorbar_max = 1.04*outl_thresh
    cvals  = [0.0, exc_thresh, good_thresh, outl_thresh, outl_thresh+0.001, colorbar_max]
    colors = ["#0000be", "#3a50de","#ffff55", "#f53216", "#bf0000", "#bf0000"]

    norm = Normalize(min(cvals),max(cvals))
    tuples = list(zip(map(norm,cvals), colors))
    cmap = LinearSegmentedColormap.from_list("", tuples, N=256)

    num_colors= 256
    high = max(cvals)

    if quantity in CBAR_MAX_DICT:
        # cap the colorbar at max_value
        num_colors = int(round(CBAR_MAX_DICT[quantity]/colorbar_max*256))
        high = CBAR_MAX_DICT[quantity]

    custom_rgb = (255 * cmap(range(num_colors))).astype('int')
    bokeh_palette = [RGB(*tuple(rgb)).to_hex() for rgb in custom_rgb]

    color_mapper = LinearColorMapper(
                palette=bokeh_palette, low=min(cvals), high=high
            )

    return norm, cmap, color_mapper

def make_simple_cmap(data, high, min_data, cmap_name="plasma", log_scale=False):

    cmap = None

    # Assign color palette based on input argument
    if cmap_name == "plasma":
        cmap = plasma
        bokeh_palette = "Plasma256"
    elif cmap_name == "magma":
        cmap = magma
        bokeh_palette = "Magma256"
    elif cmap_name == "viridis":
        cmap = viridis
        bokeh_palette = "Viridis256"
    elif cmap_name == "inferno":
        cmap = inferno
        bokeh_palette = "Inferno256"
    else:
        raise ValueError("Unknown color map")

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

    return norm, cmap, color_mapper


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
    "nu": qc.nu,
    "epsilon": qc.epsilon
}

def load_data(SET_NAME):

    DATA_FOLDER = "../../../code-data"
    with open(os.path.join(DATA_FOLDER, "labels.json")) as fhandle:
        labels_data = json.load(fhandle)
    
    if USE_AE_AVERAGE_AS_REFERENCE:
        reference_data_files = labels_data['references']['all-electron average']
        reference_short_label = "ae"
    else:
        reference_data_files = labels_data[LABELS_KEY][REFERENCE_CODE_LABEL]
        reference_short_label = labels_data[LABELS_KEY][REFERENCE_CODE_LABEL]['short_label']
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
    for code_label in labels_data[LABELS_KEY]:
        if ONLY_CODES is not None and code_label not in ONLY_CODES:
            continue
        short_labels[code_label] = labels_data[LABELS_KEY][code_label]['short_label']
        with open(os.path.join(DATA_FOLDER, labels_data[LABELS_KEY][code_label][SET_NAME])) as fhandle:
            code_results[code_label] = json.load(fhandle)
            if not code_results[code_label]['script_version'] in EXPECTED_SCRIPT_VERSION:
                raise ValueError(
                    f"This script only works with data generated at version {EXPECTED_SCRIPT_VERSION}. "
                    f"Please re-run ./get_results.py to update the data format for {code_label}! Skipping it"
                    )
                #code_results.pop(label)

    loaded_data = {
        "code_results": code_results,
        "short_labels": short_labels,
        "reference_short_label": reference_short_label,
        "compare_plugin_data": compare_plugin_data
    }

    return loaded_data


def calculate_quantities(plugin_data, compare_plugin_data, QUANTITY):
    prefactor = PREFACTOR_DICT.get(QUANTITY, 1.)

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

        # Here I normalize quantities, so that they are now per atom and not per formula unit!
        # This does not change anything for epsilon and nu, but changes for delta
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

        # Here I normalize quantities, so that they are now per atom and not per formula unit!
        # This does not change anything for epsilon and nu, but changes for delta
        CV0=compare_BM_fit_data['min_volume']/scaling_factor_comp
        CB0=compare_BM_fit_data['bulk_modulus_ev_ang3']
        CB01=compare_BM_fit_data['bulk_deriv']

        quant = quantity_for_comparison_map[QUANTITY](V0,B0,B01,CV0,CB0,CB01,prefactor,DEFAULT_wb0,DEFAULT_wb1)

        collect[configuration]["values"].append(quant)
        collect[configuration]["elements"].append(element)

    return collect

def export_json_file(SET_NAME, QUANTITY, collect, list_confs, short_labels, plugin, reference_short_label):
    # EXPORT JSON: a dictionary with key = element+config, value = measure
    data_to_export = {}
    for conf in list_confs:
        
        data_to_export.update(dict(zip(
            (f'{element}-{conf}' for element in collect[conf]["elements"]),
            collect[conf]["values"])))
    with open(f"{QUANTITY}-{SET_NAME}-{short_labels[plugin].replace(' ', '_')}-vs-{reference_short_label.replace(' ', '_')}.json", 'w') as fhandle:
        json.dump(data_to_export, fhandle)

def create_periodic_table(SET_NAME, QUANTITY, collect, list_confs, short_labels, plugin, reference_short_label, unaries, SET_MAX_SCALE):

    width = 1050
    width_cbar = 80 # needs to be manually adjusted to make the quads square...
    alpha = 0.7
    extended = True
    cbar_height = None
    cbar_standoff = 12
    cbar_fontsize = 14
    blank_color = "#c4c4c4"
    under_value = None
    under_color = "#140F0E"

    #over_value = None
    over_value = SET_MAX_SCALE

    over_color = OUTLIER_COLOR # "#140F0E"
    special_elements = None
    special_color = "#6F3023"

    # following is only used for "simple" colormap
    cmap_name = "plasma"
    log_scale = False

    options.mode.chained_assignment = None

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
    highlight_list = {}

    if SET_MAX_SCALE:
       high = SET_MAX_SCALE
    else:
       high = max_data

    non_excellent = []
    tot_count = 0
    for conf in list_confs:
        data_elements = collect[conf]["elements"]
        tot_count += len(data_elements)
        data = collect[conf]["values"]

        if len(data) != len(data_elements):
            raise ValueError("Unequal number of atomic elements and data points")

        if CMAP_TYPE == "simple":
            norm, cmap, color_mapper = make_simple_cmap(data, high, min_data, cmap_name=cmap_name, log_scale=log_scale)
        elif CMAP_TYPE == "quality":
            norm, cmap, color_mapper = make_quality_matching_cmap(QUANTITY)
        else:
            raise ValueError("Unknown colormap type!")

        color_scale = ScalarMappable(norm=norm, cmap=cmap).to_rgba(data, alpha=None)

        # Set blank color
        color_list[conf] = [blank_color] * len(elements)

        # list of line widths for highlight
        highlight_list[conf] = [0] * len(elements)

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
            
            try:                
                highlight = data_element in HIGHLIGHT[SET_NAME][QUANTITY][plugin][conf]
                highlight_list[conf][element_index] = 1.0 if highlight else 0.0
            except KeyError:
                pass
            
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
                la_dia=highlight_list["X/Diamond"],
                la_sc=highlight_list["X/SC"],
                la_bcc=highlight_list["X/BCC"],
                la_fcc=highlight_list["X/FCC"],
            )
        )

        # Plot the periodic table
        p = figure(x_range=[0,19], y_range=[11,0], tools="save", match_aspect=True)
        p.toolbar.logo = None
        p.toolbar.tools = []
        p.toolbar_location = None
        p.plot_width = width
        p.outline_line_color = None
        p.background_fill_color = None
        p.border_fill_color = None
        p.toolbar_location = "above"
        quad_args = {
            "source": source,
            "alpha": alpha,
        }
        p.quad(left="left", right="group", top="period", bottom="bottom", color="type_color_dia", **quad_args)
        p.quad(left="left", right="group", top="top", bottom="period", color="type_color_sc", **quad_args)
        p.quad(left="group", right="right", top="period", bottom="bottom", color="type_color_bcc", **quad_args)
        p.quad(left="group", right="right", top="top", bottom="period", color="type_color_fcc", **quad_args)
        # border around each box:
        p.quad(left="left", right="right", top="top", bottom="bottom", line_width=1, line_color="#b5b5b5", fill_color=None, source=source)
        # separate boxes for highlight border
        # very hacky but otherwise the colors are overlapping on the highlight border and it doesn't look nice
        quad_args_hl = {
            "source": source,
            "alpha": alpha,
            "line_color": "lime",
            "line_width": 6,
            "color": None,
        }
        p.quad(left="left", right="group", top="period", bottom="bottom", line_alpha="la_dia", **quad_args_hl)
        p.quad(left="left", right="group", top="top", bottom="period", line_alpha="la_sc", **quad_args_hl)
        p.quad(left="group", right="right", top="period", bottom="bottom", line_alpha="la_bcc", **quad_args_hl)
        p.quad(left="group", right="right", top="top", bottom="period", line_alpha="la_fcc", **quad_args_hl)

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
        p.plot_width = width - width_cbar
        p.outline_line_color = None
        p.background_fill_color = None
        p.border_fill_color = None
        p.toolbar_location = "above"
        quad_args = {
            "source": source,
            "alpha": alpha,
        #    "line_width": 1,
        #    "line_color": "#b5b5b5",
        }
        p.quad(left="left", right="group", top="top", bottom="midup", color="type_color_X2O3", **quad_args)
        p.quad(left="left", right="group", top="midup", bottom="middown", color="type_color_X2O", **quad_args)
        p.quad(left="left", right="group", top="middown", bottom="bottom", color="type_color_XO3", **quad_args)
        p.quad(left="group", right="right", top="top", bottom="midup", color="type_color_X2O5", **quad_args)
        p.quad(left="group", right="right", top="midup", bottom="middown", color="type_color_XO2", **quad_args)
        p.quad(left="group", right="right", top="middown", bottom="bottom", color="type_color_XO", **quad_args)
        # border around each box:
        p.quad(left="left", right="right", top="top", bottom="bottom", line_width=1, line_color="#b5b5b5", fill_color=None, **quad_args)
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
    p.title = f"{UNICODE_QUANTITY[QUANTITY]} for {plugin} vs. {reference_label}"
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

    # Skip the colorbar for oxides
    if unaries:
        p.add_layout(color_bar, "right")
    #p.add_layout(color_bar, "right")
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




def plot_periodic_tables(SET_NAME, QUANTITY, measures_max_and_avg, master_data_dict):

    ld = master_data_dict[SET_NAME]["loaded_data"]

    for plugin, plugin_data in ld["code_results"].items():

        print(f"Using data for method '{plugin}' (set '{SET_NAME}') compared with {ld['reference_short_label']}.")


        collect = master_data_dict[SET_NAME]["calculated_quantities"][QUANTITY][plugin]

        # Way to decide whether is a unaries or an oxides set is a bit fragile.
        if collect["X/Diamond"]["values"]:
            unaries = True
            list_confs = ["X/Diamond","X/FCC","X/BCC","X/SC"]
        else:
            unaries = False
            list_confs = ["X2O3","X2O5","X2O","XO2","XO3","XO"]

        if EXPORT_JSON:
            export_json_file(SET_NAME, QUANTITY, collect, list_confs, ld["short_labels"], plugin, ld["reference_short_label"])

        if QUANTITY in SKIP_PLOT_FOR_QUANTITIES:
            continue

        ### Determine maximum scale for the colorbar
        if QUANTITY not in SET_MAX_SCALE_DICT:
            # 0. Default (None) will take the data maximum of each periodic table separately
            SET_MAX_SCALE = None
        elif isinstance(SET_MAX_SCALE_DICT[QUANTITY], (int, float)):
            # 1. User defined a fixed maximum
            SET_MAX_SCALE = SET_MAX_SCALE_DICT[QUANTITY]
        elif SET_MAX_SCALE_DICT[QUANTITY] == 'max':
            # 2. consistent maximum for NU/EPS
            nu_max = measures_max_and_avg[plugin]["nu"]["max"]
            eps_max = measures_max_and_avg[plugin]["epsilon"]["max"]
            effective_nu_max = max(nu_max, NU_EPS_FACTOR * eps_max)

            if QUANTITY == "nu":
                SET_MAX_SCALE = effective_nu_max
            elif QUANTITY == "epsilon":
                SET_MAX_SCALE = effective_nu_max / NU_EPS_FACTOR

        elif isinstance(SET_MAX_SCALE_DICT[QUANTITY], tuple) and SET_MAX_SCALE_DICT[QUANTITY][0] == "avg":
            # 3. take the N*avg epsilon and use that for the maximum of NU and EPS
            n = SET_MAX_SCALE_DICT[QUANTITY][1]
            eps_max_scale = n*measures_max_and_avg[plugin]["epsilon"]["avg"]
            if QUANTITY == "nu":
                SET_MAX_SCALE = eps_max_scale * NU_EPS_FACTOR
            elif QUANTITY == "epsilon":
                SET_MAX_SCALE = eps_max_scale
        else:
            raise ValueError("Unknown max scale type!")

        create_periodic_table(SET_NAME, QUANTITY, collect, list_confs, ld["short_labels"], plugin, ld["reference_short_label"], unaries, SET_MAX_SCALE)


def find_code_measures_max_and_avg(master_data_dict):
    """
    For every code, we plot 4 periodic tables: unaries and oxides for nu and epsilon.

    Calculate the maximum and average of each code for each quantity over all the unaries
    and oxides to allow later for multiple colorscale options
    """

    tmp = {}

    for SET_NAME in SET_NAMES:

        ld = master_data_dict[SET_NAME]["loaded_data"]

        for QUANTITY in QUANTITIES:

            for plugin, plugin_data in ld["code_results"].items():

                if plugin not in tmp:
                    tmp[plugin] = {}
                if QUANTITY not in tmp[plugin]:
                    tmp[plugin][QUANTITY] = {"max": 0.0, "total": 0.0, "count": 0}
                
                collect = master_data_dict[SET_NAME]["calculated_quantities"][QUANTITY][plugin]

                # find the maximum and total/count in all of the collected data
                for configuration, conf_data in collect.items():
                    if len(conf_data['values']) == 0:
                        continue
                    current_conf_max_val = max(conf_data['values'])
                    tmp[plugin][QUANTITY]["max"] = max(tmp[plugin][QUANTITY]["max"], current_conf_max_val)
                    tmp[plugin][QUANTITY]["total"] += sum(conf_data['values'])
                    tmp[plugin][QUANTITY]["count"] += len(conf_data['values'])

    # build the final dict by calculating the avg
    measures_max_and_avg = {}
    for plugin in tmp:
        measures_max_and_avg[plugin] = {}
        for QUANTITY in tmp[plugin]:
            measures_max_and_avg[plugin][QUANTITY] = {
                "max": tmp[plugin][QUANTITY]["max"],
                "avg": tmp[plugin][QUANTITY]["total"]/tmp[plugin][QUANTITY]["count"],
            }
    return measures_max_and_avg


def analyze_stats(master_data_dict):

    stats = {}

    for SET_NAME in SET_NAMES:

        for QUANTITY in QUANTITIES:

            calc_q = master_data_dict[SET_NAME]["calculated_quantities"][QUANTITY]

            for plugin in calc_q:

                if plugin not in stats:
                    stats[plugin] = {}
                if QUANTITY not in stats[plugin]:
                    stats[plugin][QUANTITY] = {
                        "total":        0,
                        "excellent":    0, # 0 to excellent
                        "good":         0, # excellent to good
                        "different":    0, # good to outlier
                        "outlier":      0, # over outlier threshold
                    }

                collect = master_data_dict[SET_NAME]["calculated_quantities"][QUANTITY][plugin]

                for configuration, conf_data in collect.items():
                    vals_arr = np.array(conf_data["values"])
                    d = stats[plugin][QUANTITY]
                    d["total"] += len(vals_arr)
                    d["excellent"] += sum(vals_arr <= EXCELLENT_AGREEMENT_THRESHOLD[QUANTITY])
                    d["good"] += sum(np.logical_and(
                        EXCELLENT_AGREEMENT_THRESHOLD[QUANTITY] < vals_arr,
                        vals_arr <= GOOD_AGREEMENT_THRESHOLD[QUANTITY]
                    ))
                    d["different"] += sum(np.logical_and(
                        GOOD_AGREEMENT_THRESHOLD[QUANTITY] < vals_arr,
                        vals_arr <= OUTLIER_THRESHOLD[QUANTITY]
                    ))
                    d["outlier"] += sum(vals_arr > OUTLIER_THRESHOLD[QUANTITY])

    if PRINT_LATEX_CODE:
        # print the latex lines for the captions of S14.
    
        # map from plugin name (before @) to the latex convention
        latex_names = {
            "ABINIT": ("\\abinitlong", "ABINIT"),
            "BigDFT": ("\\bigdftlong", "BigDFT"),
            "CASTEP": ("\\casteplong", "CASTEP"),
            "CP2K/Quickstep": ("\\cptwoklong", "CP2K-quickstep"),
            "FLEUR": ("\\fleurlong", "FLEUR"),
            "GPAW": ("\\gpawlong", "GPAW"),
            "Quantum ESPRESSO": ("\\qelong", "QE"),
            "SIESTA": ("\\siestalong", "SIESTA"),
            "SIRIUS/CP2K": ("\\siriuslong", "SIRIUS-CP2K"),
            "VASP": ("\\vasplong", "VASP"),
            "WIEN2k": ("\\wientwoklong", "WIEN2k"),
        }

        print()
        print("Statistics for the captions")
        print("Copy-paste these lines into the latex code:")
        print("----")
        for plugin in stats:

            assert stats[plugin]["epsilon"]["total"] == stats[plugin]["nu"]["total"]
             
            label_long, label = latex_names[plugin.split("@")[0]]
            
            s = "\\singleapproachtemplate{" + f"{label_long}, {label}"
            s += f""", {stats[plugin]["epsilon"]["total"]}"""
            for quantity in ["epsilon", "nu"]:
                for key in ["excellent", "good", "different", "outlier"]:
                    s += f", {stats[plugin][quantity][key]}"
            s += "}"
            print(s)

        print("----")
        print()


if __name__ == "__main__":

    master_data_dict = {}

    for SET_NAME in SET_NAMES:
        ld = load_data(SET_NAME)

        master_data_dict[SET_NAME] = {
            "loaded_data": ld,
            "calculated_quantities": {}
        }
        for QUANTITY in QUANTITIES:
            master_data_dict[SET_NAME]["calculated_quantities"][QUANTITY] = {}
            for plugin, plugin_data in ld["code_results"].items():
                collect = calculate_quantities(plugin_data, ld["compare_plugin_data"], QUANTITY)
                master_data_dict[SET_NAME]["calculated_quantities"][QUANTITY][plugin] = collect

    measures_max_and_avg = find_code_measures_max_and_avg(master_data_dict)


    print("Plotting the periodic tables.")
    for SET_NAME in SET_NAMES:
        for QUANTITY in QUANTITIES:
            plot_periodic_tables(SET_NAME, QUANTITY, measures_max_and_avg, master_data_dict)

    analyze_stats(master_data_dict)
