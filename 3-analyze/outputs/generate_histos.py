#!/usr/bin/env python
import json
import os
import sys

import numpy as np
import pylab as pl
import tqdm

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

BINS = 100


def birch_murnaghan(V,E0,V0,B0,B01):
    """
    Return the energy for given volume (V - it can be a vector) according to
    the Birch Murnaghan function with parameters E0,V0,B0,B01.
    """
    r = (V0/V)**(2./3.)
    return (E0 +
            9./16. * B0 * V0 * (
            (r-1.)**3 * B01 +
            (r-1.)**2 * (6. - 4.* r)))


def intE12sq(v0w,b0w,b1w,v0f,b0f,b1f,V1,V2):
    """
    Integral of (E1(V) - E2(V))**2 in dV evaluated between volume V1 and volume V2
    """
    F1 = antiderE12sq(v0w,b0w,b1w,v0f,b0f,b1f,V1)
    F2 = antiderE12sq(v0w,b0w,b1w,v0f,b0f,b1f,V2)
    integral = F2 - F1

    return integral

def antiderE12sq(v0w,b0w,b1w,v0f,b0f,b1f,V):
    """
    Antiderivative of (E1(V) - E2(V))**2 where E1(V) and E2(V) are birch murnaghan
    functions with different parameters
    """
    antider = (81*(\
            6*b0w*b0f*(-16 + 3*b1w)*(-16 + 3*b1f)*V*v0w*(v0w/V)**(2/3)*v0f*(v0f/V)**(2/3) - \
            2*b0w*b0f*(-14 + 3*b1w)*(-16 + 3*b1f)*V*v0w*(v0w/V)**(4/3)*v0f*(v0f/V)**(2/3) - \
            2*b0w*b0f*(-16 + 3*b1w)*(-14 + 3*b1f)*V*v0w*(v0w/V)**(2/3)*v0f*(v0f/V)**(4/3) + \
            (6*b0w*b0f*(-14 + 3*b1w)*(-14 + 3*b1f)*V*v0w*(v0w/V)**(4/3)*v0f*(v0f/V)**(4/3))/5. + \
            V*(b0w*(-6 + b1w)*v0w - b0f*(-6 + b1f)*v0f)**2 - \
            (b0w*(-4 + b1w)*v0w**3 - b0f*(-4 + b1f)*v0f**3)**2/(3.*V**3) + \
            (3*b0f*(v0f/V)**(7/3)*(-2*b0w*(-14 + 3*b1f)*v0w*(-7*(-6 + b1w)*V**2 + \
            (-4 + b1w)*v0w**2) \
            - 7*b0f*(424 + 5*b1f*(-32 + 3*b1f))*V**2*\
            v0f + 2*b0f*(-4 + b1f)*(-14 + 3*b1f)*\
            v0f**3))/7. - \
            (3*b0f*(v0f/V)**(5/3)*\
            (-2*b0w*(-16 + 3*b1f)*v0w*\
            (5*(-6 + b1w)*V**2 + (-4 + b1w)*v0w**2) \
            + 10*b0f*(-6 + b1f)*(-16 + 3*b1f)*V**2*\
            v0f + b0f*(324 + 5*b1f*(-28 + 3*b1f))*\
            v0f**3))/5. +\
            (4*b0w**2*(124 + 5*(-10 + b1w)*b1w)*v0w**4 - \
            2*b0w*b0f*(-4 + b1w)*(-6 + b1f)*v0w**3*\
            v0f - 2*b0w*b0f*(-6 + b1w)*(-4 + b1f)*v0w*\
            v0f**3 + 4*b0f**2*(124 + 5*(-10 + b1f)*b1f)*v0f**4)/V + \
            (3*b0w*(v0w/V)**(7/3)*\
            (-7*b0w*(424 + 5*b1w*(-32 + 3*b1w))*V**2*\
            v0w + 2*b0w*(-4 + b1w)*(-14 + 3*b1w)*v0w**3 + \
            2*b0f*(-14 + 3*b1w)*v0f*\
            (7*(-6 + b1f)*V**2 - (-4 + b1f)*v0f**2)))/7. - \
            (3*b0w*(v0w/V)**(5/3)*\
            (10*b0w*(-6 + b1w)*(-16 + 3*b1w)*V**2*v0w + \
            b0w*(324 + 5*b1w*(-28 + 3*b1w))*v0w**3 - \
            2*b0f*(-16 + 3*b1w)*v0f*\
            (5*(-6 + b1f)*V**2 + (-4 + b1f)*v0f**2)))/5.))/256.

    return antider

def intEdV(V0,B0,B0pr,V1,V2):
    """
    integral of E(V) in dV evaluated between volumes V1 and V2
    """
    F1 = antiderE(V0,B0,B0pr,V1)
    F2 = antiderE(V0,B0,B0pr,V2)
    integral = F2 - F1

    return integral

def antiderE(V0,B0,B0pr,V):
    """
    antiderivative of the Birch Murnaghan E(V)
    """
    antider = (9*B0*V0*(-((-6 + B0pr)*V) - ((-4 + B0pr)*V0**2)/V + \
            3*(-14 + 3*B0pr)*V0*(V0/V)**(1/3) + \
            3*(-16 + 3*B0pr)*V*(V0/V)**(2/3)))/16

    return antider

def intE2dV(V0,B0,B0pr,V1,V2):
    """
    Integral of E**2(V) in dV evaluated between volume V1 and volume V2
    """
    F1 = antiderE2(V0,B0,B0pr,V1)
    F2 = antiderE2(V0,B0,B0pr,V2)
    integral = F2 - F1

    return integral

def antiderE2(V0,B0,B0pr,V):
    """
    Antiderivative of the Birch Murnaghan squared (E**2(V))
    """
    antider = (81*B0**2*V0**2*((-6 + B0pr)**2*V + \
            (4*(124 + 5*(-10 + B0pr)*B0pr)*V0**2)/V - \
            ((-4 + B0pr)**2*V0**4)/(3.*V**3) - \
            (3*(V0/V)**(2/3)* \
            (10*(-6 + B0pr)*(-16 + 3*B0pr)*V**2 + \
            (324 + 5*B0pr*(-28 + 3*B0pr))*V0**2))/(5.*V) \
            + (V0/V)**(1/3)* \
            (-3*(424 + 5*B0pr*(-32 + 3*B0pr))*V0 + \
            (6*(-4 + B0pr)*(-14 + 3*B0pr)*V0**3)/(7.*V**2))) \
            )/256.

    return antider

def epsilon2_SSR(v0w, b0w, b1w, v0f, b0f, b1f, config_string):
    """
    Calculate alternative Delta2 based on 2 EOS fits
    THE SIGNATURE OF THIS FUNCTION HAS BEEN CHOSEN TO MATCH THE ONE OF ALL THE OTHER FUNCTIONS
    RETURNING A QUANTITY THAT IS USEFUL FOR COMPARISON, THIS SIMPLIFIES THE CODE LATER.
    Even though 'config_string' is useless here.
    """

    # volume range
    Vi = 0.94 * (v0w + v0f) / 2.
    Vf = 1.06 * (v0w + v0f) / 2.
    deltaV = Vf - Vi

    intdiff2 = intE12sq(v0w,b0w,b1w,v0f,b0f,b1f,Vi,Vf)
    Eavg1 = intEdV(v0w,b0w,b1w,Vi,Vf)/deltaV
    Eavg2 = intEdV(v0w,b0w,b1w,Vi,Vf)/deltaV
    int3 = intE2dV(v0w,b0w,b1w,Vi,Vf) - \
            2*Eavg1*intEdV(v0w,b0w,b1w,Vi,Vf) + \
            deltaV*Eavg1**2 # integrate (ene - mean(ene))**2
    int4 = intE2dV(v0f,b0f,b1f,Vi,Vf) - \
            2*Eavg2*intEdV(v0f,b0f,b1f,Vi,Vf) + \
            deltaV*Eavg2**2
    delta2 = intdiff2/np.sqrt(int3*int4)

    return delta2

def rel_errors_vec_lenght(v0w, b0w, b1w, v0f, b0f, b1f, config_string, weight_b0=1, weight_b1=1, fact=1):
    """
    Returns the lenght of the vector formed by the relative error of V0, B0, B1
    THE SIGNATURE OF THIS FUNCTION HAS BEEN CHOSEN TO MATCH THE ONE OF ALL THE OTHER FUNCTIONS
    RETURNING A QUANTITY THAT IS USEFUL FOR COMPARISON, THIS SIMPLIFIES THE CODE LATER.
    Even though config_string is not usd
    """
    V0err =  2*(v0w-v0f)/(v0w+v0f)
    B0err =  2*(b0w-b0f)/(b0w+b0f)
    B1err =  2*(b1w-b1f)/(b1w+b1f)
    leng = np.sqrt(V0err**2+(weight_b0*B0err)**2+(weight_b1*B1err)**2)
    return leng*fact

def B1_rel_diff(v0w, b0w, b1w, v0f, b0f, b1f, config_string):
    """
    Returns the relative difference in the first derivative of bulk modulus.
    THE SIGNATURE OF THIS FUNCTION HAS BEEN CHOSEN TO MATCH THE ONE OF ALL THE OTHER FUNCTIONS
    RETURNING A QUANTITY THAT IS USEFUL FOR COMPARISON, THIS SIMPLIFIES THE CODE LATER.
    Even though several inputs are useless here.
    """
    return 2*(b1w-b1f)/(b1w+b1f)


def V0_rel_diff(v0w, b0w, b1w, v0f, b0f, b1f, config_string):
    """
    Returns the relative difference in the volumes.
    THE SIGNATURE OF THIS FUNCTION HAS BEEN CHOSEN TO MATCH THE ONE OF ALL THE OTHER FUNCTIONS
    RETURNING A QUANTITY THAT IS USEFUL FOR COMPARISON, THIS SIMPLIFIES THE CODE LATER.
    Even though several inputs are useless here.
    """
    return 2*(v0w-v0f)/(v0w+v0f)


def B0_rel_diff(v0w, b0w, b1w, v0f, b0f, b1f, config_string):
    """
    Returns the relative difference in the bulk modulus.
    THE SIGNATURE OF THIS FUNCTION HAS BEEN CHOSEN TO MATCH THE ONE OF ALL THE OTHER FUNCTIONS
    RETURNING A QUANTITY THAT IS USEFUL FOR COMPARISON, THIS SIMPLIFIES THE CODE LATER.
    Even though several inputs are useless here.
    """
    return 2*(b0w-b0f)/(b0w+b0f)

import functools

quantity_for_comparison_map = {
    "B0_rel_diff": B0_rel_diff,
    "V0_rel_diff": V0_rel_diff,
    "B1_rel_diff": B1_rel_diff,
    #"rel_errors_vec_lenght": functools.partial(rel_errors_vec_lenght,weight_b0=0.1,weight_b1=0.01,fact=100),
    "epsilon2_SSR": epsilon2_SSR
}


if __name__ == "__main__":
    try:
        QUANTITY = sys.argv[1]
    except IndexError:
        print(f"The first argument must be the quantity to use for comparison. Choose among {quantity_for_comparison_map.keys()}")
        sys.exit(1)

    if QUANTITY not in quantity_for_comparison_map.keys():
        print(f"The first argument must be the quantity to use for comparison. Choose among {quantity_for_comparison_map.keys()}")
        sys.exit(1)

    all_args = sys.argv[2:]

    if not all_args:
        print("The plugin's names whose results will be plotted must be listed explicitely as script arguments.")
        sys.exit(1)
    
    try:
        with open(f'results-{PLUGIN_NAME}.json') as fhandle:
            reference_plugin_data = json.load(fhandle)
    except OSError:
        print(f"No data found for your plugin '{PLUGIN_NAME}'. Did you run `./get_results.py` first?")
        sys.exit(1)
    
    print(f"Using data for plugin '{PLUGIN_NAME}' compared with {all_args}.")

    compare_plugin_data = []
    for compare_with in all_args:
        try:
            with open(f'results-{compare_with}.json') as fhandle:
                compare_plugin_data.append(json.load(fhandle))
        except OSError:
            print(f"No data found for the plugin '{compare_with}': you need the file results-{compare_with}.json.")
            sys.exit(1)

    name_file = f'histo-{QUANTITY}-{PLUGIN_NAME}'

    all_systems = set(reference_plugin_data['eos_data'].keys())
    all_systems = set(reference_plugin_data['BM_fit_data'].keys())
    #all_systems.update(compare_plugin_data['BM_fit_data'].keys())

    # Plotting
    fig = pl.figure(figsize=(18,6))

    SMALL_SIZE = 20
    MEDIUM_SIZE = 24
    BIGGER_SIZE = 28

    pl.rc('font', size=SMALL_SIZE)# controls default text sizes
    pl.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
    pl.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    pl.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    pl.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    pl.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    pl.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    
    for index, compare_plugin in enumerate(compare_plugin_data):
        
        collect = []
        collectSmall = []
        collectBig = []

        print(f"comparing with {all_args[index]}")
        progress_bar = tqdm.tqdm(sorted(all_systems))
        for element_and_configuration in progress_bar:
            progress_bar.set_description(f"{element_and_configuration:12s}")
            progress_bar.refresh()

            element, configuration = element_and_configuration.split('-')
            # Get the data for the reference plugin
            ref_BM_fit_data = reference_plugin_data['BM_fit_data'][f'{element}-{configuration}']
        
            if ref_BM_fit_data is None:
                continue
        
            V0=ref_BM_fit_data['min_volume']
            B0=ref_BM_fit_data['bulk_modulus_ev_ang3']
            B01=ref_BM_fit_data['bulk_deriv']

            # Get the data for the compare_with plugin, if specified (and if the EOS worked for the 
            # reference plugin, otherwise we don't know which E0 to use)
            try:
                compare_BM_fit_data = compare_plugin['BM_fit_data'][f'{element}-{configuration}']
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

            CV0=compare_BM_fit_data['min_volume']
            CB0=compare_BM_fit_data['bulk_modulus_ev_ang3']
            CB01=compare_BM_fit_data['bulk_deriv']

            quant = quantity_for_comparison_map[QUANTITY](V0,B0,B01,CV0,CB0,CB01,element_and_configuration)

            if quant < -0.04:
                collectSmall.append(f'{element}-{configuration}')
            if quant > 0.04:
            #if quant > 0.8:
                collectBig.append(f'{element}-{configuration}')

            collect.append(quantity_for_comparison_map[QUANTITY](V0,B0,B01,CV0,CB0,CB01,element_and_configuration))
    
        pl.hist(collect, bins=BINS, range=[-0.04, 0.04], label=f"{all_args[index]}", alpha=0.5)
        #pl.hist(collect, bins=BINS, range=[0.0, 0.6], label=f"{all_args[index]}", alpha=0.5)

        if collectBig:
            pl.annotate(f"{len(collectBig)} more for {all_args[index]}", xy=(pl.xlim()[1], (pl.ylim()[1]-pl.ylim()[0])/2/(index+1)), xytext=(pl.xlim()[1]-0.03, (pl.ylim()[1]-pl.ylim()[0])/2/(index+1)), arrowprops=dict(facecolor='black', shrink=0.05))
        if collectSmall:
            pl.annotate(f"{len(collectSmall)} more for {all_args[index]}", xy=(pl.xlim()[0], (pl.ylim()[1]-pl.ylim()[0])/2/(index+1)), xytext=(pl.xlim()[0]+0.01, (pl.ylim()[1]-pl.ylim()[0])/2/(index+1)), arrowprops=dict(facecolor='black', shrink=0.05))
   
    
    pl.legend(loc='upper right')
    pl.xlabel(f"{QUANTITY}")
    pl.ylabel("Frequency")
    pl.title(f"{PLUGIN_NAME}")
    pl.tight_layout()
    pl.savefig(f"{name_file}")
    pl.close(fig)

    print("ok")
