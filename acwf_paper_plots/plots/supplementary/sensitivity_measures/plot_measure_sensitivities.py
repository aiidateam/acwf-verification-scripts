#!/usr/bin/env python
import numpy as np
import pylab as pl
import string
from acwf_paper_plots.quantities_for_comparison import birch_murnaghan, nu, epsilon, delta

def beautify_exp(float_value):
    """Typesets a floa value as a nice LaTeX expression in exponential form."""
    stringEformat = f'{float_value:.2E}'

    assert "E" in stringEformat, "No 'E' in string, are you sure it's a float number typeset witht eh {:.2E} format?"
    assert len([c for c in stringEformat if c=="E"]) == 1, "More than one 'E' in string, are you sure it's a float number typeset witht eh {:.2E} format?"

    base, full_exponent = stringEformat.split("E")
    exp_sign, exp_value = full_exponent[0], full_exponent[1:]
    assert exp_sign in ["+", "-"], "Exponent sign (after E) neither + nor -!!"

    if exp_sign == "+": 
        exp_sign = ""
    
    while exp_value.startswith("0"):
        exp_value = exp_value[1:]

    if exp_value:
        exponent = f"{exp_sign}{exp_value}"
    else:
        # Nothing remained, so it's zero
        exponent = ""

    if exponent:
        retstring = fr"${base}\cdot 10^{{{exponent}}}$"
    else:
        retstring = fr"${base}$"
    return retstring

def beautify_flat(float_value):
    """Output the number without a scientific notation, with fixed number of digits."""
    #print(float_value)
    return f'{float_value:.2f}'

## CHANGE HERE IF YOU WANT A SCIENTIFIC NOTATION
#beautify = beautify_exp
beautify = beautify_flat

fontsize = 22

pl.rc('xtick', labelsize=fontsize)    # fontsize of the tick labels
pl.rc('ytick', labelsize=fontsize)

# A generic material, with V0, B0, B1 that is average among the set we consider
vol_form_unit = 50.61 #48.64492911765141
B0 = 0.71 #0.7236923751242582
B1 = 4.67 #4.426027837689244


if __name__ == "__main__":

    ### FIRST FIGURE, small changes
    #                            V0    B0  B1
    #relative_errors_analyzed = [(0.01, 0.1, 1), (2, 0, 0), (0, 2, 0), (0, 0, 2)] 
    ref_error_V0 = 0.4
    relative_errors_analyzed_small = [(0.06*2, 0.35*2, 2*2), (ref_error_V0, 0, 0), (0, ref_error_V0*20, 0), (0, 0, ref_error_V0*400)] 

    ### SECON FIGURE, large changes
    relative_errors_analyzed_large = [(1.8, 0, 0), (0, 143, 0), (0, -50.5, 0), (3.52, 0, 0)]

    for FILENAME, relative_errors_analyzed in [
            ('measures_sensitivity.pdf', relative_errors_analyzed_small),
            ('measures_sensitivity_large.pdf', relative_errors_analyzed_large),
        ]:
        NUM_SUBPLOTS = len(relative_errors_analyzed)
        fig, ax = pl.subplots(1, NUM_SUBPLOTS, figsize=(4.5 * NUM_SUBPLOTS,6), sharey=True)
        for index, value in enumerate(relative_errors_analyzed):
            # Perturbed values
            a_vol_form_unit = vol_form_unit+vol_form_unit*value[0]/100
            a_B0 = B0+B0*value[1]/100
            a_B1 = B1+B1*value[2]/100
            
            dense_volumes = np.linspace(
                0.94*vol_form_unit,
                1.06*vol_form_unit,
                100
            )
            
            eos_ref = birch_murnaghan(
                V=dense_volumes,
                E0=0, ## IMPORTANT! here we use the E0 of the reference plugin
                V0=vol_form_unit,
                B0=B0,
                B01=B1
            )
            
            eos = birch_murnaghan(
                V=dense_volumes,
                E0=0, ## IMPORTANT! here we use the E0 of the reference plugin
                V0=a_vol_form_unit,
                B0=a_B0,
                B01=a_B1
            )
            
            ax[index].plot(dense_volumes, eos_ref, '--b', linewidth=2, label='Reference')
            ax[index].plot(dense_volumes, eos, '-r', linewidth=1, label='Perturbed')
            if index==NUM_SUBPLOTS - 1:
                ax[index].legend(loc='upper right', fontsize=fontsize*0.8, ncol=1)
            if index==0:
                ax[index].set_ylabel("Energy [eV]",fontsize=fontsize)
            ax[index].set_xlabel(r"Volume [$\mathrm{\AA}^3$]",fontsize=fontsize)
            ax[index].annotate(f'({string.ascii_lowercase[index]})', [ax[index].get_xlim()[0] + 0.4, 0.103], fontsize=fontsize)

            eps = epsilon(a_vol_form_unit, a_B0, a_B1, vol_form_unit, B0, B1, 1, 0, 0)
            delt = delta(a_vol_form_unit, a_B0, a_B1, vol_form_unit, B0, B1, 1, 0, 0)
            v2 = nu(a_vol_form_unit, a_B0, a_B1, vol_form_unit, B0, B1, 100, 1/20, 1/400)
            #print(a_vol_form_unit, a_B0, a_B1, "|", vol_form_unit, B0, B1)

            xpos = ax[index].get_xlim()[0] + 0.5
            ypos = 0.095
            err_info = []
            if value[0]:
                err_info.append(rf"$V_0$: {value[0]}%")
            if value[1]:
                err_info.append(rf"$B_0$: {value[1]}%")
            if value[2]:
                err_info.append(rf"$B_1$: {value[2]}%")
            ax[index].annotate("\n".join(err_info), [xpos, ypos], fontsize=int(fontsize*0.8), verticalalignment='top')

            xpos = ax[index].get_xlim()[0] + 1
            ypos = 0.07
            ax[index].annotate(r'$\varepsilon$ = '+beautify(eps), [xpos, ypos], fontsize=fontsize)
            ax[index].annotate(r'$\nu$ = '+beautify(v2), [xpos, ypos-0.01], fontsize=fontsize)
            ax[index].annotate(r'$\Delta$ = '+beautify(delt) + " meV", [xpos, ypos-0.02], fontsize=fontsize)

        ax[1].set_ylim(0, 0.11)

        fig.tight_layout()
        pl.savefig(FILENAME)
        pl.close(fig)
