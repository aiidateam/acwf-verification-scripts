import numpy as np
import pylab as pl
from quantities_for_comparison import birch_murnaghan, nu, epsilon, delta

pl.rc('xtick', labelsize=18)    # fontsize of the tick labels
pl.rc('ytick', labelsize=18)

fig, ax = pl.subplots(1, 4, figsize=(18,6), sharey=True)
fig.ylabel = 'energy'

# A generic material, with V0, B0, B1 that is average among the set we consider
vol_form_unit = 50.61 #48.64492911765141
B0 = 0.71 #0.7236923751242582
B1 = 4.67 #4.426027837689244

#                             V0    B0  B1
relative_errors_analyzed = [(0.01, 0.1, 1), (2, 0, 0), (0, 2, 0), (0, 0, 2)] 

for index, value in enumerate(relative_errors_analyzed):
    
    a_vol_form_unit = vol_form_unit+vol_form_unit*value[0]/100
    a_B0 = B0+B0*value[1]/100
    a_B1 = B1+B1*value[2]/100
    
    dense_volumes = np.linspace(
          0.94*(vol_form_unit+a_vol_form_unit)/2.,
          1.06*(vol_form_unit+a_vol_form_unit)/2.,
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
    
    if index==3:
        ax[index].plot(dense_volumes, eos_ref, '-b', label='A')
        ax[index].plot(dense_volumes, eos, '-r', label='B')
        ax[index].legend(loc='upper right', fontsize=20, ncol=2)
    else:
        ax[index].plot(dense_volumes, eos_ref, '-b')
        ax[index].plot(dense_volumes, eos, '-r')
    if index==0:
        ax[index].set_ylabel("Energy",fontsize=22)
    ax[index].set_xlabel("Volume",fontsize=22)
    if index==1:
        ax[index].annotate(f'{index+1})', [48.1, 0.098], fontsize=23)
    else:
        ax[index].annotate(f'{index+1})', [47.7, 0.098], fontsize=23)

    eps = epsilon(a_vol_form_unit, a_B0, a_B1, vol_form_unit, B0, B1, 1, 0, 0)
    delt = delta(a_vol_form_unit, a_B0, a_B1, vol_form_unit, B0, B1, 1, 1, 1)
    v2 = nu(a_vol_form_unit, a_B0, a_B1, vol_form_unit, B0, B1, 1, 1/6, 1/35)

    if index==1:
        #ww=eps/100
        ax[index].annotate(r'$\varepsilon$ = '+f'{eps:.2E}', [49.2, 0.08], fontsize=23)
        ax[index].annotate(r'$\Delta$ = '+f'{delt:.2E}', [49.2, 0.07], fontsize=23)
        ax[index].annotate(r'$\nu_2$ = '+f'{v2:.2E}', [49.2, 0.06], fontsize=23)
    else:
        #ax[index].annotate(r'$\varepsilon$ = '+f'{eps.round(2)}'+r'$\cdot$10$^{-2}$', [48.6, 0.08], fontsize=23)
        ax[index].annotate(r'$\varepsilon$ = '+f'{eps:.2E}', [48.6, 0.08], fontsize=23)
        ax[index].annotate(r'$\Delta$ = '+f'{delt:.2E}', [48.6, 0.07], fontsize=23)
        ax[index].annotate(r'$\nu_2$ = '+f'{v2:.2E}', [48.6, 0.06], fontsize=23)

fig.tight_layout()
pl.savefig('Sensibility_EoSes.pdf')
pl.close(fig)
