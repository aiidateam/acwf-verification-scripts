#!/usr/bin/env python
import json
import os
import sys

import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
import csv

# Plotting
fig = pl.figure(figsize=(9,6))

TINY_SIZE = 15
SMALL_SIZE = 19
MEDIUM_SIZE = 23
BIGGER_SIZE = 27

pl.rc('font', size=SMALL_SIZE)# controls default text sizes
pl.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
pl.rc('axes', labelsize=21)    # fontsize of the x and y labels
pl.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
pl.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
pl.rc('legend', fontsize=TINY_SIZE)    # legend fontsize
pl.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


save_row = []
with open('kp_smearing_data.csv') as f:
    content=csv.reader(f)

    for row in content:
        save_row.append(row)

kp_values = np.array(save_row[0][1:])
values = np.array(save_row[1:]).T.astype('float')

labels_dict = {
    1: '80$^3$',
    2: '70$^3$',
    3: '60$^3$',
    4: '54$^3$',
    5: '48$^3$',
    6: '42$^3$',
    7: '36$^3$',
    8: '30$^3$',
    9: '24$^3$',
    10: '18$^3$'
}

#pl.plot(values[0], values[1])
#pl.plot(values[0], values[2])
#pl.plot(values[0], values[3])
#pl.plot(values[0], values[4])
#pl.plot(values[0], values[5])
#pl.plot(values[0], values[6])
#pl.plot(values[0], values[7])
#pl.plot(values[0], values[8])
#pl.plot(values[0], values[9])
#pl.plot(values[0], values[10])

for i in range(len(kp_values)):
    pl.plot(values[0], values[10-i], '-o', label=labels_dict[10-i])

#pl.show()
# Reset the xlim
#pl.xlim(x[0], x[-1])

pl.legend(loc='lower right', ncol=2)
pl.xlabel("Smearing temperature (meV)")
pl.ylabel("Total forces on displaced atom (ev/$\AA$)")
pl.title(f"Siesta Fermi-Dirac     0.0045 Ry = 61 meV")
pl.ylim(0.463, 0.475)
pl.tight_layout()
pl.savefig(f"smearing_kp_siesta.pdf")
pl.close(fig)

#with open(f"discrepancies-{what}-{plugin1}-VS-{plugin2}.txt", "w") as fhandle:
#    fhandle.write(f"## {plugin1} VS {plugin2}\n")
#    if what == 'formation-energy':
#        fhandle.write(f"## Cases with abs(formation energies) > {PRINT_THRESHOLD} eV/atom:\n")
#    elif what == 'unaries':
#        fhandle.write(f"## Cases with abs(energy difference) > {PRINT_THRESHOLD} eV/atom:\n")
#    else:
#        raise ValueError(f"Unknown value of 'what': '{what}'")
#    for dissimilarity, system, data_plugin1, data_plugin2 in abs_dissimilarities:
#        if dissimilarity < PRINT_THRESHOLD:
#            # Here I assume I already sorted them
#            break
#        fhandle.write(f"{system:30s}: {dissimilarity:.6f} ({data_plugin1} vs {data_plugin2})\n")

