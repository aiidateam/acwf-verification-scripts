########################################################
# Comparison figure for the case of Ba psp with Abinit #
########################################################
 
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.gridspec as gridspec
import json
import sys
from quantities_for_comparison import birch_murnaghan, get_volume_scaling_to_formula_unit

fig = plt.figure(figsize=[7.0,4.0],dpi=300)

size = 8
lw = 1
lwv = 0.5
lg = 2
ps = 20
ms = 10
props = dict(boxstyle='round', facecolor='white', edgecolor='white', alpha=1.0, pad=0.0)

col0 = '#808080' # Gray
col1 = '#70CA2C'
col11 = '#41ab5d'
col12 = '#d9f0a3'
col13 = '#bf812d'
col2 = '#DB8DEF'
col21 = '#dd3497'
col22 = '#fa9fb5'
col23 = '#8073ac'
col3 = '#feb24c'
col31 = '#ec7014'
col32 = '#fed976'
col33 = '#b2182b'
col4 = 'black'

gray1 = '#808080'
gray2 = '#979797'
gray3 = '#696969'

gs1 = gridspec.GridSpec(2, 2)
gs1.update(wspace=0.2, hspace=0.25, left=0.12,right=0.95,bottom=0.1,top=0.95)

##############################################################################
##############################################################################
##############################################################################
# Load files
with open('Ba-sp-v0.4-nofproj.json', 'r', encoding='utf-8') as f:
 data04 = json.load(f)
# Energy mesh in Hartree
energies = data04["energies_Ha"]

with open('Ba-sp-v0.5.fproj.json', 'r', encoding='utf-8') as g:
 data05 = json.load(g)

#####################################################

ax0 = fig.add_subplot(gs1[0,0])

shift = np.ones(len(energies))

ax0.plot(energies[:],  data04["s"]["ae_arctan_logder"],   linestyle='-' , linewidth=lw, color=col1, zorder=1, label='s-AE')
ax0.plot(energies[:],  data04["s"]["ps_arctan_logder"],   linestyle='--' , linewidth=lw+1, color=col1, zorder=1, label='s-PS')
ax0.plot(energies[:],  data04["p"]["ae_arctan_logder"] + shift*4.0,   linestyle='-' , linewidth=lw, color=col2, zorder=1, label='p-AE')
ax0.plot(energies[:],  data04["p"]["ps_arctan_logder"] + shift*4.0,   linestyle='--' , linewidth=lw+1, color=col2, zorder=1, label='p-PS')
ax0.plot(energies[:],  data04["d"]["ae_arctan_logder"] + shift*3.0,   linestyle='-' , linewidth=lw, color=col3, zorder=1, label='d-AE')
ax0.plot(energies[:],  data04["d"]["ps_arctan_logder"] + shift*3.0,   linestyle='--' , linewidth=lw+1, color=col3, zorder=1, label='d-PS')
ax0.plot(energies[:],  data04["f"]["ae_arctan_logder"] + shift*7.0,   linestyle='-' , linewidth=lw, color=col4, zorder=1, label='f-AE')
ax0.plot(energies[:],  data04["f"]["ps_arctan_logder"] + shift*7.0,   linestyle='--' , linewidth=lw+1, color=col4, zorder=1, label='f-PS')

ax0.xaxis.set_tick_params(which='major', width=lw, length=lg, labelsize=size, direction='in')
ax0.yaxis.set_tick_params(which='major', width=lw, length=lg, labelsize=size, direction='in')
ax0.xaxis.set_ticks_position('both')
ax0.yaxis.set_ticks_position('both')

plt.ylabel(r'$atan(\frac{R}{\Psi(E,R)} \frac{\partial \Psi(E,r)}{\partial r}|_R)$ [-]',fontsize=size)
plt.title('Ba - PseudoDojo v0.4',fontsize=size)

leg = plt.legend(loc=2, fontsize=size, frameon=True, handlelength=1.0, framealpha=1.0, fancybox=True, borderpad=0.3)
leg.get_frame().set_linewidth(0.0)

plt.axis([-5, 3, 2, 15])
##############################################################################
##############################################################################
##############################################################################
ax1 = fig.add_subplot(gs1[1,0])

ax1.plot(energies,  data05["s"]["ae_arctan_logder"],   linestyle='-' , linewidth=lw, color=col1, zorder=1, label='s-AE')
ax1.plot(energies,  data05["s"]["ps_arctan_logder"],   linestyle='--' , linewidth=lw+1, color=col1, zorder=1, label='s-PS')
ax1.plot(energies,  data05["p"]["ae_arctan_logder"] + shift*4.0,   linestyle='-' , linewidth=lw, color=col2, zorder=1, label='p-AE')
ax1.plot(energies,  data05["p"]["ps_arctan_logder"] + shift*4.0,   linestyle='--' , linewidth=lw+1, color=col2, zorder=1, label='p-PS')
ax1.plot(energies,  data05["d"]["ae_arctan_logder"] + shift*3.0,   linestyle='-' , linewidth=lw, color=col3, zorder=1, label='d-AE')
ax1.plot(energies,  data05["d"]["ps_arctan_logder"] + shift*3.0,   linestyle='--' , linewidth=lw+1, color=col3, zorder=1, label='d-PS')
ax1.plot(energies,  data05["f"]["ae_arctan_logder"] + shift*7.0,   linestyle='-' , linewidth=lw, color=col4, zorder=1, label='f-AE')
ax1.plot(energies,  data05["f"]["ps_arctan_logder"] + shift*7.0,   linestyle='--' , linewidth=lw+1, color=col4, zorder=1, label='f-PS')

ax1.xaxis.set_tick_params(which='major', width=lw, length=lg, labelsize=size, direction='in')
ax1.yaxis.set_tick_params(which='major', width=lw, length=lg, labelsize=size, direction='in')
ax1.xaxis.set_ticks_position('both')
ax1.yaxis.set_ticks_position('both')

plt.xlabel('Energy [Ha] ',fontsize=size)
plt.ylabel(r'$atan(\frac{R}{\Psi(E,R)} \frac{\partial \Psi(E,r)}{\partial r}|_R)$ [-] ',fontsize=size)
plt.title('Ba - PseudoDojo v0.5b',fontsize=size)

plt.axis([-5, 3, 2, 15])

##############################################################################
##############################################################################
##############################################################################


ax3 = fig.add_subplot(gs1[0,1])

# Ba - BCC which is the stable phase - reference AE average
with open('results-unaries-verification-PBE-v1-AE-average.json', 'r', encoding='utf-8') as f:
  Baae = json.load(f)
# Energy mesh in Hartree
BM_fit_ref = Baae["BM_fit_data"]["Ba-X/BCC"]
natom = Baae["num_atoms_in_sim_cell"]["Ba-X/BCC"]

# Ba - PseudoDojo v0.4
with open('results-unaries-verification-PBE-v1-abinit-PseudoDojo-0.4-PBE-SR-standard-psp8.json', 'r', encoding='utf-8') as f:
  Ba04 = json.load(f)
# Energy mesh in Hartree
BM_fit_04 = Ba04["BM_fit_data"]["Ba-X/BCC"]
eos_data04 = Ba04["eos_data"]["Ba-X/BCC"]
natom04 = Ba04["num_atoms_in_sim_cell"]["Ba-X/BCC"]

# Ba - PseudoDojo v0.5
with open('results-unaries-verification-PBE-v1-abinit-PseudoDojo-0.5b1-PBE-SR-standard-psp8.json', 'r', encoding='utf-8') as f:
  Ba05 = json.load(f)
# Energy mesh in Hartree
BM_fit_05 = Ba05["BM_fit_data"]["Ba-X/BCC"]
eos_data05 = Ba05["eos_data"]["Ba-X/BCC"]

scaling_ref_plugin = get_volume_scaling_to_formula_unit(natom04, "Ba", "X/BCC")
volumes04, energies04 = (np.array(eos_data04).T / scaling_ref_plugin).tolist()
volumes05, energies05 = (np.array(eos_data05).T / scaling_ref_plugin).tolist()
dense_volumes = np.linspace(
    min(volumes04),
    max(volumes04),
    100
)

fitted_free_energy_ref = birch_murnaghan(
     V=dense_volumes,
     E0 =BM_fit_ref['E0'],
     V0 =BM_fit_ref['min_volume'],
     B0 =BM_fit_ref['bulk_modulus_ev_ang3'],
     B01=BM_fit_ref['bulk_deriv']
)

fitted_free_energy_04 = birch_murnaghan(
     V=dense_volumes,
     E0 =BM_fit_ref['E0'] / scaling_ref_plugin, # We use the REF E0 !! 
     V0 =BM_fit_04['min_volume'] / scaling_ref_plugin,
     B0 =BM_fit_04['bulk_modulus_ev_ang3'],
     B01=BM_fit_04['bulk_deriv']
)

fitted_free_energy_05 = birch_murnaghan(
     V=dense_volumes,
     E0 =BM_fit_ref['E0'] / scaling_ref_plugin, # We use the REF E0 !!
     V0 =BM_fit_05['min_volume'] / scaling_ref_plugin,
     B0 =BM_fit_05['bulk_modulus_ev_ang3'],
     B01=BM_fit_05['bulk_deriv']
)

E0_04 = BM_fit_04['E0']
E0_05 = BM_fit_05['E0']
tmp = np.ones(len(energies04))

ax3.plot(dense_volumes,  fitted_free_energy_ref*1000,   linestyle='-' , linewidth=lw, color=col4, zorder=1, label='AE reference')
ax3.plot(dense_volumes,  fitted_free_energy_04*1000,   linestyle='-' , linewidth=lw, color=col1, zorder=1, label='PseudoDojo v0.4')
ax3.scatter(volumes04,  1000*(energies04 - E0_04*tmp),   s=ps, marker='o', edgecolor=col1, color='white', zorder=12, clip_on=True)
ax3.plot(dense_volumes,  fitted_free_energy_05*1000,   linestyle='-' , linewidth=lw, color=col3, zorder=1, label='PseudoDojo v0.5b')
ax3.scatter(volumes05,  1000*(energies05 - E0_05*tmp),   s=ps, marker='o', edgecolor=col3, color='white', zorder=12, clip_on=False)

ax3.xaxis.set_tick_params(which='major', width=lw, length=lg, labelsize=size, direction='in')
ax3.yaxis.set_tick_params(which='major', width=lw, length=lg, labelsize=size, direction='in')
ax3.xaxis.set_ticks_position('both')
ax3.yaxis.set_ticks_position('both')

plt.ylabel(r'E - TS (meV)',fontsize=size)
plt.title('Ba - BCC',fontsize=size)

leg = plt.legend(loc=0, fontsize=size, frameon=True, handlelength=1.0, framealpha=1.0, fancybox=True, borderpad=0.3)
leg.get_frame().set_linewidth(0.0)

plt.xticks([59,61,63,65,67],fontsize=size)
plt.axis([59, 67.3, 0, 8])

##############################################################################
##############################################################################
##############################################################################
ax2 = fig.add_subplot(gs1[1,1])

# BaO3 - reference AE average
with open('results-oxides-verification-PBE-v1-AE-average.json', 'r', encoding='utf-8') as f:
  BaO3ae = json.load(f)
# Energy mesh in Hartree
BM_fit_ref = BaO3ae["BM_fit_data"]["Ba-XO3"]
natom = BaO3ae["num_atoms_in_sim_cell"]["Ba-XO3"]


# BaO3 - PseudoDojo v0.4
with open('results-oxides-verification-PBE-v1-abinit-PseudoDojo-0.4-PBE-SR-standard-psp8.json', 'r', encoding='utf-8') as f:
  BaO304 = json.load(f)
# Energy mesh in Hartree
BM_fit_04 = BaO304["BM_fit_data"]["Ba-XO3"]
eos_data04 = BaO304["eos_data"]["Ba-XO3"]
natom04 = BaO304["num_atoms_in_sim_cell"]["Ba-XO3"]


# BaO3 - PseudoDojo v0.5
with open('results-oxides-verification-PBE-v1-abinit-PseudoDojo-0.5b1-PBE-SR-standard-psp8.json', 'r', encoding='utf-8') as f:
  BaO305 = json.load(f)
# Energy mesh in Hartree
BM_fit_05 = BaO305["BM_fit_data"]["Ba-XO3"]
eos_data05 = BaO305["eos_data"]["Ba-XO3"]


scaling_ref_plugin = get_volume_scaling_to_formula_unit(natom04, "Ba", "XO3")
volumes04, energies04 = (np.array(eos_data04).T / scaling_ref_plugin).tolist()
volumes05, energies05 = (np.array(eos_data05).T / scaling_ref_plugin).tolist()
dense_volumes = np.linspace(
    min(volumes04),
    max(volumes04),
    100
)

fitted_free_energy_ref = birch_murnaghan(
     V=dense_volumes,
     E0 =BM_fit_ref['E0'],
     V0 =BM_fit_ref['min_volume'],
     B0 =BM_fit_ref['bulk_modulus_ev_ang3'],
     B01=BM_fit_ref['bulk_deriv']
)

fitted_free_energy_04 = birch_murnaghan(
     V=dense_volumes,
     E0 =BM_fit_ref['E0'] / scaling_ref_plugin, # We use the REF E0 !! 
     V0 =BM_fit_04['min_volume'] / scaling_ref_plugin,
     B0 =BM_fit_04['bulk_modulus_ev_ang3'],
     B01=BM_fit_04['bulk_deriv']
)

fitted_free_energy_05 = birch_murnaghan(
     V=dense_volumes,
     E0 =BM_fit_ref['E0'] / scaling_ref_plugin, # We use the REF E0 !!
     V0 =BM_fit_05['min_volume'] / scaling_ref_plugin,
     B0 =BM_fit_05['bulk_modulus_ev_ang3'],
     B01=BM_fit_05['bulk_deriv']
)

E0_04 = BM_fit_04['E0']
E0_05 = BM_fit_05['E0']
tmp = np.ones(len(energies04))

ax2.plot(dense_volumes,  fitted_free_energy_ref*1000,   linestyle='-' , linewidth=lw, color=col4, zorder=1, label='AE reference')
ax2.plot(dense_volumes,  fitted_free_energy_04*1000,   linestyle='-' , linewidth=lw, color=col1, zorder=1, label='PseudoDojo v0.4')
ax2.scatter(volumes04,  (energies04 - E0_04*tmp)*1000,   s=ps, marker='o', edgecolor=col1, color='white', zorder=12, clip_on=True)
ax2.plot(dense_volumes,  fitted_free_energy_05*1000,   linestyle='-' , linewidth=lw, color=col3, zorder=1, label='PseudoDojo v0.5b')
ax2.scatter(volumes05,  1000*(energies05 - E0_05*tmp),   s=ps, marker='o', edgecolor=col3, color='white', zorder=12, clip_on=False)

ax2.xaxis.set_tick_params(which='major', width=lw, length=lg, labelsize=size, direction='in')
ax2.yaxis.set_tick_params(which='major', width=lw, length=lg, labelsize=size, direction='in')
ax2.xaxis.set_ticks_position('both')
ax2.yaxis.set_ticks_position('both')

plt.ylabel(r'E - TS (meV)',fontsize=size)
plt.xlabel(r'Cell volume per formula unit ($\AA^3$)',fontsize=size)
plt.title('BaO$_3$',fontsize=size)

plt.xticks([86,88,90,92,94,96,98],fontsize=size)
plt.axis([86, 98, 0, 50])

ax0.text(-0.08,0.99, '(a)', transform=ax0.transAxes, ha="center", va="center", fontsize=size, bbox=props, zorder=100)
ax1.text(-0.08,0.99, '(b)', transform=ax1.transAxes, ha="center", va="center", fontsize=size, bbox=props, zorder=100)
ax3.text(-0.12,0.99, '(c)', transform=ax3.transAxes, ha="center", va="center", fontsize=size, bbox=props, zorder=100)
ax2.text(-0.12,0.99, '(d)', transform=ax2.transAxes, ha="center", va="center", fontsize=size, bbox=props, zorder=100)

plt.setp(ax0.get_xticklabels(), visible=False)
plt.setp(ax0.get_yticklabels(), visible=False)
plt.setp(ax1.get_yticklabels(), visible=False)

plt.savefig('Ba-lodger.pdf')
