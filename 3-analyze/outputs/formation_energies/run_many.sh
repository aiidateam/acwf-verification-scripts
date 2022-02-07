for i in fleur_testPrecise_16,wien2k-finite-T_delta_k_0.1_prec3-6errors,1. \
    wien2k-finite-T_delta-k_0.045_prec3,wien2k-finite-T_delta_k_0.1_prec3-6errors,0.1 \
    wien2k-finite-T_delta-k_0.045_prec3,wien2k-finite-T_delta-k_0.06_prec3,0.1 \
    abinit-NC,wien2k-finite-T_delta-k_0.06_prec3,1. \
    castep,wien2k-finite-T_delta-k_0.06_prec3,1. \
    quantum_espresso-sssp-k0.1,wien2k-finite-T_delta-k_0.06_prec3,1. \
    vasp,wien2k-finite-T_delta-k_0.06_prec3,1. \
    abinit-NC,fleur_testPrecise_16,0.1 \
    castep,fleur_testPrecise_16,0.1 \
    quantum_espresso-sssp-k0.1,fleur_testPrecise_16,0.1 \
    vasp,fleur_testPrecise_16,0.1
do
    plugin1=`echo $i | cut -f1 -d,`
    plugin2=`echo $i | cut -f2 -d,`
    zoom=`echo $i | cut -f3 -d,`

    ./compute_formation_energies.py $plugin1 $plugin2
    ./plot_histo_formation_energies.py $plugin1 $plugin2 $zoom
done