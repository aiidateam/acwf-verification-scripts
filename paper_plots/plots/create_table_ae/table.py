from pymatgen.core.periodic_table import Element
import json
import paper_plots.quantities_for_comparison as qc 
"""
This script creates a table with the V0, B0, B1 results for WIEN2K, FLEUR and their average
for every cystal structure (the 4 unaries and the 6 oxides)
"""


def loop_over_configurations():
    """
    Main function that reads the data and, for each set, calls the function that creates 
    the table.
    """
    l=[0,0,0]
    with open('out', 'w') as w:
        with open('../../code-data/results-unaries-verification-PBE-v1-fleur.json', 'r') as d:
            fleur = json.load(d)
        with open('../../code-data/results-unaries-verification-PBE-v1-wien2k-dk_0.06.json', 'r') as d:
            wien2k = json.load(d)
        with open('../../code-data/results-unaries-verification-PBE-v1-AE-average.json', 'r') as d:
            av = json.load(d)
        for configuration in ['X/FCC', 'X/BCC', 'X/SC', 'X/Diamond']:
            sett = configuration[2:]
            create_table(w, sett, configuration, fleur, wien2k, av, l)
        with open('../../code-data/results-oxides-verification-PBE-v1-fleur.json', 'r') as d:
            fleur = json.load(d)
        with open('../../code-data/results-oxides-verification-PBE-v1-wien2k-dk_0.06.json', 'r') as d:
            wien2k = json.load(d)
        with open('../../code-data/results-oxides-verification-PBE-v1-AE-average.json', 'r') as d:
            av = json.load(d)
        for configuration in ["X2O", "XO", "X2O3", "XO2", "X2O5", "XO3"]:
            sett = configuration
            create_table(w, sett, configuration, fleur, wien2k, av, l)
        print('V0 more than 0.1% diff, B0 more than 1% diff, B1 more than 2% diff')
        print(l[0], f'({100*l[0]/960}%)', l[1], f'({100*l[1]/960}%)', l[2], f'({100*l[2]/960}%)')


def create_table(w, sett, configuration, fleur, wien2k, av, l):
    """
    Creates the table
    """
    set_name = '{table-ae-'+sett+'}'
    w.write('\\begin{center}\n')
    w.write('\\begin{longtable}{c|c|c|c|c|c|c|c|c|c|c|c|c|}\n')
    w.write(
        '\caption{' +
        f'\label{set_name} The table with all the calculated '
        '\gls{eos} parameters for the '
        f'{sett} structures obtained '
        'with \\fleur{} and \wientwok.} \\\\ \n'
        )
    w.write('& \multicolumn{3}{c|}{ \\fleur } & \multicolumn{3}{c|}{ \wientwok } & \multicolumn{3}{c|}{relative difference} & \multicolumn{3}{c|}{ average }\\\\\n')
    w.write('& V0 & B0  & B1 & V0 & B0 &  B1 & \% V0 & \% B0 &  \%  B1 & V0 & B0  & B1\\\\\n' )
    w.write('\hline\n')
    for i in range(1,96):
        element = Element.from_Z(i).name
        ref_BM_fit_data = fleur['BM_fit_data'][f'{element}-{configuration}']
        scaling_factor_f = qc.get_volume_scaling_to_formula_unit(
            fleur['num_atoms_in_sim_cell'][f'{element}-{configuration}'],
            element, configuration
        )
        V0=ref_BM_fit_data['min_volume']/scaling_factor_f
        B0=ref_BM_fit_data['bulk_modulus_ev_ang3']
        B01=ref_BM_fit_data['bulk_deriv']
        ref_BM_fit_data_2 = wien2k['BM_fit_data'][f'{element}-{configuration}']
        scaling_factor_w = qc.get_volume_scaling_to_formula_unit(
            wien2k['num_atoms_in_sim_cell'][f'{element}-{configuration}'],
            element, configuration
        )
        V0w=ref_BM_fit_data_2['min_volume']/scaling_factor_w
        B0w=ref_BM_fit_data_2['bulk_modulus_ev_ang3']
        B01w=ref_BM_fit_data_2['bulk_deriv']
        ref_BM_fit_data_3 = av['BM_fit_data'][f'{element}-{configuration}']
        scaling_factor_av = qc.get_volume_scaling_to_formula_unit(
            av['num_atoms_in_sim_cell'][f'{element}-{configuration}'],
            element, configuration
        )
        relative_diff_V0 = round(2*100*abs(V0-V0w)/(V0+V0w), 3)
        if relative_diff_V0 > 0.1:
            l[0]=l[0]+1
        relative_diff_B0 = round(2*100*abs(B0-B0w)/(B0+B0w), 3)
        if relative_diff_B0 > 1:
            l[1]=l[1]+1
        relative_diff_B1 = round(2*100*abs(B01-B01w)/(B01+B01w), 3)
        if relative_diff_B1 > 2:
            l[2]=l[2]+1

        V0=round(ref_BM_fit_data['min_volume']/scaling_factor_f, 4)
        B0=round(ref_BM_fit_data['bulk_modulus_ev_ang3'], 4)
        B01=round(ref_BM_fit_data['bulk_deriv'], 4)

        V0w=round(ref_BM_fit_data_2['min_volume']/scaling_factor_w, 4)
        B0w=round(ref_BM_fit_data_2['bulk_modulus_ev_ang3'], 4)
        B01w=round(ref_BM_fit_data_2['bulk_deriv'], 4)

        V0av=round(ref_BM_fit_data_3['min_volume']/scaling_factor_av, 4)
        B0av=round(ref_BM_fit_data_3['bulk_modulus_ev_ang3'], 4)
        B01av=round(ref_BM_fit_data_3['bulk_deriv'], 4)
    
        w.write(f'{element} & {V0} & {B0}  & {B01} & {V0w} & {B0w} & {B01w}& {relative_diff_V0} \% & {relative_diff_B0} \% &  {relative_diff_B1} \% & {V0av} & {B0av} & {B01av}\\\\\n')
    
    w.write('\end{longtable}\n')
    
    w.write('\end{center}\n')


if __name__ == "__main__":
    loop_over_configurations()
