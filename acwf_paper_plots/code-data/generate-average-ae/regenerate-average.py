import json

with open(f'../results-unaries-verification-PBE-v1-fleur.json') as fhandle:
    fleur_un = json.load(fhandle)

with open(f'../results-unaries-verification-PBE-v1-wien2k-dk_0.06.json') as fhandle:
    wien2k_un = json.load(fhandle)

with open(f'../results-oxides-verification-PBE-v1-fleur.json') as fhandle:
    fleur_ox = json.load(fhandle)

with open(f'../results-oxides-verification-PBE-v1-wien2k-dk_0.06.json') as fhandle:
    wien2k_ox = json.load(fhandle)


for set_name in ['unaries','oxides']:

    coll = {'BM_fit_data':{}}

    if set_name == 'unaries':
        all_systems = set(fleur_un['BM_fit_data'].keys())
        ref_BM_fit_data = fleur_un['BM_fit_data']
        wk_BM_fit_data = wien2k_un['BM_fit_data']
        num_at_w = wien2k_un['num_atoms_in_sim_cell']
        num_at_f = fleur_un['num_atoms_in_sim_cell']
    else:
        all_systems = set(fleur_ox['BM_fit_data'].keys())
        ref_BM_fit_data = fleur_ox['BM_fit_data']
        wk_BM_fit_data = wien2k_ox['BM_fit_data']
        num_at_w = wien2k_ox['num_atoms_in_sim_cell']
        num_at_f = fleur_ox['num_atoms_in_sim_cell']

    for element_and_configuration in all_systems:
    
        element, configuration = element_and_configuration.split('-')
        
        #ref_BM_fit_data[f'{element}-{configuration}']
        #wk_BM_fit_data[f'{element}-{configuration}']
        #num_at_w[f'{element}-{configuration}']
        #num_at_f[f'{element}-{configuration}']
        
        if num_at_w[f'{element}-{configuration}'] != num_at_f[f'{element}-{configuration}']:
            print(f'Error! Number of atoms for {element}-{configuration} is different between fleur and wien2k. Change the script.')

        V0=ref_BM_fit_data[f'{element}-{configuration}']['min_volume']
        B0=ref_BM_fit_data[f'{element}-{configuration}']['bulk_modulus_ev_ang3']
        B01=ref_BM_fit_data[f'{element}-{configuration}']['bulk_deriv']
    
        V1=wk_BM_fit_data[f'{element}-{configuration}']['min_volume']
        B1=wk_BM_fit_data[f'{element}-{configuration}']['bulk_modulus_ev_ang3']
        B11=wk_BM_fit_data[f'{element}-{configuration}']['bulk_deriv']
    
        V=(V0+V1)/2
        B=(B0+B1)/2
        Br=(B01+B11)/2
    
        if abs((B0-B1)/B) > 0.02:
            print(f"B0 difference > 2% for {element_and_configuration}")
        
        
        coll['BM_fit_data'][f'{element}-{configuration}'] = {
                'E0': 0,
                'bulk_deriv': Br,
                'bulk_modulus_ev_ang3': B,
                'min_volume': V,
                'residuals': 0 
                }
        coll['num_atoms_in_sim_cell'] = num_at_f
        coll['script_version'] = '0.0.3'
        # Keys that I do not add: 'completely_off', 'eos_data', 'failed_wfs', 'missing_outputs', 'set_name', 'stress_data', 'uuid_mapping'])
    
    with open(f'../results-{set_name}-verification-PBE-v1-AE-average.json', 'w') as dumpf:
        json.dump(coll, dumpf, indent=2, sort_keys=True)
