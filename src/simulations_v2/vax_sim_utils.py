from multi_group_simulation import MultiGroupSimulation
import numpy as np
import yaml
from load_params import load_params, load_age_sev_params

def generate_vax_unvax_multigroup_sim(orig_group_params, orig_group_names,
                                    vax_rates_by_group, orig_contact_matrix,
                                    vax_trans_mult, vax_susc_mult, 
                                    vax_age_sev_dist_files=None):
    N = len(orig_group_params)
    assert(len(orig_group_names) == N)
    assert(len(vax_rates_by_group) == N)
    assert(orig_contact_matrix.shape == (N,N))

    new_group_names = []
    for name in orig_group_names:
        new_group_names.append("Vax: {}".format(name))
        new_group_names.append("Unvax: {}".format(name))

    new_group_params = []
    new_group_vax_statuses = []

    
    new_contact_matrix = np.zeros((2*N,2*N))
    group_idx = 0
    for group_params, vax_rate in zip(orig_group_params, vax_rates_by_group):
        assert(vax_rate >= 0 and vax_rate <= 1)
        orig_popsize = group_params['population_size']
        vax_popsize = int(orig_popsize * vax_rate)
        unvax_popsize = orig_popsize - vax_popsize
        assert(vax_popsize >= 0 and vax_popsize <= orig_popsize)
        assert(unvax_popsize >= 0 and unvax_popsize <= orig_popsize)
        if vax_popsize == 0 or unvax_popsize == 0:
            print("Warning: creating subgroup with population size 0, this might generate an error")

        g_vax = group_params.copy()
        g_unvax = group_params.copy()
        g_vax['population_size'] = vax_popsize
        if vax_age_sev_dist_files != None:
            vax_age_sev_dist_file = vax_age_sev_dist_files[group_idx]
            age_sev_params = load_age_sev_params(vax_age_sev_dist_file)
            g_vax['severity_prevalence'] = age_sev_params
        

        g_unvax['population_size'] = unvax_popsize

        new_group_params.append(g_vax)
        new_group_params.append(g_unvax)

        new_group_vax_statuses.append(1)
        new_group_vax_statuses.append(0)

        # fill in values new_contact_matrix[idx,:]
        vax_group_id = group_idx * 2
        unvax_group_id = vax_group_id + 1

        for j in range(N):
            orig_contact_freq = orig_contact_matrix[group_idx,j] # how many people a day someone from group j coughs on someone in group idx

            vax_contact_freq = orig_contact_freq * (vax_popsize / orig_popsize)
            unvax_contact_freq = orig_contact_freq * (unvax_popsize / orig_popsize)

            j_vax_group_id = j * 2
            j_unvax_group_id = j * 2 + 1

            #import pdb; pdb.set_trace()
            new_contact_matrix[vax_group_id, j_vax_group_id] = vax_contact_freq
            new_contact_matrix[vax_group_id, j_unvax_group_id] = vax_contact_freq

            new_contact_matrix[unvax_group_id, j_vax_group_id] = unvax_contact_freq
            new_contact_matrix[unvax_group_id, j_unvax_group_id] = unvax_contact_freq

        group_idx += 1

    return create_multigrp_vax_sim(new_group_params, new_group_vax_statuses, new_contact_matrix,
                                    vax_trans_mult, vax_susc_mult)

def create_multigrp_vax_sim(group_params, group_vax_statuses, interaction_matrix, 
                                vax_trans_mult, vax_susc_mult):

    vax_weighted_interaction_matrix = np.zeros(interaction_matrix.shape)

    N = len(group_params)
    for i in range(N):
        for j in range(N):
            j_to_i_contact_freq = interaction_matrix[i,j]

            if group_vax_statuses[j]:
                j_to_i_contact_freq *= vax_trans_mult

            if group_vax_statuses[i]:
                j_to_i_contact_freq *= vax_susc_mult

            vax_weighted_interaction_matrix[i,j] = j_to_i_contact_freq

    for i in range(N):
        group_params[i]['expected_contacts_per_day'] = vax_weighted_interaction_matrix[i,i]

    return MultiGroupSimulation(group_params, vax_weighted_interaction_matrix)

def load_vax_group_configs(vax_config_yaml, base_config_unvax=None):
    with open(vax_config_yaml, "rb") as f:
        vax_config = yaml.load(f)
    return process_vax_config(vax_config)

def process_vax_config(vax_config):
    base_config_yaml = vax_config['_base_config']
    base_config_unvax = load_params(base_config_yaml)[1]
    
    if base_config_unvax == None:
        base_config_yaml = vax_config['_base_config']
        base_config_unvax = load_params(base_config_yaml)[1]

    base_config_vax = base_config_unvax.copy()

    base_config_unvax['test_population_fraction'] = vax_config['unvax_test_freq']
    base_config_vax['test_population_fraction'] = vax_config['vax_test_freq']

    popsize = base_config_unvax['population_size']

    vax_popsize = int(popsize * vax_config['vax_proportion'])
    unvax_popsize = popsize - vax_popsize

    base_config_vax['population_size'] = vax_popsize
    base_config_unvax['population_size'] = unvax_popsize

    default_contacts = base_config_unvax['expected_contacts_per_day']

    # 0 corresponds to vax group, 1 corresponds to unvax
    # contact_matrix[i,j] means how many times an individual in j coughs on an individual in i
    # contact_matrix[i,0] = default_contacts * (vax_population / total_population)
    #
    contact_matrix = np.zeros((2,2))  
    contact_matrix[0,0] = default_contacts * (vax_popsize / popsize)
    contact_matrix[0,1] = default_contacts * (unvax_popsize / popsize)
    contact_matrix[1,0] = default_contacts * (vax_popsize / popsize)
    contact_matrix[1,1] = default_contacts * (unvax_popsize / popsize)

    base_config_vax['expected_contacts_per_day'] = contact_matrix[0,0]
    base_config_unvax['expected_contacts_per_day'] = contact_matrix[1,1]

    return [base_config_vax, base_config_unvax], [True, False], contact_matrix
