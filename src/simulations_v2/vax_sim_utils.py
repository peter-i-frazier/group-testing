from multi_group_simulation import MultiGroupSimulation
import numpy as np
import yaml
from load_params import load_params

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

    return MultiGroupSimulation(group_params, vax_weighted_interaction_matrix)

def load_vax_group_configs(vax_config_yaml):
    with open(vax_config_yaml, "rb") as f:
        vax_config = yaml.load(f)

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
