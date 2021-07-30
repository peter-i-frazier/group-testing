from multi_group_simulation import MultiGroupSimulation
import numpy as np

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
