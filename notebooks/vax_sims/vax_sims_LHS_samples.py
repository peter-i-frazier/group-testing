import sys
import time
import yaml
import os
import numpy as np
import multiprocessing
import dill
import matplotlib.pyplot as plt
import pandas as pd
import pyDOE
from multiprocessing import Process
from scipy.stats import norm

module_path = os.path.abspath(os.path.join('../..'))
if module_path not in sys.path:
    sys.path.append(module_path + "/src/simulations_v2")
    sys.path.append(module_path + "/notebooks/pnas_paper_figs")

from util_functions import get_cum_infections

from load_params import load_params, update_sev_prevalence
from analysis_helpers import poisson_waiting_function

from multi_group_simulation import MultiGroupSimulation

from vax_sim_utils import generate_vax_unvax_multigroup_sim

UNCERTAINTY_PARAMS = ['vax_susc_mult', 'vax_transmission_mult', 'contacts_per_day_mult', 'outside_infection_rate_mult', 'cases_isolated_per_contact_trace']

UNCERTAINTY_PARAM_RANGES = {
    'vax_susc_mult': (0.097608, 0.941192), # 0.5194 +/- 1.96 * 0.2152
    'vax_transmission_mult': (0.25, 1),
    'contacts_per_day_mult': (0.9,2.7),
    'outside_infection_rate_mult': (1, 5),
    'cases_isolated_per_contact_trace': (0.5,1.5)
}

def load_calibrated_params():
    employee_base_params = load_params("./vax_sim_nominal_params/employee_nominal.yaml")[1]
    grad_base_params = load_params("./vax_sim_nominal_params/grad_research_nominal.yaml")[1]
    ug_ga_base_params = load_params("./vax_sim_nominal_params/ug_greek_athlete_nominal.yaml")[1]
    ug_other_base_params = load_params("./vax_sim_nominal_params/ug_other_nominal.yaml")[1]
    
    ug_ga_base_params['initial_ID_prevalence'] = 0.003
    ug_other_base_params['initial_ID_prevalence'] = 0.003
    
    # order is ug_ga, ug_other, grad, employees
    contact_matrix = np.matrix(
            [[0.736, 0.023, 0, 0],
                [.028, .148, 0, 0], 
                [0, 0.023, .067, 0],
                [0,0,0,1]]) * 10

    vax_rates = [0.95, 0.95, 0.95, 0.8]

    return [ug_ga_base_params, ug_other_base_params, grad_base_params, employee_base_params], \
            ["ug_ga", "ug_other", "grad", "employees"], \
            contact_matrix,\
            vax_rates


def nominal_params_vax_sim(param_modifiers={}):
    nominal_point = np.array([np.mean(UNCERTAINTY_PARAM_RANGES[key]) for key in UNCERTAINTY_PARAMS])
    if 'contacts_per_day_mult' in param_modifiers:
        nominal_point[2] = param_modifiers['contacts_per_day_mult']
    return map_lhs_point_to_vax_sim(nominal_point)

def map_lhs_point_to_vax_sim(lhs_point, param_modifiers=None):
    base_params, base_group_names, contact_matrix, vax_rates = load_calibrated_params()


    vax_susc_mult = lhs_point[0]
    vax_trans_mult = lhs_point[1]

    contact_matrix = contact_matrix * lhs_point[2]

    for params in base_params:
        params['daily_outside_infection_p'] = params['daily_outside_infection_p'] * lhs_point[3]
        

    base_params[0]['cases_isolated_per_contact'] *= lhs_point[4]
    base_params[1]['cases_isolated_per_contact'] *= lhs_point[4]
    base_params[2]['cases_isolated_per_contact'] *= lhs_point[4]

    vax_sim = generate_vax_unvax_multigroup_sim(base_params, base_group_names,
                                    vax_rates, contact_matrix,
                                    vax_trans_mult, vax_susc_mult)

    update_vax_sim_params(vax_sim, param_modifiers)

    return vax_sim


def update_vax_sim_params(vax_sim, param_modifiers):
    # order is [ug_ga_vax, ug_ga_unvax, ug_other_vax, ug_other_unvax, grad_vax, grad_unvax, employee_vax, employee_unvax]
    # which corresponds to [0,       1,            2,              3,        4,          5,            6,              7]

    test_frequency_specifiers = ['ug_ga_vax_test_frequency', 'ug_ga_unvax_test_frequency', 
            'ug_other_vax_test_frequency', 'ug_other_unvax_test_frequency',
            'grad_vax_test_frequency', 'grad_unvax_test_frequency',
            'employee_vax_test_frequency', 'employee_unvax_test_frequency']
                            
    
    for idx, frequency_specifier in enumerate(test_frequency_specifiers):
        if frequency_specifier in param_modifiers:
            vax_sim.sims[idx].test_pop_fraction = param_modifiers[frequency_specifier]


    if 'contact_tracing_delay' in param_modifiers:
        for idx in range(len(vax_sim.sims)):
            vax_sim.sims[idx].contact_tracing_delay = param_modifiers['contact_tracing_delay']

    # Changing initial prevalence
    if 'initial_ID_prevalence' in param_modifiers:
        assert len(param_modifiers['initial_ID_prevalence']) == len(vax_sim.sims)
        for idx, init_prevalence in enumerate(param_modifiers['initial_ID_prevalence']):
            vax_sim.sims[idx].init_ID_prevalence = init_prevalence

    # Changing contact tracing effectiveness
    if 'cases_isolated_per_contact_mult' in param_modifiers:
        for sim in vax_sim.sims:
            sim.cases_isolated_per_contact *= param_modifiers['cases_isolated_per_contact_mult']

    # Test delay of d days implemented by adding infectious period of d days before ID state
    if 'test_delay' in param_modifiers and param_modifiers['test_delay'] > 0:
        for sim in vax_sim.sims:
            sim.pre_ID_state = 'infectious'
            if 'max_time_pre_ID' not in param_modifiers:
                sim.max_time_pre_ID = param_modifiers['test_delay']
            else:
                sim.max_time_pre_ID = param_modifiers['max_time_pre_ID']
            sim.sample_pre_ID_times = constant_delay_function(param_modifiers['test_delay'], sim.max_time_pre_ID) 

def constant_delay_function(time, max_time):
    array = np.zeros(max_time+1)
    array[time] = 1
    return (lambda n: np.random.multinomial(n, array))

def get_cum_infections(df):
    return df[['cumulative_mild', 'cumulative_severe']].iloc[df.shape[0]-1].sum()


def run_multigroup_sim(sim, T):
    sim.run_new_trajectory(T)
    infs_by_group = []
    for group in sim.sims:
        df = group.sim_df
        infs_by_group.append(get_cum_infections(df))
    return infs_by_group


def run_multiple_trajectories(sim, T, n):
    infs_by_group_list = []
    for _ in range(n):
        infs_by_group = run_multigroup_sim(sim,T)
        infs_by_group_list.append(infs_by_group)
    return infs_by_group_list


def run_simulations(lhs_point, idx, output_folder, param_modifiers = None):
    T=112
    n=50
    vax_sim = map_lhs_point_to_vax_sim(point, param_modifiers)
    list_of_infs_by_group = run_multiple_trajectories(vax_sim, T, n)
    with open(output_folder + "lhs_point_{}.dill".format(idx), "wb") as f:
        dill.dump(lhs_point, f)
    with open(output_folder + "list_of_infs_by_group_{}.dill".format(idx), "wb") as f:
        dill.dump(list_of_infs_by_group, f)

if __name__ == "__main__":


    lb = list()
    ub = list()

    for param in UNCERTAINTY_PARAMS:
        lb.append(UNCERTAINTY_PARAM_RANGES[param][0])
        ub.append(UNCERTAINTY_PARAM_RANGES[param][1])


    np.random.seed(2021)

    dim = len(UNCERTAINTY_PARAMS)
    num_samples = 200
    lhs_points = pyDOE.lhs(dim, samples=num_samples)

    for i in range(dim):
        lhs_points[:, i] = (1 - lhs_points[:,i]) * lb[i] + lhs_points[:,i] * ub[i]

    param_modifiers = {'ug_ga_vax_test_frequency': 1/7, 'ug_ga_unvax_test_frequency': 2/7,
            'ug_other_vax_test_frequency': 1/7, 'ug_other_unvax_test_frequency': 2/7,
            'grad_vax_test_frequency': 1/7, 'grad_unvax_test_frequency': 2/7,
            'employee_vax_test_frequency': 1/7, 'employee_unvax_test_frequency': 2/7}
   
    processes = []
    timestamp = time.time()
    output_folder = "./lhs_vax_sims:{}/".format(timestamp)
    os.mkdir(output_folder)
    for i in range(lhs_points.shape[0]):
        point = lhs_points[i,:]

        p = multiprocessing.Process(target = run_simulations, args = (point, i, output_folder, param_modifiers))
        p.start()
        processes.append(p)
    print("done launching {} processes".format(len(processes)))
    for p in processes:
        p.join()

    print("processes finished running")
            









