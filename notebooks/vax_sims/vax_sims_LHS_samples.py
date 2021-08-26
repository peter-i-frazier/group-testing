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
    'contacts_per_day_mult': (1.4,3.2),
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

    vax_rates = [0.9, 0.9, 0.9, 0.8]

    return [ug_ga_base_params, ug_other_base_params, grad_base_params, employee_base_params], \
            ["ug_ga", "ug_other", "grad", "employees"], \
            contact_matrix,\
            vax_rates

def map_lhs_point_to_vax_sim(lhs_point):
    base_params, base_group_names, contact_matrix, vax_rates = load_calibrated_params()


    vax_susc_mult = lhs_point[0]
    vax_trans_mult = lhs_point[1]

    contact_matrix = contact_matrix * lhs_point[2]

    for params in base_params:
        params['daily_outside_infection_p'] = params['daily_outside_infection_p'] * lhs_point[3]

    base_params[0]['cases_isolated_per_contact_trace'] = lhs_point[4]
    base_params[1]['cases_isolated_per_contact_trace'] = lhs_point[4]
    base_params[2]['cases_isolated_per_contact_trace'] = lhs_point[4]

    return generate_vax_unvax_multigroup_sim(base_params, base_group_names,
                                    vax_rates, contact_matrix,
                                    vax_trans_mult, vax_susc_mult)


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


def run_simulations(lhs_point, idx, output_folder):
    T=112
    n=50
    vax_sim = map_lhs_point_to_vax_sim(point)
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

   
    processes = []
    timestamp = time.time()
    output_folder = "./lhs_vax_sims:{}/".format(timestamp)
    os.mkdir(output_folder)
    for i in range(lhs_points.shape[0]):
        point = lhs_points[i,:]

        p = multiprocessing.Process(target = run_simulations, args = (point, i, output_folder))
        p.start()
        processes.append(p)
    print("done launching {} processes".format(len(processes)))
    for p in processes:
        p.join()

    print("processes finished running")
            









