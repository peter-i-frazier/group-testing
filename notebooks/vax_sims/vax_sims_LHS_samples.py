import sys
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
    'vax_susc_mult': (0.1, 0.5),
    'vax_transmission_mult': (0.4, 0.8),
    'contacts_per_day_mult': (1.4,3.2),
    'outside_infection_rate_mult': (1, 3),
    'cases_isolated_per_contact_trace': (0.5,1.5)
}

lb = list()
ub = list()

for param in UNCERTAINTY_PARAMS:
    lb.append(UNCERTAINTY_PARAM_RANGES[param][0])
    ub.append(UNCERTAINTY_PARAM_RANGES[param][1])


np.random.seed(2021)

dim = len(UNCERTAINTY_PARAMS)
num_samples = 2000
lhs_points = pyDOE.lhs(dim, samples=num_samples)

for i in range(dim):
    lhs_points[:, i] = (1 - lhs_points[:,i]) * lb[i] + lhs_points[:,i] * ub[i]

def load_calibrated_params():
    employee_base_params = load_params("./vax_sim_nominal_params/employee_nominal.yaml")[1]
    grad_base_params = load_params("./vax_sim_nominal_params/grad_research_nominal.yaml")[1]
    ug_ga_base_params = load_params("./vax_sim_nominal_params/ug_greek_athlete_nominal.yaml")[1]
    ug_other_base_params = load_params("./vax_sim_nominal_params/ug_other_nominal.yaml")[1]
    
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



for i in range(lhs_points.shape[0]):
# for i in range(20):
    point = lhs_points[i,:]
    vax_sim = map_lhs_point_to_vax_sim(point)

    vax_sim.run_new_trajectory(25)
    import pdb; pdb.set_trace()

            









