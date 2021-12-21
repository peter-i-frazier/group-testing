import sys
import os
import numpy as np
import multiprocessing
import pickle
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import numpy as np
from statsmodels.api import OLS, add_constant
import time
from scipy.stats import norm
import datetime as dt

module_path = os.path.abspath(os.path.join('.'))
sys.path.append(module_path + "/../../notebooks/vax_sims/")
from vax_sims_LHS_samples import *
from plot_utils import *
sys.path.append(module_path + "/../../src/simulations_v2")
from stochastic_simulation import StochasticSimulation


UNCERTAINTY_PARAMS = ['vax_susc_mult', 'vax_transmission_mult', 'contacts_per_day_mult', 'outside_infection_rate_mult',
                              'cases_isolated_per_contact_trace', 'initial_ID_prevalence']

UNCERTAINTY_PARAM_RANGES = {
            'vax_susc_mult': (0.097608, 0.941192), # 0.5194 +/- 1.96 * 0.2152
                'vax_transmission_mult': (0.25, 1),
                    'contacts_per_day_mult': (0.9,2.7),
                        'outside_infection_rate_mult': (1, 5),
                            'cases_isolated_per_contact_trace': (0.5,1.5),
                                'initial_ID_prevalence': (0.003, 0.0054)
                                }

def load_posterior_df():
    df = pd.read_csv('../vax_sims/posterior_csvs/21_10_08_14:44_posteriors.csv')
    return df

def load_virtual_params():
    df = load_posterior_df()
    MLE_point = list(df.sort_values('posterior', ascending=False).head(1)[UNCERTAINTY_PARAMS].iloc[0])

    group_params = load_params("../vax_sims/vax_sim_nominal_params/ug_greek_athlete_nominal.yaml")[1]

    vax_susc_mult = MLE_point[0]
    vax_trans_mult = MLE_point[1]

    # Contacts per day should be 7.59 (7.36 + 0.23) from the contact matrix
    omicron_multiplier = 5
    social_distance_multiplier = 1
    group_params['expected_contacts_per_day'] = 7.59 * vax_susc_mult * vax_trans_mult * MLE_point[2] * omicron_multiplier * social_distance_multiplier

    group_params['daily_outside_infection_p'] = group_params['daily_outside_infection_p'] * MLE_point[3] * 10
    group_params['cases_isolated_per_contact'] *= MLE_point[4]
    group_params['initial_ID_prevalence'] = 0.04

    group_params['population_size'] = 8018
    group_params['test_population_fraction'] = 1/7

    return group_params

def update_virtual_params(update_dict, base_params=None):
    if base_params is None:
        base_params = load_virtual_params()

    if 'population_size' in update_dict.keys():
        base_params['population_size'] = update_dict['population_size']

    if 'test_population_fraction' in update_dict.keys():
        base_params['test_population_fraction'] = update_dict['test_population_fraction']

    if 'cases_isolated_per_contact' in update_dict.keys():
        base_params['cases_isolated_per_contact'] = update_dict['cases_isolated_per_contact']

    if ('omicron_multiplier' in update_dict.keys()) and ('social_distance_multiplier' in update_dict.keys()):
        base_params['expected_contacts_per_day'] *= update_dict['omicron_multiplier'] * update_dict['social_distance_multiplier'] / 5

    if 'initial_R_fraction' in update_dict.keys():
        base_params['initial_R_count'] = int(np.round(update_dict['initial_R_fraction'] * base_params['population_size']))

    return base_params

def run_simulation(n=100, T=35):
    sim_results = list()
    for _ in range(n):
        sim_results.append(sim.run_new_trajectory(T))
    return sim_results

if __name__ == '__main__':
    print(load_virtual_params())

