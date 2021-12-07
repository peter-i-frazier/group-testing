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

def sample_from_prior():
    return_point = [0] * 6
    while min(return_point) <= 0:
        idx = 0
        for param in UNCERTAINTY_PARAMS:
            mean = np.mean(UNCERTAINTY_PARAM_RANGES[param])
            sd = (UNCERTAINTY_PARAM_RANGES[param][1] - UNCERTAINTY_PARAM_RANGES[param][0])/(2*1.96)
            return_point[idx] = np.random.normal(mean, sd)
            idx += 1
    return return_point

def map_lhs_point_to_vax_sim(lhs_point, param_modifiers=None, vax_rates=None, omicron_multiplier=1):
    base_params, base_group_names, contact_matrix, default_vax_rates = load_calibrated_params()
    if vax_rates == None:
        vax_rates = default_vax_rates

    vax_susc_mult = lhs_point[0]
    vax_trans_mult = lhs_point[1]

    contact_matrix = contact_matrix * lhs_point[2] * omicron_multiplier

    for params in base_params:
        params['daily_outside_infection_p'] = params['daily_outside_infection_p'] * lhs_point[3]
        params['contact_tracing_delay'] = 3


    base_params[0]['cases_isolated_per_contact'] *= lhs_point[4]
    base_params[1]['cases_isolated_per_contact'] *= lhs_point[4]
    base_params[2]['cases_isolated_per_contact'] *= lhs_point[4]
    
    base_params[0]['initial_ID_prevalence'] = lhs_point[5]
    base_params[1]['initial_ID_prevalence'] = lhs_point[5]
    base_params[2]['initial_ID_prevalence'] = lhs_point[5]
    
    vax_sim = generate_vax_unvax_multigroup_sim(base_params, base_group_names,
                                    vax_rates, contact_matrix,
                                    vax_trans_mult, vax_susc_mult)

    update_vax_sim_params(vax_sim, param_modifiers)

    return vax_sim


PARAMS_POST_MOVEIN = {'ug_ga_vax_test_frequency': 2/7, 'ug_ga_unvax_test_frequency': 2/7,
            'ug_other_vax_test_frequency': 1/7, 'ug_other_unvax_test_frequency': 2/7,
            'grad_vax_test_frequency': 1/7, 'grad_unvax_test_frequency': 2/7,
            'employee_vax_test_frequency': 1/7, 'employee_unvax_test_frequency': 2/7,
            'test_delay': 1, 'max_time_pre_ID': 2}

PARAMS_PRE_MOVEIN = {'ug_ga_vax_test_frequency': 1/7, 'ug_ga_unvax_test_frequency': 2/7,
            'ug_other_vax_test_frequency': 1/7, 'ug_other_unvax_test_frequency': 2/7,
            'grad_vax_test_frequency': 1/7, 'grad_unvax_test_frequency': 2/7,
            'employee_vax_test_frequency': 1/7, 'employee_unvax_test_frequency': 2/7,
            'test_delay': 2, 'max_time_pre_ID': 2}


def get_cum_infections(df):
    return df[['cumulative_mild', 'cumulative_severe']].iloc[df.shape[0]-1].sum()


def get_cum_inf_trajectory(df):
    return np.sum(df[['cumulative_mild', 'cumulative_severe']], axis=1)


def run_new_trajectory(sim, T, change_t, override_premovein_params=None):
    sim.reset_initial_state()
    assert(override_premovein_params != None)
    if override_premovein_params != None:
        update_vax_sim_params(sim, override_premovein_params)
    else:
        update_vax_sim_params(sim, PARAMS_PRE_MOVEIN)
        
    for t in range(T):
        sim.step()
        #if t == change_t:
        #    update_vax_sim_params(sim, PARAMS_POST_MOVEIN)

    for single_group_sim in sim.sims:
        single_group_sim.update_severity_levels()

    sim_df = sim.sims[0].sim_df
    for sim in sim.sims[1:]:
        sim_df = sim_df.add(sim.sim_df)
    return sim_df

CHANGE_T = 6 # start of simulation is aug 23, assume changes occur on Aug 29
def run_multigroup_sim(sim, T, override_premovein_params=None):
    run_new_trajectory(sim, T, CHANGE_T, override_premovein_params)
    inf_trajs_by_group = []
    for group in sim.sims:
        df = group.sim_df
        inf_trajs_by_group.append(get_cum_inf_trajectory(df))
    return inf_trajs_by_group


def get_centre_point():
    centre = {}
    for param in UNCERTAINTY_PARAM_RANGES:
        lb, ub = UNCERTAINTY_PARAM_RANGES[param]
        centre[param] = (lb + ub) / 2
    return centre


def run_multiple_trajs(sim, T, n, override_premovein_params=None):
    infs_by_group_list = []
    for _ in range(n):
        infs_by_group = run_multigroup_sim(sim,T, override_premovein_params)
        infs_by_group_list.append(infs_by_group)
    return infs_by_group_list

def get_timestamp():
    return str(time.time()).split('.')[0]


def sample_and_save(vax_rates, 
                    omicron_multiplier,
                    save_folder, nsamples=100, T=112):
    gc.collect()


    test_policy_premovein = PARAMS_PRE_MOVEIN
    test_policy_postmovein = PARAMS_POST_MOVEIN

    test_policy_premovein['test_delay'] = 1
    test_policy_premovein['max_time_pre_ID'] = 2

    test_policy_postmovein['test_delay'] = 1
    test_policy_postmovein['max_time_pre_ID'] = 2


    df = load_posterior_df()
    posterior_sample = np.random.multinomial(nsamples, list(df['posterior']))
    posterior_point_idxs = []
    for idx, count in enumerate(posterior_sample):
        for _ in range(count):
            posterior_point_idxs.append(idx)

    assert(len(posterior_point_idxs) == nsamples)

    posterior_points = []
    inf_trajs_by_group = []

    for idx in posterior_point_idxs:

        posterior_point = list(df[UNCERTAINTY_PARAMS].iloc[idx])

        posterior_sim = map_lhs_point_to_vax_sim(posterior_point, test_policy_premovein, vax_rates, omicron_multiplier=omicron_multiplier)
        infs_by_group = run_multigroup_sim(posterior_sim, T, override_premovein_params = test_policy_premovein)

        posterior_points.append(posterior_point)
        inf_trajs_by_group.append(infs_by_group)


    dill_path = save_folder + 'test_policy_{}_vax_rates_{}.dill'.format(test_policy_idx, vax_rates_idx)

    pickle.dump([posterior_point_idxs, posterior_point, inf_trajs_by_group], open(dill_path, 'wb'))



import multiprocessing as mp
import gc
from joblib import Parallel, delayed
import multiprocessing
import pandas as pd


def load_posterior_df():
    df = pd.read_csv('posterior_csvs/21_10_08_14:44_posteriors.csv')
    return df





vax_rates_to_try = [1,1,1,1]
omicron_multipliers = np.linspace(1,3.5,5)

param_modifiers = PARAMS_PRE_MOVEIN
override_params = PARAMS_POST_MOVEIN

if __name__ == "__main__":

    param_modifiers = PARAMS_PRE_MOVEIN
    override_params = PARAMS_POST_MOVEIN
    nsamples = 1000
    T = 112 
    
    save_folder = '../../notebooks/vax_sims/posterior_test_frequency_sims_{}/'.format(get_timestamp())
    os.mkdir(save_folder)



    print("kicking off processes now, saving results in {}".format(save_folder))
    num_cores = multiprocessing.cpu_count()
    results = Parallel(n_jobs=num_cores)(delayed(sample_and_save)(vax_rates_to_try, 
                                                                    omicron_multiplier, 
                                                                    save_folder, 
                                                                    nsamples=nsamples, 
                                                                    T=T) for omicron_multiplier in omicron_multipliers)
