import pandas as pd
import dill
import numpy as np

from statsmodels.api import OLS, add_constant
from util_functions import *

import os
import sys
module_path = os.path.abspath(os.path.join('../..'))
if module_path not in sys.path:
        sys.path.append(module_path + "/src/simulations_v2")
from load_params import load_params, update_sev_prevalence
from analysis_helpers import poisson_waiting_function

# list of parameter names that are varied in the uncertainty analysis
UNCERTAINTY_PARAMS_LIST = ['asymp_prob_mult', 'inital_prev_mult', 'R0', 'outside_inf_mult', 'daily_self_report_prob',
                       'ct_mult', 'ct_testing_ratio', 'test_sensitivity', 'test_noncompliance', 'E_time', 'ID_time',
                         'Sy_time', 'virtual_noncompliance', 'intermittent_non-compliance', 'virtual_r0_mult',
                                                                            'virtual_pop_size']

ADDITIONAL_VIRTUAL_PARAMS = ['virtual_noncompliance', 'intermittent_non-compliance', 'virtual_r0_mult', 'virtual_pop_size']


PARAM_BOUNDS = {
    'asymp_prob_mult': (24/47, 70/47), # Our nominal estimate for US population: 47%
    'inital_prev_mult': (0.5, 1.5),
    'R0': (1,4),
    'outside_inf_mult': (0.5, 1.5),
    'daily_self_report_prob': (0.22, 0.5),
    'ct_mult': (1,2),
    'ct_testing_ratio': (0.5, 1.5),
    'test_sensitivity': (0.4, 0.8),
    'test_noncompliance': (0.05, 0.15),
    'E_time': (1,3),
    'ID_time': (2,4),
    'Sy_time': (11,13),
    'virtual_noncompliance': (0.25, 0.75),
    'intermittent_non-compliance': (0.25,0.75),
    'virtual_r0_mult': (0.97, 1.5),
    'virtual_pop_size': (0,1), # Slider from min to max
}

def get_centre_point():
    centre = {}
    for param in PARAM_BOUNDS:
        lb, ub = PARAM_BOUNDS[param]
        centre[param] = (lb + ub) / 2
    return centre



def get_test_FNR(sensitivity, compliance):
    if 1 - (sensitivity * compliance) > 1:
        print(sensitivity, compliance)
    return 1 - (sensitivity * compliance)

def uncertainty_point_to_params_dict(uncertainty_point_dict):

    uncertainty_point = []
    for param in  UNCERTAINTY_PARAMS_LIST + ADDITIONAL_VIRTUAL_PARAMS:
        uncertainty_point.append(uncertainty_point_dict[param])
    
    res_params_list, res_interaction_matrix, res_group_names = get_nominal_params()
    virtual_persistent_noncompliance = uncertainty_point[12]
    virtual_ug_pop = 4500 * (1 - uncertainty_point[15]) + 7950 * uncertainty_point[15]
    virtual_gs_other_pop = 4770 * (1 - uncertainty_point[15]) + 5850 * uncertainty_point[15]
    virtual_params_list, virtual_interaction_matrix, virtual_group_names = get_virtual_params(virtual_persistent_noncompliance, virtual_ug_pop, virtual_gs_other_pop)

    # Asmptomatic Prob Mult
    for params in res_params_list:
        params['severity_prevalence'] = update_sev_prevalence(params['severity_prevalence'], uncertainty_point[0] * params['severity_prevalence'][0])
        if params['severity_prevalence'][0] > 1:
            params['severity_prevalence'] = [1,0,0,0]
    for params in virtual_params_list:
        params['severity_prevalence'] = update_sev_prevalence(params['severity_prevalence'], uncertainty_point[0] * params['severity_prevalence'][0])
        if params['severity_prevalence'][0] > 1:
            params['severity_prevalence'] = [1,0,0,0]
    
    # Initial Prevalence Mult
    for params in res_params_list:
        params['initial_ID_prevalence'] *= uncertainty_point[1]
    for params in virtual_params_list:
        params['initial_ID_prevalence'] *= uncertainty_point[1]
    
    # R0 adjustment
    res_interaction_matrix *= uncertainty_point[2]/2.5
    virtual_interaction_matrix *= uncertainty_point[2] * uncertainty_point[14]/2.5
    
    # Outside inf mult
    for params in res_params_list:
        params['daily_outside_infection_p'] *= uncertainty_point[3]
    for params in virtual_params_list:
        params['daily_outside_infection_p'] *= uncertainty_point[3]
    
    # Daily self-report prob
    for params in res_params_list:
        params['severe_symptoms_daily_self_report_p'] = uncertainty_point[4]
    for params in virtual_params_list:
        params['severe_symptoms_daily_self_report_p'] = uncertainty_point[4]
        
    # CT mult
    for params in res_params_list:
        params['cases_isolated_per_contact'] *= uncertainty_point[5]
    for params in virtual_params_list:
        params['cases_isolated_per_contact'] *= uncertainty_point[5]
    
    # CT testing ratio
    for params in res_params_list:
        params['contact_trace_testing_frac'] = uncertainty_point[6]
    for params in virtual_params_list:
        params['contact_trace_testing_frac'] = uncertainty_point[6]
    
    # Test sensitivity and Test compliance (note: non-compliance is provided in uncertainty point)
    for params in res_params_list:
        params['test_protocol_QFNR'] = get_test_FNR(uncertainty_point[7], 1-uncertainty_point[8])
    for params in virtual_params_list:
        params['test_protocol_QFNR'] = get_test_FNR(uncertainty_point[7], 1-(uncertainty_point[13]))

    # E_time, ID_time, Sy_time
    for params in res_params_list:
        params['exposed_time_function'] = poisson_waiting_function(7, uncertainty_point[9])
        params['ID_time_function'] = poisson_waiting_function(8, uncertainty_point[10])
        params['SyID_mild_time_function'] = poisson_waiting_function(20, uncertainty_point[11])
        params['SyID_severe_time_function'] = poisson_waiting_function(20, uncertainty_point[11])
    
    for params in virtual_params_list:
        params['exposed_time_function'] = poisson_waiting_function(7, uncertainty_point[9])
        params['ID_time_function'] = poisson_waiting_function(8, uncertainty_point[10])
        params['SyID_mild_time_function'] = poisson_waiting_function(20, uncertainty_point[11])
        params['SyID_severe_time_function'] = poisson_waiting_function(20, uncertainty_point[11])
        
    return (res_params_list, res_interaction_matrix, res_group_names), \
            (virtual_params_list, virtual_interaction_matrix, virtual_group_names)

def calculate_pessimistic_scenario(results):
    # the keys in dict(results.params) specify whether this is for residential
    # or virtual vs. residential
    lr_results = dict(results.params)
    range_dict = dict()
    params = set(lr_results.keys()) - set(['const'])
    for param in params:
        range_dict[param] = (PARAM_BOUNDS[param][1] + PARAM_BOUNDS[param][0])/2

    
    
    sum_squares = 0
    for param in params:
        sum_squares += ((lr_results[param]*range_dict[param])/2) ** 2

    # calculate pessimistic scenario based on available params
    pess_scenario = dict()
    for param in params:
        pess_scenario[param] = np.mean(PARAM_BOUNDS[param]) + \
            ((lr_results[param] * (range_dict[param])**2) / 2) / np.sqrt(sum_squares)

    
    # add default virtual params if not present
    default_virtual_param_vals = {param:(PARAM_BOUNDS[param][1] + PARAM_BOUNDS[param][0])/2 for param in ADDITIONAL_VIRTUAL_PARAMS}
    for virtual_param, val in default_virtual_param_vals.items():
        if virtual_param not in params:
            pess_scenario[virtual_param] = val

    return pess_scenario

def residential_regression(scenario_data):
    residential_columns = scenario_data.columns[0:12]
    residential_target = 'res_cornell_inf_50'
    X_res = scenario_data[residential_columns]
    Y_res_outcomes = np.array(scenario_data[[residential_target]])

    X = add_constant(X_res)
    model = OLS(Y_res_outcomes,X)
    results = model.fit()
    return results

def virtual_vs_residential_regression(scenario_data):
    residential_columns = scenario_data.columns[0:12]
    residential_target = 'res_cornell_inf_50'
    Y_res_outcomes = np.array(scenario_data[[residential_target]])

    virtual_columns = scenario_data.columns[0:16]
    virtual_infs = 'vir_cornell_inf_50'
    Y_vir_outcomes = np.array(scenario_data[[virtual_infs]])

    X = scenario_data[virtual_columns]

    Y_target = Y_res_outcomes - Y_vir_outcomes

    X = add_constant(X)
    model = OLS(Y_target,X)
    results = model.fit()
    return results


def params_dict_to_uncertainty_point():
    # will have to make assumptions about the parameters in the dictionary that don't get
    # varied in the LHS experiments
    pass



def get_stats(inf_matrix):
    cornell_inf = np.array(inf_matrix)[:,:-1].sum(axis=1)
    ithaca_inf = np.array(inf_matrix)[:,-1]
    return np.quantile(cornell_inf, [0.1,0.5,0.9]), np.quantile(ithaca_inf, [0.1,0.5,0.9])


def load_sim_output(sim_output_files):
    scenario_data = pd.DataFrame(columns=UNCERTAINTY_PARAMS_LIST+\
            ['res_cornell_inf_10','res_cornell_inf_50','res_cornell_inf_90','res_ithaca_inf_10',
                'res_ithaca_inf_50','res_ithaca_inf_90']+\
            ['vir_cornell_inf_10','vir_cornell_inf_50','vir_cornell_inf_90',
                'vir_ithaca_inf_10','vir_ithaca_inf_50','vir_ithaca_inf_90'])

    for fname in sim_output_files:
        with open(fname, 'rb') as fhandle:
            [uncertainty_point, res_inf_matrix, res_hosp_matrix, virtual_inf_matrix, virtual_hosp_matrix] = dill.load(fhandle)

        new_row = dict()
        for index, col_name in enumerate(UNCERTAINTY_PARAMS_LIST):
            if type(uncertainty_point) == dict:
                new_row[col_name] = uncertainty_point[col_name]
            else:
                new_row[col_name] = uncertainty_point[index]

        res_cornell_inf_quantiles, res_ithaca_inf_quantiles = get_stats(res_inf_matrix)
        new_row['res_cornell_inf_10'] = res_cornell_inf_quantiles[0]
        new_row['res_cornell_inf_50'] = res_cornell_inf_quantiles[1]
        new_row['res_cornell_inf_90'] = res_cornell_inf_quantiles[2]
        new_row['res_ithaca_inf_10'] = res_ithaca_inf_quantiles[0]
        new_row['res_ithaca_inf_50'] = res_ithaca_inf_quantiles[1]
        new_row['res_ithaca_inf_90'] = res_ithaca_inf_quantiles[2]

        if virtual_inf_matrix != None:
            vir_cornell_inf_quantiles, vir_ithaca_inf_quantiles = get_stats(virtual_inf_matrix)
            new_row['vir_cornell_inf_10'] = vir_cornell_inf_quantiles[0]
            new_row['vir_cornell_inf_50'] = vir_cornell_inf_quantiles[1]
            new_row['vir_cornell_inf_90'] = vir_cornell_inf_quantiles[2]
            new_row['vir_ithaca_inf_10'] = vir_ithaca_inf_quantiles[0]
            new_row['vir_ithaca_inf_50'] = vir_ithaca_inf_quantiles[1]
            new_row['vir_ithaca_inf_90'] = vir_ithaca_inf_quantiles[2]
        else:
            new_row['vir_cornell_inf_10'] = None
            new_row['vir_cornell_inf_50'] = None
            new_row['vir_cornell_inf_90'] = None
            new_row['vir_ithaca_inf_10'] = None
            new_row['vir_ithaca_inf_50'] = None
            new_row['vir_ithaca_inf_90'] = None

            

        scenario_data = scenario_data.append(new_row, ignore_index=True)
    return scenario_data 


