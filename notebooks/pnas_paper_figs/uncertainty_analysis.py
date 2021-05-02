import pandas as pd
import dill
import numpy as np

from statsmodels.api import OLS, add_constant

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

def calculate_pessimistic_scenario(results):
    # the keys in dict(results.params) specify whether this is for residential
    # or virtual vs. residential
    lr_results = dict(results.params)
    range_dict = dict()
    params = set(lr_results.keys()) - set(['const'])
    for param in params:
        range_dict[param] = (PARAM_BOUNDS[param][1] - PARAM_BOUNDS[param][0])/2

    sum_squares = 0
    for param in params:
        sum_squares += ((lr_results[param]*range_dict[param])/2) ** 2

    pess_scenario = dict()
    for param in params:
        pess_scenario[param] = np.mean(PARAM_BOUNDS[param]) + \
            ((lr_results[param] * (range_dict[param])**2) / 2) / np.sqrt(sum_squares)
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

def uncertainty_point_to_params_dict():
    pass

def params_dict_to_uncertainty_point():
    # will have to make assumptions about the parameters in the dictionary that don't get
    # varied in the LHS experiments
    pass



def get_stats(inf_matrix):
    cornell_inf = np.array(inf_matrix)[:,:-1].sum(axis=1)
    ithaca_inf = np.array(inf_matrix)[:,-1]
    return np.quantile(cornell_inf, [0.1,0.5,0.9]), np.quantile(ithaca_inf, [0.1,0.5,0.9])


def load_lhs_output(sim_output_files):
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
            new_row[col_name] = uncertainty_point[index]

        res_cornell_inf_quantiles, res_ithaca_inf_quantiles = get_stats(res_inf_matrix)
        new_row['res_cornell_inf_10'] = res_cornell_inf_quantiles[0]
        new_row['res_cornell_inf_50'] = res_cornell_inf_quantiles[1]
        new_row['res_cornell_inf_90'] = res_cornell_inf_quantiles[2]
        new_row['res_ithaca_inf_10'] = res_ithaca_inf_quantiles[0]
        new_row['res_ithaca_inf_50'] = res_ithaca_inf_quantiles[1]
        new_row['res_ithaca_inf_90'] = res_ithaca_inf_quantiles[2]

        vir_cornell_inf_quantiles, vir_ithaca_inf_quantiles = get_stats(virtual_inf_matrix)
        new_row['vir_cornell_inf_10'] = vir_cornell_inf_quantiles[0]
        new_row['vir_cornell_inf_50'] = vir_cornell_inf_quantiles[1]
        new_row['vir_cornell_inf_90'] = vir_cornell_inf_quantiles[2]
        new_row['vir_ithaca_inf_10'] = vir_ithaca_inf_quantiles[0]
        new_row['vir_ithaca_inf_50'] = vir_ithaca_inf_quantiles[1]
        new_row['vir_ithaca_inf_90'] = vir_ithaca_inf_quantiles[2]
        

        scenario_data = scenario_data.append(new_row, ignore_index=True)
    return scenario_data 


