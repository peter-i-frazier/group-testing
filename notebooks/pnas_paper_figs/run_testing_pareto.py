import sys
import os
import numpy as np
import multiprocessing
import dill
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import norm
import itertools as it

from util_functions import *
from uncertainty_analysis import *
from sim_helper_functions import *

def get_tests_per_day(test_policy, params_list):
    assert len(test_policy) == len(params_list)
    tests_per_day = 0
    for index, params in enumerate(params_list[:6]):
        tests_per_day += np.ceil(params['population_size'] * test_policy[index])
    return tests_per_day

def equal_freq_policy(tests_per_day, params_list):
    on_campus_pop = 0
    for index in range(6):
        on_campus_pop += params_list[index]['population_size']
    equal_policy = list()
    for _ in range(6):
        equal_policy.append(tests_per_day / on_campus_pop)
    equal_policy.append(1/30)
    equal_policy.append(0)
    return equal_policy

def run_simulation(uncertainty_point,
                    test_policy,
                    filename,
                    point_id=None,
                    nreps=50, T=112, base_seed=100000):
    np.random.seed(base_seed)
    # get params
    (res_params_list, res_interaction_matrix, res_group_names),\
            (virtual_params_list, virtual_interaction_matrix, virtual_group_names) \
            = uncertainty_point_to_params_dict(uncertainty_point)

    print("running sim with id {}, nreps {}, T {}, filename {}".format(point_id, nreps, T, filename))
    # run simulations
    # Residential Simulation

    # Running res sims
    print('running residential sim with id {}'.format(point_id))
    res_tests_per_day, res_inf_matrix, res_hosp_matrix = evaluate_testing_policy(res_params_list, res_interaction_matrix, res_group_names, test_policy, T, nreps)

    virtual_tests_per_day, virtual_inf_matrix, virtual_hosp_matrix = (None, None, None)

    # save output
    file = open(filename, mode='wb')
    dill.dump([uncertainty_point, res_tests_per_day, res_inf_matrix, res_hosp_matrix], file)
    file.close()
    return


def run_new_process(uncertainty_point, test_policy, filename, point_id, nreps=50, T=112, base_seed = 1000000):
    p = Process(target = run_simulation, args = (uncertainty_point, test_policy, filename, point_id, nreps, T, base_seed))
    p.start()
    return p


def run_sims_new_process(uncertainty_point_dicts, test_policies, output_fnames, nreps=50, T=112, wait_for_processes_to_join=True):
    seed = np.random.randint(1000000, 10000000)
    idx = 0
    processes = []
    for uncertainty_dict, output_fname, test_policy in zip(uncertainty_point_dicts, output_fnames, test_policies):
        idx += 1
        p = run_new_process(uncertainty_dict, test_policy, output_fname, idx, nreps=nreps, T=T, base_seed=seed+idx)
        processes.append(p)

    print("launched {} processes".format(len(processes)))
    if wait_for_processes_to_join:
        for p in processes:
            p.join()
        print("done running processes")
        return processes
    else:
        return processes



if __name__ == "__main__":
    np.random.seed(2021)

    lhs_output_sim_files = []
    for i in range(2000):
        fname = '/home/aaj54/group-testing/notebooks/apr_29_scenarios/point_{}.dill'.format(i)
        lhs_output_sim_files.append(fname)

    scenario_data = load_sim_output(lhs_output_sim_files)
    res_results = residential_regression(scenario_data)
    res_pessimistic = calculate_pessimistic_scenario(res_results)

    groups_tested_2x_week = list()
    for ngroups in range(7):
        for entry in it.combinations(range(6), ngroups):
            groups_tested_2x_week.append(entry)

    test_policies = list()
    for groups_tested in groups_tested_2x_week:
        baseline_policy = [1/7,1/7,1/7,1/7,1/7,1/7,1/30,0]
        for group in groups_tested:
            baseline_policy[group] = 2/7
        test_policies.append(baseline_policy)

#    test_policies = [[1/7,1/7,1/7,1/7,1/7,1/7,1/14,0],
#    [2/7,1/7,1/7,1/7,1/7,1/7,1/30,0],
#    [2/7,2/7,1/7,1/7,1/7,1/7,1/30,0],
#    [2/7,2/7,1/7,1/7,2/7,1/7,1/30,0],
#    [2/7,1/7,2/7,1/7,1/7,1/7,1/30,0],
#    [2/7,1/7,1/7,2/7,1/7,1/7,1/30,0],
#    [2/7,1/7,1/7,1/7,2/7,1/7,1/30,0],
#    [2/7,1/7,1/7,1/7,1/7,2/7,1/30,0],
#    [2/7,2/7,2/7,1/7,2/7,1/7,1/30,0],
#    [2/7,2/7,1/7,2/7,2/7,1/7,1/30,0]]

    equal_test_policies = list()

    centre_points_list = list()
    centre_equal_points_list = list()
    pess_points_list = list()
    pess_equal_points_list = list()

    for test_policy in test_policies:
        centre_points_list.append(get_centre_point())
        centre_equal_points_list.append(get_centre_point())
        pess_points_list.append(res_pessimistic.copy())
        pess_equal_points_list.append(res_pessimistic.copy())

        params_list = uncertainty_point_to_params_dict(get_centre_point())[0][0]
        equal_test_policies.append(equal_freq_policy(get_tests_per_day(test_policy, params_list), params_list))

    centre_folder = './jun_23_sims/test_pareto_centre_{}/'.format(get_timestamp())
    centre_equal_folder = './jun_23_sims/test_pareto_centre_equal_{}/'.format(get_timestamp())
    pess_folder = './jun_23_sims/test_pareto_pess_{}/'.format(get_timestamp())
    pess_equal_folder = './jun_23_sims/test_pareto_pess_equal_{}/'.format(get_timestamp())

    os.mkdir(centre_folder)
    os.mkdir(centre_equal_folder)
    os.mkdir(pess_folder)
    os.mkdir(pess_equal_folder)

    centre_fnames = [centre_folder + 'point_{}.dill'.format(count) for count in \
                   range(len(centre_points_list))]
    centre_equal_fnames = [centre_equal_folder + 'point_{}.dill'.format(count) for count in \
                   range(len(centre_equal_points_list))]
    pess_fnames = [pess_folder + 'point_{}.dill'.format(count) for count in range(len(pess_points_list))]
    pess_equal_fnames = [pess_equal_folder + 'point_{}.dill'.format(count) for count in range(len(pess_equal_points_list))]

    processes = []

    processes.extend(run_sims_new_process(centre_points_list, test_policies, centre_fnames, nreps=100))
    processes.extend(run_sims_new_process(pess_points_list, test_policies, pess_fnames, nreps=100))
    print(equal_test_policies[0])
    processes.extend(run_sims_new_process(centre_equal_points_list, equal_test_policies, centre_equal_fnames, nreps=100))
    processes.extend(run_sims_new_process(pess_equal_points_list, equal_test_policies, pess_equal_fnames, nreps=100))

    file = open(centre_folder + 'test_policies.dill', mode='wb')
    dill.dump([test_policies], file)
    file.close()

    file = open(centre_equal_folder + 'test_policies.dill', mode='wb')
    dill.dump([equal_test_policies], file)
    file.close()

    file = open(pess_folder + 'test_policies.dill', mode='wb')
    dill.dump([test_policies], file)
    file.close()

    file = open(pess_equal_folder + 'test_policies.dill', mode='wb')
    dill.dump([equal_test_policies], file)
    file.close()

    print("finished launching processes, waiting for them to finish")
    for p in processes:
        p.join()
        

