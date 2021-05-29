import sys
import os
import numpy as np
import multiprocessing
import dill
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import norm

from util_functions import *
from uncertainty_analysis import *
from sim_helper_functions import *

def run_simulation(uncertainty_point,
                    test_policy,
                    filename,
                    point_id=None,
                    nreps=50, T=112):
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


def run_new_process(uncertainty_point, test_policy, filename, point_id, nreps=50, T=112):
    p = Process(target = run_simulation, args = (uncertainty_point, test_policy, filename, point_id, nreps, T))
    p.start()
    return p


def run_sims_new_process(uncertainty_point_dicts, test_policies, output_fnames, nreps=50, T=112, wait_for_processes_to_join=True):
    idx = 0
    processes = []
    for uncertainty_dict, output_fname, test_policy in zip(uncertainty_point_dicts, output_fnames, test_policies):
        idx += 1
        p = run_new_process(uncertainty_dict, test_policy, output_fname, idx, nreps=nreps, T=T)
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
        fname = '/home/aaj54/group-testing/notebooks/apr_24_scenarios/point_{}.dill'.format(i)
        lhs_output_sim_files.append(fname)

    scenario_data = load_sim_output(lhs_output_sim_files)
    res_results = residential_regression(scenario_data)
    res_pessimistic = calculate_pessimistic_scenario(res_results)

    test_policies = [[1/7,1/7,1/7,1/7,1/7,1/7,1/14,0],
    [2/7,1/7,1/7,1/7,1/7,1/7,1/30,0],
    [2/7,2/7,1/7,1/7,1/7,1/7,1/30,0],
    [2/7,2/7,1/7,1/7,2/7,1/7,1/30,0],
    [2/7,1/7,2/7,1/7,1/7,1/7,1/30,0],
    [2/7,1/7,1/7,2/7,1/7,1/7,1/30,0],
    [2/7,1/7,1/7,1/7,2/7,1/7,1/30,0],
    [2/7,1/7,1/7,1/7,1/7,2/7,1/30,0],
    [2/7,2/7,2/7,1/7,2/7,1/7,1/30,0],
    [2/7,2/7,1/7,2/7,2/7,1/7,1/30,0]]

    centre_points_list = list()
    pess_points_list = list()
    for _ in test_policies:
        centre_points_list.append(get_centre_point())
        pess_points_list.append(res_pessimistic.copy())

    centre_folder = './test_pareto_centre_{}/'.format(get_timestamp())
    pess_folder = './test_pareto_pess_{}/'.format(get_timestamp())

    os.mkdir(centre_folder)
    os.mkdir(pess_folder)

    centre_fnames = [centre_folder + 'point_{}.dill'.format(count) for count in \
                   range(len(centre_points_list))]

    pess_fnames = [pess_folder + 'point_{}.dill'.format(count) for count in range(len(pess_points_list))]

    processes = []

    processes.extend(run_sims_new_process(centre_points_list, test_policies, centre_fnames, nreps=50))
    processes.extend(run_sims_new_process(pess_points_list, test_policies, pess_fnames, nreps=50))

    print("finished launching processes, waiting for them to finish")
    for p in processes:
        p.join()
        

