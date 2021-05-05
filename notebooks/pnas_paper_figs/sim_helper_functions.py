import numpy as np
from uncertainty_analysis import uncertainty_point_to_params_dict
from multiprocessing import Process
import dill

import os
import sys
module_path = os.path.abspath(os.path.join('../..'))
if module_path not in sys.path:
        sys.path.append(module_path + "/src/simulations_v2")
from load_params import load_params, update_sev_prevalence
from analysis_helpers import poisson_waiting_function
from multi_group_simulation import MultiGroupSimulation
from util_functions import *


def get_point_on_line(centre, pess, mult, direction = None):
    if direction == None:
        direction = {}
        for param in centre:
            direction[param] = pess[param] - centre[param]

    point = {}
    for param in centre:
        point[param] = centre[param] + mult * direction[param]

    return point


def get_points_on_line(centre, pess, npoints=13, mult_lb=-0.1, mult_ub=1.1):
    mults = np.linspace(mult_lb, mult_ub, npoints)
    direction = {}
    for param in centre:
        direction[param] = pess[param] - centre[param]

    return [get_point_on_line(centre, pess, mult, direction) for mult in mults]


def run_sensitivity_analysis_sims(centre, pess, param_to_vary, output_folder, npoints=13, mult_ub=1.1, mult_lb = -0.1, nreps=50):
    mults = np.linspace(mult_lb, mult_ub, npoints)
    output_fnames = [output_folder + "/param_to_vary_{}__mult_{}".format(param_to_vary, mult) for mult in mults]

    sim_points = get_points_on_line(centre, pess, npoints, mult_lb, mult_ub)
    run_sims_new_process(sim_points, output_fnames, nreps=nreps, run_residential_only=True)

    


def run_new_process(uncertainty_point, filename, point_id, nreps=50, T=112, run_only_residential=True):
    p = Process(target = run_simulation, args = (uncertainty_point, filename, point_id, nreps, T, run_only_residential))
    p.start()
    return p

def run_sims_new_process(uncertainty_point_dicts, output_fnames, nreps=50, T=112, run_only_residential=True, wait_for_processes_to_join=True):
    idx = 0
    processes = []
    for uncertainty_dict, output_fname in zip(uncertainty_point_dicts, output_fnames):
        idx += 1
        p = run_new_process(uncertainty_dict, output_fname, idx, nreps=nreps, T=T)
        processes.append(p)

    print("launched {} processes".format(len(processes)))
    if wait_for_processes_to_join:
        for p in processes:
            p.join()
        print("done running processes")
        return processes
    else:
        return processes


def run_simulation(uncertainty_point, 
                    filename, 
                    point_id=None,
                    nreps=50, T=112, 
                    run_only_residential=False):
    # get params
    (res_params_list, res_interaction_matrix, res_group_names),\
            (virtual_params_list, virtual_interaction_matrix, virtual_group_names) \
            = uncertainty_point_to_params_dict(uncertainty_point)

    print("running sim with id {}, nreps {}, T {}, filename {}, run_only_residential {}".format(point_id, nreps, T, filename, run_only_residential))
    # run simulations
    # Residential Simulation
    res_test_policy = [2/7,2/7,1/7,1/7,2/7,1/7,1/30,0]
    virtual_test_policy = [0, 2/7,1/7,0,1/7, 2/7,1/7,1/30, 0]
    
    # Running res sims
    print('running residential sim with id {}'.format(point_id))
    res_tests_per_day, res_inf_matrix, res_hosp_matrix = evaluate_testing_policy(res_params_list, res_interaction_matrix, res_group_names, res_test_policy, T, nreps)
    
    # Running virtual sims
    if not run_only_residential:
        print('running virtual sim with id {}'.format(point_id))
        virtual_tests_per_day, virtual_inf_matrix, virtual_hosp_matrix = evaluate_testing_policy(virtual_params_list, virtual_interaction_matrix, virtual_group_names, virtual_test_policy, T, nreps)
    else:
        virtual_tests_per_day, virtual_inf_matrix, virtual_hosp_matrix = (None, None, None)
    
    # save output
    file = open(filename, mode='wb')
    dill.dump([uncertainty_point, res_inf_matrix, res_hosp_matrix, virtual_inf_matrix, virtual_hosp_matrix], file)
    file.close()    
    return


