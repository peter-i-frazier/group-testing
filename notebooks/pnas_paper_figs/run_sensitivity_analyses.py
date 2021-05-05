import sys
import os
import numpy as np
import multiprocessing
import dill
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


from util_functions import *
from uncertainty_analysis import *
from sim_helper_functions import *


def launch_sensitivity_analysis(centre, pess, param, base_folder, nreps=50):
    param_lb = min(centre[param], pess[param])
    param_ub = max(centre[param], pess[param])

    centre_lb = centre.copy()
    centre_lb[param] = param_lb
    
    centre_ub = centre.copy()
    centre_ub[param] = param_ub

    centre_points = get_points_on_line(centre_lb, centre_ub, mult_lb=-1.1, mult_ub=1.1)

    pess_lb = pess.copy()
    pess_lb[param] = param_lb
    pess_ub = pess.copy()
    pess_ub[param] = param_ub

    pess_points = get_points_on_line(pess_lb, pess_ub, mult_lb=-1.1, mult_ub=1.1)

    folder_name = '{}/{}/'.format(base_folder, param)
    os.mkdir(folder_name)

    centre_fnames = [folder_name + 'centre_mult_{}.dill'.format(mult) for mult in \
                   range(len(centre_points))] 

    pess_fnames = [folder_name + 'pess_mult_{}.dill'.format(mult) for mult in \
                   range(len(pess_points))]

    uncertainty_points = centre_points + pess_points
    fnames = centre_fnames + pess_fnames
    #import pdb; pdb.set_trace()
    processes = run_sims_new_process(uncertainty_points, fnames, nreps=nreps, run_only_residential=True,
            wait_for_processes_to_join=False)
    return processes


    

if __name__ == "__main__":
    lhs_output_sim_files = []
    for i in range(2000):
        fname = '/home/jmc678/covid_data/group-testing/notebooks/apr_29_scenarios/point_{}.dill'.format(i)
        lhs_output_sim_files.append(fname)


    scenario_data = load_sim_output(lhs_output_sim_files)
    res_results = residential_regression(scenario_data)
    res_pessimistic = calculate_pessimistic_scenario(res_results)
    centre = get_centre_point()

    base_folder = './sensitivity_sims_{}/'.format(get_timestamp())
    os.mkdir(base_folder)

    processes = []

    for param in UNCERTAINTY_PARAMS_LIST:
        processes.extend(launch_sensitivity_analysis(centre, res_pessimistic, param, base_folder,nreps=50))

    print("finished launching processes, waiting for them to finish")
    for p in processes:
        p.join()
        

