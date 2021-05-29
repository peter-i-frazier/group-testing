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


def get_direction(pess, centre):
    direction = dict()
    for param in centre.keys():
        direction[param] = pess[param] - centre[param]
    return direction

def generate_new_params(centre, direction, mult):
    new_params = dict()
    for param in centre.keys():
        new_params[param] = centre[param] + mult * direction[param]
    return new_params

def launch_sensitivity_analysis(centre, direction, base_folder, mult_list=np.linspace(-1.1, 1.1, 23), nreps=50):

    points = [generate_new_params(centre, direction, mult) for mult in mult_list]

    folder_name = '{}/{}/'.format(base_folder, 'pess_sensitivity')
    os.mkdir(folder_name)

    fnames = [folder_name + 'mult_{}.dill'.format(mult) for mult in mult_list] 

    #import pdb; pdb.set_trace()
    processes = run_sims_new_process(points, fnames, nreps=nreps, run_only_residential=False,
            wait_for_processes_to_join=False)
    return processes


MULT_RANGE = np.linspace(-1.1, 1.1, 23)    
if __name__ == "__main__":
    lhs_output_sim_files = []
    for i in range(2000):
        fname = '/home/jmc678/covid_data/group-testing/notebooks/apr_29_scenarios/point_{}.dill'.format(i)
        lhs_output_sim_files.append(fname)

    scenario_data = load_sim_output(lhs_output_sim_files)
    res_results = residential_regression(scenario_data)
    #res_results = virtual_vs_residential_regression(scenario_data)
    res_pessimistic = calculate_pessimistic_scenario(res_results)
    centre = get_centre_point()

    direction = get_direction(res_pessimistic, centre)

    # toggle whether/not to use poisson contact tracing
    os.environ['use_poisson_contact_tracing'] = 'True'

    base_folder = './may_29_sims/pess_res_sensitivity_sims_{}/'.format(get_timestamp())
    os.mkdir(base_folder)

    processes = []

    for mult in MULT_RANGE:
    #for param in PARAMS_LIST:
        processes.extend(launch_sensitivity_analysis(centre, direction, base_folder, mult_list = MULT_RANGE, nreps=100))

    print("finished launching processes, waiting for them to finish")
    for p in processes:
        p.join()
        

