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


def run_contour_plot_sims(base_point, x_variable = 'virtual_pop_size', y_variable = 'virtual_noncompliance', base_folder_name='./virtual_contour_plot_sims', x_lb=None, x_ub=None, y_lb=None, y_ub=None, npoints=13, nreps=50):

    assert x_variable in get_centre_point().keys()
    assert y_variable in get_centre_point().keys()

    if x_lb == None:
        x_lb = PARAM_BOUNDS[x_variable][0]
    if x_ub == None:
        x_ub = PARAM_BOUNDS[x_variable][1]
    if y_lb == None:
        y_lb = PARAM_BOUNDS[y_variable][0]
    if y_ub == None:
        y_ub = PARAM_BOUNDS[y_variable][1]

    x_values = np.linspace(x_lb, x_ub, npoints)
    y_values = np.linspace(y_lb, y_ub, npoints)

    folder_name = "{}_{}/".format(base_folder_name, get_timestamp())
    os.mkdir(folder_name)
    
    uncertainty_points = []
    filenames = []
    for x_val in x_values:
        for y_val in y_values:
            point = base_point.copy()
            point[x_variable] = x_val
            point[y_variable] = y_val
            uncertainty_points.append(point)
            fname = folder_name + "{}_{}_{}_{}.dill".format(
                                    x_variable, x_val, y_variable, y_val)
            filenames.append(fname)

    processes = run_sims_new_process(uncertainty_points, filenames, 
                    nreps=nreps,
                    wait_for_processes_to_join=False, run_only_residential=False)
    
    return processes
    

if __name__ == "__main__":

    lhs_output_sim_files = []
    for i in range(2000):
        fname = '/home/aaj54/group-testing/notebooks/apr_24_scenarios/point_{}.dill'.format(i)
        lhs_output_sim_files.append(fname)

    scenario_data = load_sim_output(lhs_output_sim_files)
    res_results = residential_regression(scenario_data)
    res_pessimistic = calculate_pessimistic_scenario(res_results)
    centre = get_centre_point()
    
    nreps=100
    npoints=13
    
    # procs = run_contour_plot_sims(centre, x_variable = 'virtual_pop_size', y_variable = 'virtual_noncompliance',
    #         base_folder_name='./virtual_contour_plot_sims', npoints=npoints, nreps=nreps)

    procs = run_contour_plot_sims(res_pessimistic, x_variable = 'test_noncompliance', y_variable = 'test_sensitivity',
            base_folder_name='./test_comp_sens_pess_contour_sims', x_lb = 0.01, x_ub = 0.5, npoints=npoints, nreps=nreps)

    for p in procs:
        p.join()

