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


def run_contour_plot_sims(centre, pess, mult_lb, mult_ub, 
        test_policy_mult_lb, test_policy_mult_ub, npoints=13, nreps=50):

    x_axis_base_points = get_points_on_line(centre, pess, mult_ub=mult_ub, 
            mult_lb=mult_lb, npoints=npoints)

    test_policy_mults = np.linspace(test_policy_mult_lb, test_policy_mult_ub, 
                                    npoints)

    folder_name = "jun_23_sims/scenario_test_freq_contour_{}/".format(get_timestamp())
    os.mkdir(folder_name)
    
    x_axis_indices = range(len(x_axis_base_points))

    uncertainty_points = []
    test_mults_for_sim_runs = []
    filenames = []
    for x_idx in x_axis_indices:
        for test_policy_mult in test_policy_mults:
            point = x_axis_base_points[x_idx].copy()
            uncertainty_points.append(point)
            test_mults_for_sim_runs.append(test_policy_mult)
            fname = folder_name + "test_mult_{}_x_axis_point_idx_{}.dill".format(
                                    x_idx, test_policy_mult)
            filenames.append(fname)

    processes = run_sims_new_process(uncertainty_points, filenames, 
                    test_policy_multipliers=test_mults_for_sim_runs,
                    nreps=nreps,
                    wait_for_processes_to_join=False)
    
    return processes
    

if __name__ == "__main__":
    lhs_output_sim_files = []
    for i in range(2000):
        fname = '/home/aaj54/group-testing/notebooks/apr_29_scenarios/point_{}.dill'.format(i)
        lhs_output_sim_files.append(fname)

    scenario_data = load_sim_output(lhs_output_sim_files)
    res_results = residential_regression(scenario_data)
    res_pessimistic = calculate_pessimistic_scenario(res_results)
    centre = get_centre_point()

    nreps=100
    npoints=13
    procs = run_contour_plot_sims(centre, res_pessimistic, mult_lb=-1.1, mult_ub=1.1,
            test_policy_mult_lb=0.5, test_policy_mult_ub=1.5, npoints=npoints,
            nreps=nreps)

    for p in procs:
        p.join()

