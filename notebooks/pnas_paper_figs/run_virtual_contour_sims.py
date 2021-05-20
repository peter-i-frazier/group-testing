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


def run_virtual_contour_plot_sims(x_variable = 'virtual_pop_size', y_variable = 'virtual_noncompliance', mult_lb=1, mult_ub=1, npoints=13, nreps=50):

    assert x_variable in get_centre_point().keys()
    assert y_variable in get_centre_point().keys()

    x_values = np.linspace(PARAM_BOUNDS[x_variable][0]*mult_lb, PARAM_BOUNDS[x_variable][1]*mult_ub, npoints)
    y_values = np.linspace(PARAM_BOUNDS[y_variable][0]*mult_lb, PARAM_BOUNDS[y_variable][1]*mult_ub, npoints)

    folder_name = "./virtual_contour_plot_sims_{}/".format(get_timestamp())
    os.mkdir(folder_name)
    
    uncertainty_points = []
    filenames = []
    for x_val in x_values:
        for y_val in y_values:
            point = get_centre_point()
            point[x_variable] = x_val
            point[y_variable] = y_val
            uncertainty_points.append(point)
            fname = folder_name + "{}_{}_{}_{}.dill".format(
                                    x_variable, x_val, y_variable, y_val)
            filenames.append(fname)

    processes = run_sims_new_process(uncertainty_points, filenames, 
                    nreps=nreps,
                    wait_for_processes_to_join=False)
    
    return processes
    

if __name__ == "__main__":
    
    nreps=50
    npoints=13
    procs = run_virtual_contour_plot_sims(x_variable = 'virtual_pop_size', y_variable = 'virtual_noncompliance',
            npoints=npoints, nreps=nreps)

    for p in procs:
        p.join()

