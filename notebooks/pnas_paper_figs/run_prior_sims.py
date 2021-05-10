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

def sample_from_prior():
    mu = list()
    var = list()
    param_order = list()
    for param, bounds in PARAM_BOUNDS.items():
        param_order.append(param)
        mu.append(np.mean(bounds))
        sd = (bounds[1] - np.mean(bounds)) / norm.ppf(0.975)
        var.append(sd ** 2)
    point = np.random.multivariate_normal(np.array(mu), np.diag(var))
    point_df = dict()
    for index, value in enumerate(point):
        point_df[param_order[index]] = value
        if value < 0:
            print(param_order[index], value)
    return point_df

def launch_prior_sims(npoints = 200, nreps=50):

    points = [sample_from_prior() for _ in range(npoints)]

    for point in points:
        print(point)

    timestamp = get_timestamp()
    folder_name = './prior_sims/timestamp_{}/'.format(timestamp)
    os.mkdir(folder_name)

    fnames = [folder_name + 'prior_{}.dill'.format(count) for count in \
                   range(len(points))] 

    processes = run_sims_new_process(points, fnames, nreps=nreps, run_only_residential=False,
            wait_for_processes_to_join=False)
    return processes


    

if __name__ == "__main__":
    np.random.seed(2021)
    processes = []

    processes.extend(launch_prior_sims(npoints = 200, nreps=50))

    print("finished launching processes, waiting for them to finish")
    for p in processes:
        p.join()
        

