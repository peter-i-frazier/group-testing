import sys
import numpy as np
import yaml
import time
import os
import multiprocessing

from fall_realistic import base_params as realistic_params
from fall_realistic import base_params_testing as realistic_params_testing

from fall_slightly_optimistic import base_params as optimistic_params
from fall_slightly_optimistic import base_params_testing as optimistic_params_testing

from fall_slightly_pessimistic import base_params as pessimistic_params
from fall_slightly_pessimistic import base_params_testing as pessimistic_params_testing



from analysis_helpers import run_multiple_trajectories
#import dill

BASE_DIRECTORY="/nfs01/covid_sims/"

def run_background_sim(output_dir, sim_params, ntrajectories=150, time_horizon=112):
    try:
        dfs = run_multiple_trajectories(sim_params, ntrajectories, time_horizon)
        
        # record output
        for idx, df in enumerate(dfs):
            df_file_name = "{}/{}.csv".format(output_dir, idx)
            df.to_csv(df_file_name)
    except Exception as e:
        error_msg = "Encountered error: {}".format(str(e))
        print(error_msg)
        f = open(output_dir + "/error.txt", "w")
        f.write(error_msg)
        f.close()

if __name__ == "__main__":
    sim_timestamp = time.time()
    sim_id = "nominal_sims.{}".format(sim_timestamp)
    print("Simulation ID: {}".format(sim_id))
    sim_main_dir = BASE_DIRECTORY + str(sim_id)
    os.mkdir(sim_main_dir)
    print("created main directory {}".format(sim_main_dir))

    sub_dir = "{}/{}".format(sim_main_dir, "realistic_params")
    os.mkdir(sub_dir)
    fn_args = (sub_dir, realistic_params, 500, 112)
    proc = multiprocessing.Process(target = run_background_sim, args=fn_args)
    proc.start()

    sub_dir = "{}/{}".format(sim_main_dir, "realistic_params_testing")
    os.mkdir(sub_dir)
    fn_args = (sub_dir, realistic_params_testing, 500, 112)
    proc = multiprocessing.Process(target = run_background_sim, args=fn_args)
    proc.start()

    sub_dir = "{}/{}".format(sim_main_dir, "optimistic_params_testing")
    os.mkdir(sub_dir)
    fn_args = (sub_dir, optimistic_params_testing, 500, 112)
    proc = multiprocessing.Process(target = run_background_sim, args=fn_args)
    proc.start()

    sub_dir = "{}/{}".format(sim_main_dir, "optimistic_params")
    os.mkdir(sub_dir)
    fn_args = (sub_dir, optimistic_params, 500, 112)
    proc = multiprocessing.Process(target = run_background_sim, args=fn_args)
    proc.start()

    sub_dir = "{}/{}".format(sim_main_dir, "pessimistic_params_testing")
    os.mkdir(sub_dir)
    fn_args = (sub_dir, pessimistic_params_testing, 500, 112)
    proc = multiprocessing.Process(target = run_background_sim, args=fn_args)
    proc.start()

    sub_dir = "{}/{}".format(sim_main_dir, "pessimistic_params")
    os.mkdir(sub_dir)
    fn_args = (sub_dir, pessimistic_params, 500, 112)
    proc = multiprocessing.Process(target = run_background_sim, args=fn_args)
    proc.start()
