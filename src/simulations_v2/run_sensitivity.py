import sys
import numpy as np
import yaml
import time
import os
import multiprocessing
from analysis_helpers import run_multiple_trajectories
import dill
import argparse
from load_params import load_params
from plotting_util import plot_from_folders

BASE_DIRECTORY="/nfs01/covid_sims/"

VALID_PARAMS_TO_VARY = [
    'contact_tracing_isolations'
    'expected_contacts_per_day',
    'mild_symptoms_daily_self_report_p',
    'severe_symptoms_daily_self_report_p',
    'initial_ID_prevalence',
    'exposed_infection_p',
    'asymptomatic_p',
    'contact_tracing_delay',
    'test_protocol_QFNR',
    'test_population_fraction',
    'daily_outside_infection_p'
    ]

def run_background_sim(output_dir, sim_params, ntrajectories=150, time_horizon=112):
    dfs = run_multiple_trajectories(sim_params, ntrajectories, time_horizon)
    
    # record output
    for idx, df in enumerate(dfs):
        df_file_name = "{}/{}.csv".format(output_dir, idx)
        df.to_csv(df_file_name)


def update_params(sim_params, param_to_vary, param_val):
    # VERY TEMPORARY HACK TO GET SENSITIVITY SIMS FOR ASYMPTOMATIC %
    if param_to_vary == 'asymptomatic_p':
        assert(sim_params['mild_severity_levels'] == 1)
        curr_prevalence_dist = sim_params['severity_prevalence']
        assert(param_val >= 0 and param_val <= 1)
        new_dist = [param_val]
        remaining_mass = sum(curr_prevalence_dist[1:])

        # need to scale so that param_val + x * remaning_mass == 1
        scale = (1 - param_val) / remaining_mass
        idx = 1
        while idx < len(curr_prevalence_dist):
            new_dist.append(curr_prevalence_dist[idx] * scale)
            idx += 1
        assert(np.isclose(sum(new_dist), 1))
        sim_params['severity_prevalence'] = np.array(new_dist)

    # VERY TEMPORARY HACK TO GET SENSITIVITY SIMS WORKING FOR CONTACT RECALL %
    elif param_to_vary == 'contact_tracing_constant':
        num_isolations = sim_params['cases_isolated_per_contact']
        base_recall = sim_params['contact_recall']
        new_isolations = num_isolations * param_val / base_recall
        new_quarantines = max(7 - new_isolations, 0)
        sim_params['cases_isolated_per_contact'] = new_isolations
        sim_params['cases_quarantined_per_contact'] = new_quarantines
    elif param_to_vary == 'contact_tracing_isolations':
        sim_params['cases_isolated_per_contact'] = param_val
        sim_params['cases_quarantined_per_contact'] = 7 - param_val
    else:
        sim_params[param_to_vary] = param_val


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run multiple simulations using multiprocessing')
    parser.add_argument('-o', '--outputdir', default=BASE_DIRECTORY, 
                        help='directory to store simulation output')
    parser.add_argument('-V', '--verbose', action='store_true', help='include verbose output')
    parser.add_argument('-s', '--scenarios', nargs='+', required=True,
            help='list of YAML config files specifying base sets of scenario parameters to use')

    # TODO: make a list somewhere of all allowable values of param-to-vary
    parser.add_argument('-p', '--param-to-vary', 
            help='which param should be varied in the corresponding sensitivity sims', required=True)
    
    parser.add_argument('-v', '--values', required=True, nargs='+',
            help='what values should the varying parameter take')

    parser.add_argument('-n', '--ntrajectories', default=500,
            help='how many trajectories to simulate for each (scenario, value) pair')
    parser.add_argument('-t', '--time-horizon', default=112,
            help='how many days to simulate for each trajectory')

    parser.add_argument('-f', '--fig-dir',
            help='specify folder where plots should be saved')

    args = parser.parse_args()

    scenarios = {}
    for scenario_file in args.scenarios:
        scn_name, scn_params = load_params(scenario_file)
        scenarios[scn_name] = scn_params

    param_to_vary = args.param_to_vary
    if param_to_vary not in VALID_PARAMS_TO_VARY:
        print("Received invalid parameter to vary: {}".format(param_to_vary))
        exit()
    param_values = [float(v) for v in args.values]

    timestamp = time.time()
    sim_id = "{timestamp}-{param_to_vary}".format(
                timestamp=str(time.time()).split('.')[0], 
                param_to_vary = param_to_vary)
    
    basedir = args.outputdir
    ntrajectories = int(args.ntrajectories)
    time_horizon = int(args.time_horizon)
    verbose = args.verbose
    if verbose:
        print("Simulation ID: {}".format(sim_id))

    if not os.path.isdir(basedir):
        print("Directory {} does not exist. Please create it.".format(basedir))
        exit()

    sim_main_dir = basedir + "/" + str(sim_id)
    os.mkdir(sim_main_dir)
    if verbose:
        print("Output directory {} created".format(sim_main_dir))


   
    jobs = []
    scn_dirs = {} 
    for scn_name in scenarios:
        sim_scn_dir = sim_main_dir + "/" + scn_name
        os.mkdir(sim_scn_dir)
        scn_dirs[scn_name] = sim_scn_dir

        if verbose:
            print("Output scenario-directory {} created".format(sim_sub_dir))

        for param_val in param_values:
            # create the relevant subdirectory
            sim_sub_dir = "{}/{}-{}".format(sim_scn_dir, param_to_vary, param_val)
            os.mkdir(sim_sub_dir)
            if verbose: 
                print("Created directory {} to save output".format(sim_sub_dir))

            # instantiate relevant sim params
            sim_params = scenarios[scn_name].copy()
            
            update_params(sim_params, param_to_vary, param_val)        
            
            dill.dump(sim_params, open("{}/sim_params.dill".format(sim_sub_dir), "wb"))
            # start new process
            fn_args = (sim_sub_dir, sim_params, ntrajectories, time_horizon)
            #run_background_sim(*fn_args)
            proc = multiprocessing.Process(target = run_background_sim, args=fn_args)
            #proc.daemon = True
            jobs.append(proc)
            proc.start()
            if verbose:
                print("starting process for {} value {}".format(param_to_vary, param_val))
        
    print("Running simulations for {} scenarios and {} parameter values across {} separate processes.".format(len(scenarios), len(param_values), len(jobs)))
    print("Results being saved in output directory {}.".format(sim_main_dir))
    print("Waiting for simulations to finish...")
    for p in jobs:
        p.join()

            
    print("Simulations done. Generating plots now...")
    if args.fig_dir == None:
        fig_dir = sim_main_dir
    else:
        fig_dir = args.fig_dir
    plot_from_folders(scn_dirs, param_to_vary, fig_dir)
    print("Saved plots to directory {}".format(fig_dir))







