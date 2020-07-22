import sys
import itertools
import numpy as np
import yaml
import time
import os
import multiprocessing
from analysis_helpers import run_multiple_trajectories
import dill
import argparse
from load_params import load_params
from plotting_util import plot_from_folder

BASE_DIRECTORY= os.path.abspath(os.path.join('')) + "/sim_output/"

VALID_PARAMS_TO_VARY = [
    'contact_tracing_isolations',
    'expected_contacts_per_day',
    'symptomatic_daily_self_report_p',
    'asymptomatic_daily_self_report_p',
    'initial_ID_prevalence',
    'exposed_infection_p',
    'asymptomatic_p',
    'contact_tracing_delay',
    'test_protocol_QFNR',
    'test_population_fraction',
    'daily_outside_infection_p',
    'initial_E_count',
    'initial_pre_ID_count',
    'initial_ID_count',
    'initial_SyID_mild_count',
    'initial_SyID_severe_count',
    'population_size'
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

    elif param_to_vary == 'symptomatic_daily_self_report_p':
        sim_params['severe_symptoms_daily_self_report_p'] = param_val
    elif param_to_vary == 'asymptomatic_daily_self_report_p':
        sim_params['mild_symptoms_daily_self_report_p'] = param_val


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

def iter_param_variations(base_params, params_to_vary, param_values):
    # iterator that generates all parameter configurations corresponding to
    # all combinations of parameter values across the different params_to_vary.
    # each return value is a tuple (param_specifier, params) where params is the parameter
    # dictionary object, and param_specifier is a smaller dict specifying the varying
    # params and the value they are taking righ tnow 
    base_params = base_params.copy()
    params_list = [param_values[param] for param in params_to_vary]
    for param_tuple in itertools.product(*params_list):
        param_specifier = {}
        for param, value in zip(params_to_vary, param_tuple):
            update_params(base_params, param, value)
            param_specifier[param] = value
            
        yield param_specifier, base_params


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run multiple simulations using multiprocessing')
    parser.add_argument('-o', '--outputdir', default=BASE_DIRECTORY, 
                        help='directory to store simulation output')
    parser.add_argument('-V', '--verbose', action='store_true', help='include verbose output')
    parser.add_argument('-s', '--scenarios', nargs='+', required=True,
            help='list of YAML config files specifying base sets of scenario parameters to use')

    parser.add_argument('-p', '--param-to-vary', action='append',
            help='which param(s) should be varied in the corresponding sensitivity sims', required=True)
    
    parser.add_argument('-v', '--values', required=True, nargs='+', action='append',
            help='what values should the varying parameter(s) take')

    parser.add_argument('-n', '--ntrajectories', default=500,
            help='how many trajectories to simulate for each (scenario, value) pair')
    parser.add_argument('-t', '--time-horizon', default=112,
            help='how many days to simulate for each trajectory')

    parser.add_argument('-f', '--fig-dir',
            help='specify folder where plots should be saved')

    args = parser.parse_args()

    if len(args.values) != len(args.param_to_vary):
        raise(Exception("Number of parameters specified doesn't match number of value ranges specified"))

    scenarios = {}
    for scenario_file in args.scenarios:
        scn_name, scn_params = load_params(scenario_file)
        scenarios[scn_name] = scn_params

    params_to_vary = args.param_to_vary
    param_values = {}
    for param_to_vary, values in zip(params_to_vary, args.values):
        if param_to_vary not in VALID_PARAMS_TO_VARY:
            print("Received invalid parameter to vary: {}".format(param_to_vary))
            exit()
        if param_to_vary == 'contact_tracing_delay':
            param_values[param_to_vary] = [int(v) for v in values]
        else:
            param_values[param_to_vary] = [float(v) for v in values]

    

    if len(params_to_vary) == 1:
        sim_id = "{timestamp}-{param_to_vary}".format(
                    timestamp=str(time.time()), 
                    param_to_vary = params_to_vary[0])
    else:
        sim_id = "{timestamp}-multiparam".format(
                    timestamp=str(time.time()).split('.')[0]) 
    
    print("Using Simulation ID: {}".format(sim_id))

    
    basedir = args.outputdir
    ntrajectories = int(args.ntrajectories)
    time_horizon = int(args.time_horizon)
    verbose = args.verbose

    if not os.path.isdir(basedir):
        print("Directory {} does not exist. Please create it.".format(basedir))
        exit()

    sim_main_dir = basedir + "/" + str(sim_id)
    os.mkdir(sim_main_dir)
    print("Output directory {} created".format(sim_main_dir))


   
    jobs = []
    scn_dirs = {} 
    for scn_name, scn_params in scenarios.items():
        sim_scn_dir = sim_main_dir + "/" + scn_name
        os.mkdir(sim_scn_dir)
        scn_dirs[scn_name] = sim_scn_dir
        dill.dump(scn_params, open("{}/scn_params.dill".format(sim_scn_dir), "wb"))
        job_counter = 0
        if verbose:
            print("Output scenario-directory {} created".format(sim_sub_dir))

        for param_specifier, sim_params in iter_param_variations(scn_params, params_to_vary, param_values):
            # create the relevant subdirectory
            sim_sub_dir = "{}/simulation-{}".format(sim_scn_dir, job_counter)
            job_counter += 1
            os.mkdir(sim_sub_dir)
            with open('{}/param_specifier.yaml'.format(sim_sub_dir), 'w') as outfile:
                yaml.dump(param_specifier, outfile, default_flow_style=False)
            if verbose: 
                print("Created directory {} to save output".format(sim_sub_dir))

            
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

    print("Running simulations for {} scenarios and {} varying parameters across {} separate processes.".format(len(scenarios), len(params_to_vary), len(jobs)))
    print("Results being saved in output directory {}.".format(sim_main_dir))
    print("Waiting for simulations to finish...")
    for p in jobs:
        p.join()

            
    if len(params_to_vary) > 1:
        print("Simulations done. Not auto-generating plots because > 1 parameter was varied")
        print("Exiting now...")
        exit()
    print("Simulations done. Generating plots now...")
    if args.fig_dir == None:
        fig_dir = sim_main_dir
    else:
        fig_dir = args.fig_dir
    plot_from_folder(sim_main_dir, fig_dir)
    print("Saved plots to directory {}".format(fig_dir))







