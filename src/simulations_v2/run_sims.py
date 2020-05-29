import sys
import numpy as np
import yaml
import time
import os
import multiprocessing
from analysis_helpers import run_multiple_trajectories
import dill
from params import ParamConfig
import argparse

BASE_DIRECTORY="/nfs01/covid_sims/"

def run_background_sim(output_dir, sim_params, ntrajectories=150, time_horizon=112):
    try:
        dfs = run_multiple_trajectories(sim_params, ntrajectories, time_horizon)
        
        # record output
        for idx, df in enumerate(dfs):
            df_file_name = "{}/{}.csv".format(output_dir, idx)
            df.to_csv(df_file_name)
    except Exception as e:
        print(e)
        error_msg = "Encountered error: {}".format(str(e))
        print(error_msg)
        f = open(output_dir + "/error.txt", "w")
        f.write(error_msg)
        f.close()


def update_params(sim_params, param_to_vary, param_val):
    # VERY TEMPORARY HACK TO GET SENSITIVITY SIMS FOR ASYMPTOMATIC %
    if param_to_vary == 'asymptomatic_p':
        assert(sim_params['mild_severity_levels'] == 1)
        curr_prevalence_dist = params['severity_prevalence']
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
    parser.add_argument('-w', '--when', choices=['fall', 'june'], 
                        default='fall',
                        help='what time period should the simulations use')
    parser.add_argument('-a', '--assumption', choices=['pessimistic', 'optimistic', 'nominal'],
                        default='nominal',
                        help='what assumption-level should the simulations use')
    parser.add_argument('-n', '--notest', action='store_true', help='turn off asymptomatic testing')
    parser.add_argument('-o', '--outputdir', default=BASE_DIRECTORY, 
                        help='directory to store simulation output')
    parser.add_argument('-s', '--silent', action='store_true', help='turn off script output')
    parser.add_argument('config',
            help='YAML config file specifying which parameters to vary across the different processes')

    args = parser.parse_args()

    sim_config = yaml.load(open(args.config))

    params = ParamConfig.load_config(args.when, not args.notest, args.assumption)


    timestamp = time.time()
    sim_id = "{timestamp}-{when}-{assumption}-{testing}-{param_to_vary}".format(
                timestamp=str(time.time()).split('.')[0], 
                when=args.when, 
                assumption=args.assumption,
                testing= "notest" if args.notest else "withtest",
                param_to_vary = sim_config['param_to_vary'])
    if not args.silent:
        print("Simulation ID: {}".format(sim_id))

    basedir = args.outputdir

    if not os.path.isdir(basedir):
        print("Directory {} does not exist. Please create it.".format(basedir))
        exit()

    sim_main_dir = basedir + "/" + str(sim_id)
    try:
        os.mkdir(sim_main_dir)
        if not args.silent:
            print("Output directory {} created".format(sim_main_dir))
    except FileExistsError:
        print("Output directory {} already exists".format(sim_main_dir))
        exit()
    except FileNotFoundError:
        print("Directory {} cannot be created.".format(sim_main_dir))
        exit()


    if 'base_params_to_update' in sim_config and sim_config['base_params_to_update'] != None:
        for param, val in sim_config['base_params_to_update'].items():
            update_params(params, param, val)
    
    param_to_vary = sim_config['param_to_vary']
    param_values = sim_config['parameter_values']
    if 'ntrajectories' in sim_config:
        ntrajectories = sim_config['ntrajectories']
    else:
        ntrajectories = 500

    if 'time_horizon' in sim_config:
        time_horizon = sim_config['time_horizon']
    else:
        time_horizon = 112 # 16 weeks

    if len(param_values) == 0:
        print("Empty list of parameters given; nothing to do")
        exit()

    for param_val in param_values:
        # create the relevant subdirectory
        sim_sub_dir = "{}/{}-{}".format(sim_main_dir, param_to_vary, param_val)
        os.mkdir(sim_sub_dir)
        if not args.silent: 
            print("Created directory {} to save output".format(sim_sub_dir))
        # instantiate relevant sim params
        sim_params = params.copy()

        
        update_params(sim_params, param_to_vary, param_val)        
        
        dill.dump(sim_params, open("{}/sim_params.dill".format(sim_sub_dir), "wb"))
        if not args.silent:
            print("Saved sim_params to dill file")
        # start new process
        fn_args = (sim_sub_dir, sim_params, ntrajectories, time_horizon)
        proc = multiprocessing.Process(target = run_background_sim, args=fn_args)
        #proc.daemon = True
        proc.start()
        if not args.silent:
            print("starting process for {} value {}".format(param_to_vary, param_val))
            print("process PID = {}".format(proc.pid))

    if not args.silent:
        print("Waiting for processes to finish...")

        








