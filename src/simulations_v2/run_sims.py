import sys
import numpy as np
import yaml
import time
import os
import multiprocessing
from analysis_helpers import run_multiple_trajectories
import dill

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


def update_params(sim_params, param_to_vary, param_val):
    # VERY TEMPORARY HACK TO GET SENSITIVITY SIMS FOR ASYMPTOMATIC %
    if param_to_vary == 'asymptomatic_p':
        assert(sim_params['mild_severity_levels'] == 1)
        curr_prevalence_dist = params['severity_prevalence']
        assert(param_val >= 0 and param_val <= 1)
        curr_prevalence_dist[0] = param_val
        remaining_mass = sum(curr_prevalence_dist[1:])

        # need to scale so that param_val + x * remaning_mass == 1
        scale = (1 - param_val) / remaining_mass
        idx = 1
        while idx < len(curr_prevalence_dist):
            curr_prevalence_dist[idx] = curr_prevalence_dist[idx] * scale
        assert(np.isclose(curr_prevalence_dist, 1))
        sim_params['severity_prevalence'] = curr_prevalence_dist

    # VERY TEMPORARY HACK TO GET SENSITIVITY SIMS WORKING FOR CONTACT RECALL %
    elif param_to_vary == 'contact_tracing_constant':
        num_isolations = sim_params['cases_isolated_per_contact']
        base_recall = sim_params['contact_recall']
        new_isolations = num_isolations * param_val / base_recall
        new_quarantines = max(7 - new_isolations, 0)
        sim_params['cases_isolated_per_contact'] = new_isolations
        sim_params['cases_quarantined_per_contact'] = new_quarantines
    elif param_to_vary == 'contact_tracing_isolations':
        sim_param['cases_isolated_per_contact'] = param_val
        sim_param['cases_quarantined_per_contact'] = 7 - param_val
    else:
        sim_params[param_to_vary] = param_val


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage: python {} yaml-config-file 'fall'|'base'|'fall_old_severity' (optional:simulation-name) (optional:basedir)".format(sys.argv[0]))
        exit()

    # TODO: should we change this to a better naming convention?
    #  e.g., fall, fall-optimistic, fall-pessimistic, summer, etc.
    if sys.argv[2] == 'fall':
        from fall_params import base_params
    elif sys.argv[2] == 'base':
        from base_params import base_params
    elif sys.argv[2] == 'fall_old_severity':
        from fall_params import base_params_tuesday_severity as base_params
    elif sys.argv[2] == 'fall_slightly_pessimistic_testing':
        from fall_slightly_pessimistic import base_params_testing as base_params
    elif sys.argv[2] == 'fall_slightly_optimistic_testing':
        from fall_slightly_optimistic import base_params_testing as base_params
    elif sys.argv[2] == 'fall_realistic_testing':
        from fall_realistic import base_params_testing as base_params
    elif sys.argv[2] == 'fall_slightly_pessimistic':
        from fall_slightly_pessimistic import base_params as base_params
    elif sys.argv[2] == 'fall_slightly_optimistic':
        from fall_slightly_optimistic import base_params as base_params
    elif sys.argv[2] == 'fall_realistic':
        from fall_realistic import base_params as base_params
    elif sys.argv[2] == 'june_slightly_pessimistic_testing':
        from june_slightly_pessimistic import base_params_testing as base_params
    elif sys.argv[2] == 'june_slightly_optimistic_testing':
        from june_slightly_optimistic import base_params_testing as base_params
    elif sys.argv[2] == 'june_realistic_testing':
        from june_realistic import base_params_testing as base_params
    elif sys.argv[2] == 'june_slightly_pessimistic':
        from june_slightly_pessimistic import base_params as base_params
    elif sys.argv[2] == 'june_slightly_optimistic':
        from june_slightly_optimistic import base_params as base_params
    elif sys.argv[2] == 'june_realistic':
        from june_realistic import base_params as base_params
    else:
        print("Error: second argument must be 'fall' or 'base' or 'fall_old_severity', but got {}".format(sys.argv[2]))


    sim_config = yaml.load(open(sys.argv[1]))


    params = base_params.copy()
    if 'base_params_to_update' in sim_config and sim_config['base_params_to_update'] != None:
        for param, val in sim_config['base_params_to_update'].items():

            # got rid of the following sanity check because we are adding a lot of parameters on the fly
            # should eventually have some kind of sanity check here when the codebase stabilizes

            #if param not in params and param != 'contact_trace_testing_frac' and param != 'init_ID_prevalence_stochastic':
            #    print("Configuration attempting to modify non-existent parameter {}".format(
            #        param))
            #    exit()
    
            update_params(params, param, val)

    sim_timestamp = time.time()
    if len(sys.argv) > 3:
        sim_id = "{}.{}".format(sys.argv[3], sim_timestamp)
    else:
        sim_id = str(sim_timestamp)
    print("Simulation ID: {}".format(sim_id))


    if len(sys.argv) > 4:
        basedir = sys.argv[4]
    else:
        basedir = BASE_DIRECTORY

    if not os.path.isdir(basedir):
        print("Directory {} does not exist. Please create it.".format(basedir))
        exit()

    sim_main_dir = basedir + str(sim_id)
    try:
        os.mkdir(sim_main_dir)
        print("Output directory {} created".format(sim_main_dir))
    except FileExistsError:
        print("Output directory {} already exists".format(sim_main_dir))
        exit()
    except FileNotFoundError:
        print("Directory {} cannot be created.".format(sim_main_dir))
        exit()


    
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
        sim_sub_dir = "{}/{}.{}".format(sim_main_dir, param_to_vary, param_val)
        os.mkdir(sim_sub_dir)
        print("Created directory {} to save output".format(sim_sub_dir))
        # instantiate relevant sim params
        sim_params = params.copy()

        
        update_params(sim_params, param_to_vary, param_val)        
        
        dill.dump(sim_params, open("{}/sim_params.dill".format(sim_sub_dir), "wb"))
        print("Saved sim_params to dill file")
        # start new process
        fn_args = (sim_sub_dir, sim_params, ntrajectories, time_horizon)
        proc = multiprocessing.Process(target = run_background_sim, args=fn_args)
        #proc.daemon = True
        proc.start()
        print("starting process for {} value {}".format(param_to_vary, param_val))
        print("process PID = {}".format(proc.pid))

    print("Waiting for processes to finish...")

        








