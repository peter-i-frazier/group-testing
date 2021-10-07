import dill as dill
import time
import pandas as pd
from os import listdir
from os.path import isfile, join
from fall_2021_projection_utils import aggregate_trajs, \
                                       compute_lognormal_loglik
from fall_2021_projection_utils import UNCERTAINTY_PARAMS
import numpy as np

max_files = 10

if __name__ == "__main__":
    dill_output_paths = ["fall_2021_prior_samples:1633381945", "fall_2021_prior_samples:1633456006"]

    files = []
    for path in dill_output_paths:
        files_in_path = [join(path, f) for f in listdir(path) if isfile(join(path, f)) and 'with_trajectories' in f ]
        files = files + files_in_path

    sampled_points = {}
    aggregated_trajs = {}


    print("preparing to load {} files out of a total of {} available".format(max_files, len(files)))

    if max_files < len(files):
        files = files[0:max_files]

    count = 0
    for f in files:
        with open(f, "rb") as fhandle:
            [point, inf_trajs_by_group] = dill.load(fhandle)
            sampled_points[f] = point
            aggregated_trajs[f] = aggregate_trajs(inf_trajs_by_group)
        count += 1
        if count % 10 == 0:
            print("Loaded {} points".format(count))
            #break

    print("done loading {} files".format(len(files)))
    
    logliks = {}
    for f in sampled_points:
        logliks[f] = compute_lognormal_loglik(aggregated_trajs[f])


    param_vals = {}
    param_logliks = {}

    for idx, param in enumerate(UNCERTAINTY_PARAMS):
        param_vals[param] = []
        param_logliks[param] = []
        for f in sampled_points:
            param_vals[param].append(sampled_points[f][idx])
            param_logliks[param].append(logliks[f])

    df = pd.DataFrame(param_vals)
    df['log_likelihood'] = param_logliks[UNCERTAINTY_PARAMS[0]]
    df['combined_spread_mult'] = df['vax_transmission_mult'] * df['vax_susc_mult'] * \
					  df['contacts_per_day_mult']
    min_loglik = min(df[df['log_likelihood'] > -np.inf]['log_likelihood'])
    df.replace(-np.inf, min_loglik, inplace=True)                                        

    likelihoods = [np.exp(loglik) for loglik in param_logliks[UNCERTAINTY_PARAMS[0]]]
    df['likelihood'] = likelihoods
    normalizer = df['likelihood'].sum()
    df['posterior'] = df['likelihood'] / normalizer

    csv_name = time.strftime("posterior_csvs/%y_%m_%d_%H:%M_posteriors.csv", time.localtime())
    df.to_csv(csv_name)
    print("Saved CSV file {}".format(csv_name))
