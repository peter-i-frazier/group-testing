import numpy as np
from scipy.stats import poisson

ACTUAL_TRAJ = [2,2,11,23,42,57,55,39,31,48,42,25,29,4,9,8,13,14,8,8,1,7,8,5,3,6,3,0,3,4,9,3,8,1,0,0]


UNCERTAINTY_PARAMS = ['vax_susc_mult', 'vax_transmission_mult', 'contacts_per_day_mult', 'outside_infection_rate_mult',
                      'cases_isolated_per_contact_trace', 'initial_ID_prevalence']

UNCERTAINTY_PARAM_RANGES = {
    'vax_susc_mult': (0.097608, 0.941192), # 0.5194 +/- 1.96 * 0.2152
    'vax_transmission_mult': (0.25, 1),
    'contacts_per_day_mult': (0.9,2.7),
    'outside_infection_rate_mult': (1, 5),
    'cases_isolated_per_contact_trace': (0.5,1.5),
    'initial_ID_prevalence': (0.003, 0.0054)
}

def compute_log_likelihood(simulated_cumulative_trajs, eps=1e-5):
    true_positives_by_day = ACTUAL_TRAJ
    simulated_mean_traj = np.mean(np.array(simulated_cumulative_trajs), axis=0)
    
    num_days = min(len(true_positives_by_day), len(simulated_mean_traj))
    
    simulated_positives_by_day = [simulated_mean_traj[0]]
    day = 1
    while day < num_days:
        simulated_positives_by_day.append(simulated_mean_traj[day] - simulated_mean_traj[day - 1])
        day += 1
    
    loglik = 0
    for true_positives, simulated_positives in zip(true_positives_by_day, simulated_positives_by_day):
        loglik += np.log(poisson.pmf(true_positives, simulated_positives) + eps)
    
    return loglik


def aggregate_trajs(inf_trajs_by_group):
    aggregated_trajs = []
    for trajs in inf_trajs_by_group:
        num_groups = len(trajs)
        assert(num_groups >= 1)
        aggregated_traj = trajs[0]
        group_idx = 1
        while group_idx < num_groups:
            aggregated_traj += trajs[group_idx]
            group_idx += 1
        aggregated_trajs.append(aggregated_traj)
    return aggregated_trajs

