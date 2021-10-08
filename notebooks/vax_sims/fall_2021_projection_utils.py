import numpy as np
from scipy.stats import poisson
from scipy.stats import norm


ACTUAL_TRAJ = [2,2,11,23,42,57,55,39,31,48,42,25,29,4,9,8,13,14,8,8,1,7,8,5,3,6,3,0,3,4,9,3,8,1]#,0,0]


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


"""
Code for transforming trajectory data & computing parameter likelihoods based on the Log-Normal assumption
"""
import numpy as np

def convert_cum_traj_to_daily_count(traj):
    prev_cum_cases = 0
    daily_count = []
    for current_cum_cases in traj:
        daily_count.append(current_cum_cases - prev_cum_cases)
        prev_cum_cases = current_cum_cases
    return daily_count


def convert_daily_count_traj_to_weekly_count(daily_count_traj):
    weekly_count = []
    lower_idx = 0
    upper_idx = 7 
    while lower_idx <= len(daily_count_traj):
        weekly_count.append(sum(daily_count_traj[lower_idx:upper_idx]))
        lower_idx = upper_idx
        upper_idx += 7
    assert(np.abs(sum(weekly_count) - sum(daily_count_traj)) < 1e-5)
    return weekly_count


def get_weekly_counts(cum_traj):
    daily_count = convert_cum_traj_to_daily_count(cum_traj)
    return convert_daily_count_traj_to_weekly_count(daily_count)


def compute_lognormal_loglik(cum_trajs):
    actual_weekly_traj = convert_daily_count_traj_to_weekly_count(ACTUAL_TRAJ)
    num_weeks = len(actual_weekly_traj)
    
    # estimate lognormal parameters
    log_counts_by_week = {week_idx:[] for week_idx in range(num_weeks)}
    for cum_traj in cum_trajs:
        weekly_counts = get_weekly_counts(cum_traj)
        for week_idx in range(num_weeks):
            log_counts_by_week[week_idx].append(np.log(weekly_counts[week_idx] + 1e-10))
    
    means_by_week = {week_idx: np.mean(log_counts_by_week[week_idx]) for week_idx in range(num_weeks)}
    stddevs_by_week = {week_idx: np.std(log_counts_by_week[week_idx]) for week_idx in range(num_weeks)}
    
    # compute log likelihood for observed data
    loglik = 0
    for week_idx, actual_positives in enumerate(actual_weekly_traj):
        log_positives = np.log(actual_positives)
        #scaled_log_positives = (log_positives - means_by_week[week_idx]) / stddevs_by_week[week_idx]
        #print(log_positives, means_by_week[week_idx], stddevs_by_week[week_idx], scaled_log_positives)
        loglik += norm.logpdf(log_positives, loc=means_by_week[week_idx], 
                                  scale=stddevs_by_week[week_idx])
    
    return loglik


"""
Code for transforming trajectory data & computing parameter likelihoods based on the Poisson assumption
"""

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


def get_simulated_positives_by_day(simulated_cumulative_trajs):
    simulated_mean_traj = np.mean(np.array(simulated_cumulative_trajs), axis=0)
    
    true_positives_by_day = ACTUAL_TRAJ
    num_days = min(len(true_positives_by_day), len(simulated_mean_traj))
    
    simulated_positives_by_day = [simulated_mean_traj[0]]
    day = 1
    while day < num_days:
        simulated_positives_by_day.append(simulated_mean_traj[day] - simulated_mean_traj[day - 1])
        day += 1
    return simulated_positives_by_day

"""

def get_positives_by_week_from_traj(simulated_positives_by_day):
    simulated_cum_positives_by_week = []
    lower_idx = 0
    upper_idx = 7 
    while lower_idx <= len(simulated_positives_by_day):
        simulated_cum_positives_by_week.append(sum(simulated_positives_by_day[0:upper_idx]))
        lower_idx = upper_idx
        upper_idx += 7
    #assert(np.abs(sum(simulated_cum_positives_by_week) - sum(simulated_positives_by_day)) < 1e-5)

    simulated_positives_by_week = [simulated_cum_positives_by_week[0]]
    week = 1
    while week < len(simulated_cum_positives_by_week):
        simulated_positives_by_week.append(simulated_cum_positives_by_week[week] - simulated_cum_positives_by_week[week - 1])
        week += 1
    return simulated_positives_by_week
    """



def get_positives_by_week(simulated_cumulative_trajs):
    simulated_positives_by_day = get_simulated_positives_by_day(simulated_cumulative_trajs)
    simulated_positives_by_week = []
    lower_idx = 0
    upper_idx = 7 
    while lower_idx <= len(simulated_positives_by_day):
        simulated_positives_by_week.append(sum(simulated_positives_by_day[lower_idx:upper_idx]))
        lower_idx = upper_idx
        upper_idx += 7
    assert(np.abs(sum(simulated_positives_by_week) - sum(simulated_positives_by_day)) < 1e-5)
    return simulated_positives_by_week


def get_true_positives_by_week():
    positives_by_day = ACTUAL_TRAJ
    positives_by_week = []
    lower_idx = 0
    upper_idx = 7 
    while lower_idx <= len(positives_by_day):
        positives_by_week.append(sum(positives_by_day[lower_idx:upper_idx]))
        lower_idx = upper_idx
        upper_idx += 7
    assert(sum(positives_by_week) == sum(positives_by_day))
    return positives_by_week

def compute_log_likelihood_by_week(simulated_cumulative_trajs, eps=1e-5):
    simulated_positives_by_week = get_positives_by_week(simulated_cumulative_trajs)
    true_positives_by_week = get_true_positives_by_week()
    loglik = 0
    for true_positives, simulated_positives in zip(true_positives_by_week, simulated_positives_by_week):
        loglik += np.log(poisson.pmf(true_positives, simulated_positives) + eps)
    
    return loglik


def aggregate_trajs_student_only(inf_trajs_by_group):
    aggregated_trajs = []
    for trajs in inf_trajs_by_group:
        num_groups = len(trajs) - 2 # subtract off the two employee groups at the end
        assert(num_groups >= 1)
        aggregated_traj = trajs[0]
        group_idx = 1
        while group_idx < num_groups:
            aggregated_traj += trajs[group_idx]
            group_idx += 1
        aggregated_trajs.append(aggregated_traj)
    return aggregated_trajs


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

