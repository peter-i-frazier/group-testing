import sys
import os
module_path = os.path.abspath(os.path.join('../..'))
sys.path.append(module_path)
import numpy as np
import pandas as pd

from multi_group_simulation import MultiGroupSimulation
from load_params import load_params
from analysis_helpers import binomial_exit_function

NUM_GROUPS=4
PARAMS_FOLDER = "/home/aaj54/group-testing/src/simulations_v2/params/spring_2021_calibration/"

def load_actual_trajectories():
    actual_counts_by_group = \
            {group_idx: pd.read_csv(PARAMS_FOLDER + "actual_counts_group_{}.csv".format(group_idx)) for group_idx in [1,2,3,4]}
    cumulative_actual_counts = actual_counts_by_group[1][['cum_case_count']].copy()
    for idx in [2,3,4]:
        cumulative_actual_counts = cumulative_actual_counts.add(actual_counts_by_group[idx][['cum_case_count']])
    return cumulative_actual_counts, actual_counts_by_group


def compute_score(sim_df, actual_df):
    # return MSE of the day/day error between sim_df and actual_df
    errors = np.array(sim_df['QI'] - actual_df['cum_case_count'])
    squared_errs = [error ** 2 for error in errors]
    return np.mean(squared_errs)

def compute_avg_score(sim_traj, actual):
    errors = np.array(np.array(sim_traj) - np.array(actual))
    squared_errs = [error ** 2 for error in errors]
    return np.mean(squared_errs)

class SpringCalibration:
    # 125 days from january 21 through may 25
    # day 78 corresponds to april 9, which is the week we saw increased testing frequency according to the data
    # day 65 corresponds to march 27
    def __init__(self, exposed_infection_p, sim_length=124, change_t=65, initial_prev_days=18):

        # transpose of 4x4 matrix in pnas calibration doc
        # needs to be updated
        interaction_matrix = np.loadtxt(PARAMS_FOLDER + "/interaction_matrix.csv")
        # interaction_matrix = np.array([[0.9464, 0.0943, 0.    , 0.    ],
        #                                [0.2036, 0.3205, 0.    , 0.0189],
        #                                [0.    , 0.    , 0.5156, 0.0094],
        #                                [0.0236, 0.0044, 0.0273, 0.0849]])

        # need to double check
        self.empirical_test_frequencies = [0.373, 0.283, 0.141, 0.141]
        # need to double check
        self.post_apr_9_mba_test_freq = 0.237
        self.change_t = change_t
        self.exposed_infection_p = exposed_infection_p
        self.sim_length = sim_length
        self.initial_prev_days = initial_prev_days

        self.cumulative_actual_counts, self.actual_counts_by_group = load_actual_trajectories()

        self.original_initial_prev = list()
        self.original_outside_prev = list()

        group_names = []
        group_params = []
        for idx in range(NUM_GROUPS):
            name, params = load_params(PARAMS_FOLDER + "/group_{}_students_spring_2021.yaml".format(idx+1))
            #name, params = load_params("./group_{}_students_spring_2021.yaml".format(idx+1))
            params['test_population_fraction'] = self.empirical_test_frequencies[idx]
            params['exposed_infection_p'] = self.exposed_infection_p
            params['sample_QI_exit_function'] = binomial_exit_function(0)
            self.original_initial_prev.append(params['initial_ID_prevalence'])
            params['initial_ID_prevalence'] = params['initial_ID_prevalence'] / initial_prev_days
            self.original_outside_prev.append(params['daily_outside_infection_p'])
            params['daily_outside_infection_p'] += params['initial_ID_prevalence'] / initial_prev_days
            group_names.append(name)
            group_params.append(params)


        self.multigroup_sim = MultiGroupSimulation(group_params, interaction_matrix, group_names)


    def reset_initial_state(self):
        for idx in range(NUM_GROUPS):
            self.multigroup_sim.sims[idx].test_pop_fraction = self.empirical_test_frequencies[idx]

        self.multigroup_sim.reset_initial_state()


    def increase_mba_test_freq(self):
        self.multigroup_sim.sims[2].test_pop_fraction = self.post_apr_9_mba_test_freq


    def update_outside_inf_rate(self):
        for idx, group in enumerate(self.multigroup_sim.sims):
            group.daily_outside_infection_p = self.original_outside_prev[idx]

    def run_new_trajectory(self):
        T = self.sim_length
        self.reset_initial_state()
        for t in range(T):
            if t == self.change_t:
                self.increase_mba_test_freq()
            if t == self.initial_prev_days:
                self.update_outside_inf_rate()
            self.multigroup_sim.step()

        for sim in self.multigroup_sim.sims:
            sim.update_severity_levels()

        sim_df = self.multigroup_sim.sims[0].sim_df
        for sim in self.multigroup_sim.sims[1:]:
            sim_df = sim_df.add(sim.sim_df)
        return sim_df, [sim.sim_df for sim in self.multigroup_sim.sims]



#    def score_trajectory(self, cumulative_sim_df, sim_dfs_by_group):
#        cumulative_score = compute_score(cumulative_sim_df, self.cumulative_actual_counts)
#        scores_by_group = [compute_score(sim_dfs_by_group[group_idx-1], self.actual_counts_by_group[group_idx])
#                                for group_idx in [1,2,3,4]]
#        return cumulative_score, np.mean(scores_by_group)


    def score_avg_trajectories(self, cumulative_traj, avg_trajs):
        cumulative_score = compute_avg_score(cumulative_traj, self.cumulative_actual_counts)
        scores_by_group = [compute_avg_score(avg_trajs[group_idx-1], self.actual_counts_by_group[group_idx]['cum_case_count']) for group_idx in [1,2,3,4]]
        return cumulative_score, scores_by_group


    def run_and_score_trajectories(self, ntrajs):
        cumulative_scores = []
        avg_group_scores = []
        sim_outputs = list()
        cum_trajs = list()
        trajs_by_group = list()
        for _ in range(4):
            trajs_by_group.append(list())
        for _ in range(ntrajs):
            cum_df, dfs_by_group = self.run_new_trajectory()
            sim_outputs.append((cum_df, dfs_by_group))
            cum_trajs.append(cum_df['QI'])
            for idx, df in enumerate(dfs_by_group):
                trajs_by_group[idx].append(df['QI'])
            # cum_score, avg_group_score = self.score_trajectory(cum_df, dfs_by_group)
            # cumulative_scores.append(cum_score)
            # avg_group_scores.append(avg_group_score)
        avg_cum_traj = np.mean(cum_trajs, axis=0)
        avg_trajs = list()
        for idx, trajs in enumerate(trajs_by_group):
            avg_trajs.append(np.mean(trajs, axis=0))
        return (self.score_avg_trajectories(avg_cum_traj, avg_trajs), sim_outputs)

