"""
implements a Multi-Group Population-Level stochastic simulation
makes the following assumptions:
    * inter-group and intra-group interactions are perfectly mixed, but the
        rates may vary
    * no inter-group contact tracing
"""
from stochastic_simulation import StochasticSimulation
import numpy as np

class MultiGroupSimulation:
    def __init__(self, group_params,
                        interaction_matrix,
                        group_names=[]):
        """
        group_params: A list of dictionaries of length N.  Each dictionary is
                    used as the config for a different individual-group
                    StochasticSimulation object.
        interaction_matrix: A N x N matrix such that the (i,j) element indicates
                    the rate at which members of group i are exposed to members of
                    group j.
                    Specifically, each free member of group i has Poisson(lambda[i,j])
                    contacts each day with a free member of group j
        group_names: An optional list of strings of length N indicating the name
                    of each group
        """
        self.sims = [StochasticSimulation(params) for params in group_params]
        self.N = len(group_params)
        self.interaction_matrix = interaction_matrix

        self.original_interaction_matrix = interaction_matrix
        self.original_daily_contacts = [sim.daily_contacts_lambda for sim in self.sims]
        self.lockdown_in_effect = False
        self.simulate_lockdown = False


    def configure_lockdown(self, 
                            post_lockdown_interaction_matrix,
                            new_cases_threshold, # observed new cases found by testing / self-reporting that trigger lockdown
                                                # specified as proportion of total population (i.e. the value should be between 0 and 1)
                            new_cases_time_window # number of days over which new cases are computed for the previous threshold
                            ):
        self.simulate_lockdown = True
        self.post_lockdown_interaction_matrix = post_lockdown_interaction_matrix
        self.new_cases_threshold = new_cases_threshold
        self.new_cases_time_window = new_cases_time_window
        self.new_case_counts = [0] * new_cases_time_window


    def step_lockdown_status(self):
        assert(self.simulate_lockdown)
        self.update_case_counts()
        if sum(self.new_case_counts) >= self.new_cases_threshold * self.get_total_population():
            self.lockdown_in_effect = True
            self.interaction_matrix = self.post_lockdown_interaction_matrix
            for i in range(self.N):
                self.sims[i].daily_contacts_lambda = self.post_lockdown_interaction_matrix[i,i]


    def update_case_counts(self):
        new_cases_today = 0
        for sim in self.sims:
            new_cases_today += sim.new_QS_from_last_test
            new_cases_today += sim.new_QI_from_last_test
            new_cases_today += sim.new_QI_from_self_reports

        #shift case count array down
        self.new_case_counts.pop(0)
        self.new_case_counts.append(new_cases_today)


    def get_interaction_mtx(self):
        return self.interaction_matrix
    

    def get_total_population(self):
        return sum([sim.pop_size for sim in self.sims])


    def set_interaction_mtx(self, interaction_mtx):
        self.interaction_matrix = interaction_mtx


    def reset_initial_state(self):
        self.lockdown_in_effect = False
        self.interaction_matrix = self.original_interaction_matrix

        if self.simulate_lockdown:
            self.new_case_counts = [0] * self.new_cases_time_window
            for sim, contacts in zip(self.sims, self.original_daily_contacts):
                sim.daily_contacts_lambda = contacts


        for sim in self.sims:
            sim.reset_initial_state()
    

    def run_new_trajectory(self, T):
        self.reset_initial_state()
        lockdown_statuses = []
        for _ in range(T):
            self.step()
            if self.simulate_lockdown:
                self.step_lockdown_status()
            lockdown_statuses.append(self.lockdown_in_effect)

        sim_df = self.sims[0].sim_df
        for sim in self.sims[1:]:
            sim_df = sim_df.add(sim.sim_df)
        return lockdown_statuses, sim_df


    def get_free_total(self, i):
        # get the free-total count from group i
        free_total = self.get_free_infectious(i)

        if self.sims[i].pre_ID_state == 'detectable':
            free_total += sum(self.sims[i].pre_ID)

        free_total += self.sims[i].S + self.sims[i].R + sum(self.sims[i].E)
        return free_total


    def get_free_infectious(self, i):
        # get the free-infectious total from group j

        if self.sims[i].pre_ID_state == 'infectious':
            free_infectious = sum(self.sims[i].pre_ID)
        else:
            free_infectious = 0

        free_infectious += sum(self.sims[i].ID)
        free_infectious += sum(self.sims[i].SyID_mild)
        free_infectious += sum(self.sims[i].SyID_severe)

        return free_infectious


    def get_quarantine_susceptible(self, i):
        return self.sims[i].QS


    def get_quarantine_infected(self, i):
        return self.sims[i].QI


    def step(self):
        # do inter-group interactions first, so that no updates happen after each sim adds
        # a row to their dataframe
        new_E_holder = [0]* self.N
        for i in range(self.N):
            for j in range(self.N):
                if i == j:
                    continue
                free_susceptible_i = self.sims[i].S

                interactions_lambda_i_j = self.interaction_matrix[i,j]

                free_infectious_j = self.get_free_infectious(j)
                free_total_j = self.get_free_total(j)

                poisson_param = free_susceptible_i * interactions_lambda_i_j * \
                                    free_infectious_j / free_total_j

                n_susceptible_infectious_contacts = np.random.poisson(poisson_param)

                new_E = np.random.binomial(n_susceptible_infectious_contacts, self.sims[i].exposed_infection_p)
                new_E_holder[i] += new_E

        for i in range(self.N):
            self.sims[i].add_new_infections(new_E_holder[i])

        # do individual-group steps
        for sim in self.sims:
            sim.step()
