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

    def reset_initial_state(self):
        for sim in self.sims:
            sim.reset_initial_state()

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


























