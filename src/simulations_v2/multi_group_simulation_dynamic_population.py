"""
implements a Multi-Group Population-Level stochastic simulation
makes the following assumptions:
    * inter-group and intra-group interactions are perfectly mixed, but the
        rates may vary
    * no inter-group contact tracing
"""
from stochastic_simulation import StochasticSimulation
import numpy as np

class MultiGroupSimulationDynamicPopulation:
    def __init__(self, sims, # a list of already-initialized DynamicPopulationSims
                            # in the movein period, assumes that inter-group contacts only happen
                            # for the subset of the population not in self-isolation
                        interaction_matrix, # contact matrix
                        group_names=[]):
        
        self.sims = sims
        self.N = len(sims)
        self.interaction_matrix = interaction_matrix

    def get_interaction_mtx(self):
        return self.interaction_matrix
    

    def get_total_population(self):
        return sum([sim.pop_size for sim in self.sims])


    def set_interaction_mtx(self, interaction_mtx):
        self.interaction_matrix = interaction_mtx

    def run_new_trajectory(self, T):
        for t in range(T):
            self.step()

    def get_free_total(self, i):
        sim_obj_i = self.sims[i].get_sim_obj()
        # get the free-total count from group i
        free_total = self.get_free_infectious(i)

        if sim_obj_i.pre_ID_state == 'detectable':
            free_total += sum(sim_obj_i.pre_ID)

        free_total += sim_obj_i.S + sim_obj_i.R + sum(sim_obj_i.E)
        return free_total


    def get_free_infectious(self, i):
        # get the free-infectious total from group j
        sim_obj_i = self.sims[i].get_sim_obj()
        if sim_obj_i.pre_ID_state == 'infectious':
            free_infectious = sum(sim_obj_i.pre_ID)
        else:
            free_infectious = 0

        free_infectious += sum(sim_obj_i.ID)
        free_infectious += sum(sim_obj_i.SyID_mild)
        free_infectious += sum(sim_obj_i.SyID_severe)

        return free_infectious


    def step(self):
        # do inter-group interactions first, so that no updates happen after each sim adds
        # a row to their dataframe
        new_E_holder = [0]* self.N
        for i in range(self.N):
            sim_obj_i = self.sims[i].get_sim_obj()
            for j in range(self.N):
                if i == j:
                    continue
                
                free_susceptible_i = sim_obj_i.S

                interactions_lambda_i_j = self.interaction_matrix[i,j]

                free_infectious_j = self.get_free_infectious(j)
                free_total_j = self.get_free_total(j)
                
                if free_total_j == 0:
                    continue

                poisson_param = free_susceptible_i * interactions_lambda_i_j * \
                                    free_infectious_j / free_total_j

                n_susceptible_infectious_contacts = np.random.poisson(poisson_param)

                new_E = np.random.binomial(n_susceptible_infectious_contacts, sim_obj_i.exposed_infection_p)
                new_E_holder[i] += new_E

        for i in range(self.N):
            sim_obj_i = self.sims[i].get_sim_obj()
            sim_obj_i.add_new_infections(new_E_holder[i])

        # do individual-group steps
        for sim in self.sims:
            sim.step()
