from base_population import BasePopulation
import numpy as np
import random

"""
Implements the population dynamics for the individual-interaction population model
which models stochastic day-by-day interactions between individuals in the population

In addition to the methods and attributes inherited from BasePopulation,
the IndividualInteractionPopulation class has the following attributes:
    * interaction_lambda: parameter to a Poisson() distribution
        governing the total number of interactions that occur each day.
        The distribution works as follows:
        - A number N is sampled from Poisson(lambda) each day
        - N pairs of agent IDs (i,j) are sampled uniformly at random (with replacement)
        - Each (i,j) is an interaction which occurs on the specified day
        - If i = j then the interaction has no effect
        - If i & j share the same infection status then the interaction has no effect
        - If one of (i,j) is infected and the other is uninfected, then the disease is
          passed along with some probability, governed by the following parameter.
    * interaction_infection_pct: probability that the infection is passed along
                                for an interaction (i,j) where one agent is infected
                                and the other is uninfected
    * days_since_last_interaction: a matrix of size n_agents x n_agents such that the (i,j)
                                    component specifies how many days it has been since 
                                    (i,j) had an interaction
"""

class IndividualInteractionPopulation(BasePopulation):


    def __init__(self, n_agents,
                        disease_length,
                        quarantine_length,
                        days_until_symptomatic,
                        interaction_frequency_lambda,
                        interaction_infection_pct,
                        initial_days_since_interaction=300,
                        initial_prevalence=0):

        self.interaction_lambda = interaction_frequency_lambda
        self.interaction_infection_pct = interaction_infection_pct
        self.days_since_last_interaction = np.full((n_agents, n_agents), 
                                                initial_days_since_interaction)

        super().__init__(n_agents, 
                disease_length, 
                quarantine_length,
                days_until_symptomatic,
                initial_prevalence)


    def get_summary_data(self):
        return {'num_quarantined':self.get_num_quarantined(),
                'num_infected': self.get_num_infected()}
    
    def respond_to_test_results(self, test_results, day_of_test):
        pass

    def resolve_infection(self, agent_id):
        # Do nothing when an infection is resolved -- the agent_id goes back
        # to being a normal agent, as if they were never infected
        pass
    

    def resolve_quarantine(self, agent_id):
        # Do nothing when an agent leaves quarantine
        pass


    def step(self):
        # Implement the infection dynamics for the stochastic interactions which occur
        # each day

        scaled_lambda = self.interaction_lambda / float(self.n_agents ** 2)

        n_agents_free = self.n_agents - self.get_num_quarantined()

        # sample interaction counts which occur between unquarantined ("free") agents
        interaction_matrix = np.random.poisson(scaled_lambda, (n_agents_free, n_agents_free))


        self.days_since_last_interaction = self.days_since_last_interaction + \
                                            np.ones(self.days_since_last_interaction.shape)
        
        # i_free_idx and j_free_idx are indices over unquarantined agents
        # for use with interaction_matrix
        # i and j correspond directly to agent_ids

        i_free_idx = 0
        for i in self.agents:

            if self.quarantine_status[i]:
                continue

            j_free_idx = 0

            for j in self.agents:

                # ensure we loop over each pair (i,j) only once, by only considering j <= i
                if j > i:
                    continue

                if self.quarantine_status[j]:
                    continue
                
                n_interactions = interaction_matrix[i_free_idx, j_free_idx] + interaction_matrix[j_free_idx, i_free_idx]

                # update interaction count if i&j interact today
                if n_interactions > 0:
                    self.days_since_last_interaction[i,j] = 0
                    self.days_since_last_interaction[j,i] = 0

                # sample probability of passing along disease only if exactly one
                # individual is infected
                if self.infection_status[i] + self.infection_status[j] == 1:
                    for _ in range(n_interactions):
                        if random.random() < self.interaction_infection_pct:
                            self.infect_agent(i)
                            self.infect_agent(j)

                j_free_idx += 1

            
            i_free_idx += 1

        assert(i_free_idx == n_agents_free)

        






































