"""
Implements the agent-based version of the simulation described
here: 
https://docs.google.com/document/d/17tm2uPOUxrcfXG-W2paw02K6pSx_iDMXY51n5RSjTMI/edit#
"""
"""
TODO:
    * implement self-reporting from symptomatic individuals
    * implement adaptive testing
"""

import numpy as np

from agent import Agent
from surveillance_testing import SurveillanceTesting
 
class MultiAgentSim:

    def __init__(self, 
            n_agents, 
            init_infection_p,
            mean_test_delay=1,
            mean_contact_trace_delay=1,
            contact_trace_time_window = 7,
            adaptive_testing_delay=2,
            mean_adaptive_testing_delay=2,
            adaptive_testing_time_window=14,
            adaptive_testing_recall_pct=0.75,
            contact_trace_recall_pct = 0.5,
            test_FPR=0,
            test_FNR=0.1, # realized FNR will be 1 - (1 - test_FNR) * detectability
            test_schedule_proportions={(3,6): 0.3, (2,5): 0.3, (1,4): 0.3},
            non_compliance_params=(1,10), # expected non-compliance of 9%
            infectivity_alpha=1,
            use_contact_trace=True,
            use_testing=True,
            use_adaptive_testing=True):
        self.use_contact_trace = use_contact_trace
        self.use_testing = use_testing
        self.use_adaptive_testing = use_adaptive_testing
        self.agent_ids = list(range(n_agents))
        self.agents = {i:Agent(infectivity_alpha=infectivity_alpha) for i in self.agent_ids}

        for i in self.agents:
            if np.random.uniform() < init_infection_p:
                self.agents[i].start_infection(0) # maybe in future want to initialize half-way-through infections
    
        self.contact_inner_products = np.matrix(
                [[np.inner(self.agents[i].contact_vec, self.agents[j].contact_vec)
                    for j in self.agents] for i in self.agents])

        self.testing = SurveillanceTesting(self.agents, 
                                            test_FNR,
                                            mean_test_delay,
                                            mean_contact_trace_delay,
                                            test_schedule_proportions, 
                                            non_compliance_params, 
                                            contact_trace_time_window,
                                            contact_trace_recall_pct,
                                            adaptive_testing_time_window,
                                            mean_adaptive_testing_delay,
                                            adaptive_testing_recall_pct,
                                            self.contact_inner_products)
        self.curr_time_period = 0


    def step(self):
        t = self.curr_time_period
        self.step_interactions(t)
        if self.use_testing:
            new_recorded_positives = self.testing.step_test(t)
            if self.use_contact_trace:
                self.testing.step_contact_trace(t, new_recorded_positives)
            if self.use_adaptive_testing:
                self.testing.step_adaptive_testing(t, new_recorded_positives)
            self.step_isolation_removals(t)
        self.curr_time_period += 1

    def step_isolation_removals(self, t):
        for agent in self.agents.values():
            if agent.is_in_isolation and not any(agent.past_three_results):
                agent.remove_from_isolation()


    def sample_contacts(self, i):
        inner_product_vec = np.array(self.contact_inner_products[i,:].T).flatten()
        inner_product_vec[i] = 0
        normalizer = sum(inner_product_vec)
        probability_vec = inner_product_vec / normalizer

        n_contacts = self.agents[i].sample_num_contacts()

        contact_ids = np.random.choice(self.agent_ids, n_contacts, 
                                        replace=True, p=probability_vec)
        return contact_ids
    

    def get_infected_agents(self):
        return [i for i in self.agents if self.agents[i].has_infection()]


    def step_interactions(self, t):
        infected_agents = self.get_infected_agents()
        for i in infected_agents:
            contact_ids = set(self.sample_contacts(i))
            for j in contact_ids:
                if not self.agents[j].is_free_and_susceptible():
                    continue

                infectivity = self.agents[i].get_infectivity(t)

                if np.random.uniform() < infectivity:
                    self.agents[j].start_infection(t)


        


























