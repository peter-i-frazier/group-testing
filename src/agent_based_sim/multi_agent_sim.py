"""
Implements the agent-based version of the simulation described
here: 
https://docs.google.com/document/d/17tm2uPOUxrcfXG-W2paw02K6pSx_iDMXY51n5RSjTMI/edit#
"""
"""
TODO:
    * contact_trace_recall_rate
"""

import numpy as np

from agent import Agent
from surveillance_testing import SurveillanceTesting, debug
 
class MultiAgentSim:

    def __init__(self, 
            n_agents, 
            init_infection_p,
            deterministic_agent_count=False,
            surveillance_test_delay_distn=[0.2,0.7,0.1],
            contact_tracing_delay_distn=[0.2,0.5,0.3],
            contact_tracing_recall_window = 4,
            contact_tracing_recall_rate=1,
            adaptive_testing_delay_distn=[0.2,0.3,0.3,0.2],
            adaptive_testing_time_window=10,
            adaptive_testing_recall_rate=0.75,
            test_FPR=0,
            test_FNR=0.1, # realized FNR will be 1 - (1 - test_FNR) * detectability
            test_schedule_proportions={(3,6): 0.3, (2,5): 0.3, (1,4): 0.3},
            non_compliance_params=(1,10), # expected non-compliance of 9%
            nb_r_multiplier=1,
            use_contact_trace=True,
            use_testing=True,
            use_adaptive_testing=True,
            use_pessimistic_detectability_curve=False,
            record_dataset=False,
            use_norm_over_innerprod=False,
            agent_dimensionality=20,
            normalize_agent_vector=True):
        self.record_dataset = record_dataset
        if self.record_dataset:
            self.recorded_contacts = {}
            self.recorded_isolation_statuses = {}
        self.use_contact_trace = use_contact_trace
        self.use_testing = use_testing
        self.use_adaptive_testing = use_adaptive_testing
        self.agent_ids = list(range(n_agents))
        self.agents = {i:Agent(nb_r_multiplier=nb_r_multiplier, 
                                use_pessimistic_detectability_curve=use_pessimistic_detectability_curve,
                                contact_vec_dim=agent_dimensionality,
                                normalize_agent_vector = normalize_agent_vector) 
                                for i in self.agent_ids}

        for i in self.agents:
            if deterministic_agent_count and (i+1) / n_agents <= init_infection_p:
                self.agents[i].start_infection(0)
            elif not deterministic_agent_count and  np.random.uniform() < init_infection_p:
                self.agents[i].start_infection(0) # maybe in future want to initialize half-way-through infections
    
        if use_norm_over_innerprod:
            self.contact_inner_products = np.matrix(
                    [[1 / (1e-16 + np.linalg.norm(self.agents[i].contact_vec - self.agents[j].contact_vec)**16) 
                        for j in self.agents] for i in self.agents])
        else:
            self.contact_inner_products = np.matrix(
                    [[np.inner(self.agents[i].contact_vec, self.agents[j].contact_vec)
                        for j in self.agents] for i in self.agents])

        self.testing = SurveillanceTesting(self.agents, 
                                            test_FNR,
                                            surveillance_test_delay_distn,
                                            contact_tracing_delay_distn,
                                            contact_tracing_recall_rate,
                                            test_schedule_proportions, 
                                            non_compliance_params, 
                                            adaptive_testing_time_window,
                                            adaptive_testing_delay_distn,
                                            adaptive_testing_recall_rate,
                                            self.contact_inner_products)
        self.curr_time_period = 0


    def step(self):
        t = self.curr_time_period
        
        if self.record_dataset:
            self.recorded_isolation_statuses[t] = {}
            for agent_id in self.agents:
                if self.agents[agent_id].is_in_isolation:
                    self.recorded_isolation_statuses[t][agent_id] = self.agents[agent_id].is_in_isolation

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
            if agent.is_in_isolation and not any(agent.past_three_results) and not agent.is_isolated_for_followup:
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
    

    def get_free_infected_agents(self):
        return [i for i in self.agents if self.agents[i].has_infection() and not self.agents[i].is_in_isolation]


    def get_infected_agents(self):
        return [i for i in self.agents if self.agents[i].has_infection()]


    def step_interactions(self, t):
        infected_agents = self.get_free_infected_agents()
        total_contacts = 0

        if self.record_dataset:
            self.recorded_contacts[t] = {}

        for i in infected_agents:
            contact_ids = set([j for j in self.sample_contacts(i) if not self.agents[j].is_in_isolation])
            self.agents[i].record_contacts(contact_ids)
            total_contacts += len(contact_ids)
            
            if self.record_dataset:
                self.recorded_contacts[t][i] = contact_ids

            for j in contact_ids:
                if self.agents[j].is_in_isolation:
                    continue

                infectivity = self.agents[i].get_infectivity(t)

                if np.random.uniform() < infectivity:
                    self.agents[j].start_infection(t)
        debug("interactions at time {}: there were {} free & infected agents, and they interacted with {} free individuals".format(t,
            len(infected_agents), total_contacts))
            

        


























