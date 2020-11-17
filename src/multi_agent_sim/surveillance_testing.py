import numpy as np

def sample_delay(distn):
    values = list(range(len(distn)))
    return np.random.choice(values, 1, p=distn)[0]

class SurveillanceTesting:

    def __init__(self, agents, infection_dynamics, quarantine_manager, params):
        self.agents = agents
        self.infection_dynamics = infection_dynamics
        self.quarantine_manager = quarantine_manager
        self.testing_window = params['st_testing_window']
        self.missed_test_rate = params['st_missed_test_rate']

        # initialize agent test schedules
        uniform_dist = [1/self.testing_window] * self.testing_window
        for agent in self.agents.values():
            test_day = sample_delay(uniform_dist)
            agent.add_param('test_day', test_day)
        
        self.num_surveillance_tests = 0

    def step_surveillance_tests(self, day):
        curr_day = day % self.testing_window
        positive_agent_ids = []
        for agent_id, agent in self.agents.items():
            if agent.get_param('test_day') != curr_day:
                continue
            if np.random.uniform() <= self.missed_test_rate:
                continue
            if self.infection_dynamics.sample_test_result(agent_id, day):
                positive_agent_ids.append(agent_id)
                self.quarantine_manager.isolate_agent(agent_id, day)
            self.num_surveillance_tests += 1

        return positive_agent_ids

