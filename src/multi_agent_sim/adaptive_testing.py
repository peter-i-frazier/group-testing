import numpy as np

def sample_delay(distn):
    values = list(range(len(distn)))
    return np.random.choice(values, 1, p=distn)[0]

class AdaptiveTesting:

    def __init__(self, agents, quarantine_manager, social_network, infection_dynamics, params):
        self.agents = agents
        self.social_network = social_network
        self.infection_dynamics = infection_dynamics
        self.quarantine_manager = quarantine_manager
        self.delay_dist = params['at_delay_distribution']
        # when an agent triggers an adaptive test, we look at their average contacts per day 
        # and multiply that avg # contacts by this multiplier, and then cast a net whose size is equal to the product
        # meaning that we take the (product) closest individuals to the triggering agent in the social network
        self.net_size_contact_multiplier = params['at_net_size_contact_multiplier']
        self.recall_rate = params['at_recall_rate']
        self.num_adaptive_tests = 0

        self.tests_to_run_by_day = {}
        self.test_days_stepped = []


    # this needs to be called before step_adaptive_tests() for the same day
    # is called
    def trigger_adaptive_test(self, agent_id, day):
        assert(day not in self.test_days_stepped)
        mean_contacts = self.social_network.get_mean_contacts(agent_id)

        net_size = int(mean_contacts * self.net_size_contact_multiplier)

        closest_agents = self.social_network.get_n_closest_agents(agent_id, net_size)

        for agent in closest_agents:
            delay = sample_delay(self.delay_dist)
            test_day = day + delay
            if test_day not in self.tests_to_run_by_day:
                self.tests_to_run_by_day[test_day] = []
            self.tests_to_run_by_day[test_day].append(agent)


    def step_adaptive_tests(self, day):
        self.test_days_stepped.append(day)
        if day not in self.tests_to_run_by_day:
            return []

        agents_to_test = self.tests_to_run_by_day[day]
        self.num_adaptive_tests += len(agents_to_test)
        positive_agents = []

        for agent in agents_to_test:
            if self.infection_dynamics.sample_test_result(agent, day):
                positive_agents.append(agent)
                self.quarantine_manager.isolate_agent(agent, day)

        return positive_agents



