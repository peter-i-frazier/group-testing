from infection_dynamics import InfectionDynamics
from social_network import SocialNetwork
from agent import Agent
from quarantine_manager import QuarantineManager

class MultiAgentSim:

    def __init__(self, n_agents, infection_dynamics_params, social_network_params):
        self.n_agents = n_agents
        self.agents = {idx:Agent(idx) for idx in range(n_agents)}
        self.infection = InfectionDynamics(self.agents, infection_dynamics_params)
        self.network = SocialNetwork(self.agents, social_network_params)
        self.quarantine = QuarantineManager(self.agents)
        self.curr_day = 0


    def step(self):
        self.curr_day += 1

        # start by simulating contacts & transmissions for infected agents

        inf_agent_ids = self.infection.get_infected_agent_ids(self.curr_day)

        for agent_id in inf_agent_ids:
            if self.quarantine.is_agent_isolated(agent_id):
                continue

            contact_ids = self.network.sample_contacts(agent_id)
            infectee_ids = self.infection.sample_transmissions(agent_id, contact_ids, self.curr_day)

            for infectee_id in infectee_ids:
                if not self.quarantine.is_agent_isolated(infectee_id):
                    self.infection.infect_agent(infectee_id, self.curr_day)

        # next, sample self-reports from each agent
        self_report_ids = self.infection.sample_self_report_isolations(inf_agent_ids, self.curr_day)
        for self_report_id in self_report_ids:
            self.quarantine.isolate_agent(self_report_id, self.curr_day)

        self.quarantine.step_isolation_removals(self.curr_day)


