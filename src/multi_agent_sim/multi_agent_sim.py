from infection_dynamics import InfectionDynamics
from social_network import SocialNetwork
from agent import Agent
from quarantine_manager import QuarantineManager
from contact_tracing import ContactTracing
from adaptive_testing import AdaptiveTesting
from surveillance_testing import SurveillanceTesting

class MultiAgentSim:

    def __init__(self, 
            main_params, 
            infection_dynamics_params,
            social_network_params,
            ct_params,
            at_params,
            st_params):
        self.n_agents = main_params['n_agents']
        self.agents = {idx:Agent(idx) for idx in range(self.n_agents)}
        self.infection = InfectionDynamics(self.agents, infection_dynamics_params)
        self.network = SocialNetwork(self.agents, social_network_params)
        self.quarantine = QuarantineManager(self.agents)
        self.use_ct = main_params['use_contact_tracing']
        self.ct = ContactTracing(self.agents, self.quarantine, ct_params)
        self.curr_day = 0

        self.use_at = main_params['use_adaptive_testing']
        self.at = AdaptiveTesting(self.agents, self.quarantine, self.network, self.infection, at_params)

        self.use_st = main_params['use_surveillance_testing']
        self.st = SurveillanceTesting(self.agents, self.infection, self.quarantine, st_params)


    def step(self):
        self.curr_day += 1

        # step 1: simulate contacts & transmission dynamics for the day

        # iterate over all free & infected agents
        inf_agent_ids = self.infection.get_infected_agent_ids(self.curr_day)
        for agent_id in inf_agent_ids:
            if self.quarantine.is_agent_isolated(agent_id):
                continue

            # sample other agents w/ whom agent_id comes into contact
            contact_ids = self.network.sample_contacts(agent_id)

            # sample subset of contact agents who become infected based on their interaction w/ agent_id
            infectee_ids = self.infection.sample_transmissions(agent_id, contact_ids, self.curr_day)

            # record contacts for contact tracing
            if self.use_ct:
                self.ct.record_contacts(agent_id, contact_ids, self.curr_day)

            # start infection dynamics for newly infected agents
            for infectee_id in infectee_ids:
                if not self.quarantine.is_agent_isolated(infectee_id):
                    self.infection.infect_agent(infectee_id, self.curr_day)

        # next, sample self-reports from each agent
        self_report_ids = self.infection.sample_self_report_isolations(inf_agent_ids, self.curr_day)
        for self_report_id in self_report_ids:

            # when an agent self-reports they immediately go into quarantine
            if not self.quarantine.is_agent_isolated(self_report_id):
                self.quarantine.isolate_agent(self_report_id, self.curr_day)

            # we also begin the contact trace process for them
            if self.use_ct:
                self.ct.add_agent_to_trace_queue(self_report_id, self.curr_day)

            if self.use_at:
                self.at.trigger_adaptive_test(self_report_id, self.curr_day)

        # run surveillance tests
        st_positive_ids = self.st.step_surveillance_tests(self.curr_day)

        if self.use_ct:
            for agent_id in st_positive_ids:
                self.ct.add_agent_to_trace_queue(agent_id, self.curr_day)

        if self.use_at:
            for agent_id in st_positive_ids:
                self.at.trigger_adaptive_test(agent_id, self.curr_day)

        # process quarantine removals
        self.quarantine.step_isolation_removals(self.curr_day)

        # process contact-traces for the day
        if self.use_ct:
            self.ct.step_trace_queue(self.curr_day)

        if self.use_at:
            self.at.step_adaptive_tests(self.curr_day)


