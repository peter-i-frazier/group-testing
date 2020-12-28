"""
This class is responsible for tracking the infection status over all agents
in the population across time.  
"""

from math import exp
import numpy as np
from scipy.stats import gamma as gamma_dist

MAX_INF_LENGTH = 30


OPT_DETECT_PARAMS = (1.25,100,0.02)
PESS_DETECT_PARAMS =  (1.5,30,0.04)

# helper function to discretize gamma distribution, used for forward generation density
def gen_disc_gamma(max_x, alpha_shape, beta_rate):
    pdf_vals = []
    for x in range(1, max_x+1):
        pdf_vals.append(gamma_dist.pdf(x * beta_rate, alpha_shape))
    normalizer = sum(pdf_vals)
    return [pdf_val / normalizer for pdf_val in pdf_vals]


class InfectionDynamics:
    def __init__(self, agents, params):
        self.agents = agents

        # read infection dynamics parameters
        
        self.forward_gen_alpha = params.get('forward_gen_alpha', 8)
        self.forward_gen_beta_hyperparams = params.get('forward_gen_beta_hyperparams', (5,5))

        self.detectability_curve_type = params.get('detectability_curve_type', 'optimistic')
        assert(self.detectability_curve_type in ['optimistic', 'pessimistic'])
        if self.detectability_curve_type == 'optimistic':
            self.detectability_params = OPT_DETECT_PARAMS
        else:
            self.detectability_params = PESS_DETECT_PARAMS

        # on infection day t, an agent self-reports with probability
        #       multiplier * f(t - delay)
        # where f is the forward-gen density for that agent.
        # the self-reporting agent immediatley has a test performed, and they isolate iff test comes back positive
        # (the word 'delay' is a bit of a misnomer in this context)
        self.self_reporting_multiplier = params.get('self_reporting_multiplier', 0.8)
        self.self_reporting_delay = params.get('self_reporting_delay', 3)

        self.outside_infection_rate = params.get('outside_infection_rate', 0)
        
        init_infection_p = params.get('init_infection_rate', 0.001)
        deterministic = params.get('use_deterministic_infection_counts', False)

        # initialize agent infections
        n_agents = len(self.agents)
        for i in self.agents:
            if deterministic and i / n_agents <= init_infection_p:
                self.infect_agent(i, 0)
            elif (not deterministic) and np.random.uniform() <= init_infection_p:
                self.infect_agent(i, 0)


    def step_outside_infections(self, day):
        for agent_id in self.agents:
            if self.agents[agent_id].get_param('day_infection_started') != None:
                continue
            if np.random.uniform() <= self.outside_infection_rate:
                self.infect_agent(agent_id, day)


    def sample_transmissions(self, infector_id, contact_ids, day):
        infector = self.agents[infector_id]
        day_inf_started = infector.get_param('day_infection_started')
        assert(day_inf_started != None)

        if day - day_inf_started > MAX_INF_LENGTH:
            return []
        else:
            day_since_started_idx = day - day_inf_started - 1

        disc_gamma = infector.get_param('forward_gen_disc_gamma')
        transmission_p = disc_gamma[day_since_started_idx]
        return [contact_id for contact_id in contact_ids if np.random.uniform() <= transmission_p]


    def sample_self_reports(self, agent_ids, curr_day):
        self_reports = []
        for agent_id in agent_ids:
            agent = self.agents[agent_id]
            if agent.get_param('day_infection_started') == None:
                continue
            if curr_day - agent.get_param('day_infection_started') > MAX_INF_LENGTH:
                continue

            self_report_day_idx = curr_day - agent.get_param('day_infection_started') - self.self_reporting_delay
            if self_report_day_idx < 0:
                continue
            disc_gamma = agent.get_param('forward_gen_disc_gamma')
            self_report_p = self.self_reporting_multiplier * disc_gamma[self_report_day_idx]
            if np.random.uniform() < self_report_p:
                self_reports.append(agent_id)

        return self_reports
    

    def sample_self_report_isolations(self, agent_ids, curr_day):
        self_reports = self.sample_self_reports(agent_ids, curr_day)
        return [agent_id for agent_id in self_reports if self.sample_test_result(agent_id, curr_day)]


    def sample_test_result(self, agent_id, curr_day):
        if self.agents[agent_id].get_param('day_infection_started') == None:
            return False
        A, B, C = self.detectability_params
        day_idx = curr_day - self.agents[agent_id].get_param('day_infection_started') - 1
        if day_idx >= MAX_INF_LENGTH:
            return False
        disc_gamma = self.agents[agent_id].get_param('forward_gen_disc_gamma')
        transmission_p = disc_gamma[day_idx]

        detectability_p = 1 / (1 + A * exp(-B * (transmission_p - C)))
        return np.random.uniform() <= detectability_p


    def infect_agent(self, agent_id, day):
        if self.agents[agent_id].get_param('day_infection_started') != None:
            return

        self.agents[agent_id].add_param('day_infection_started', day)
        self.agents[agent_id].add_param('forward_gen_alpha', self.forward_gen_alpha)

        # sample the forward-gen beta param from a gamma prior; used saved hyperparams
        hyper_alpha, hyper_beta = self.forward_gen_beta_hyperparams
        beta = np.random.gamma(hyper_alpha, 1/hyper_beta)
        self.agents[agent_id].add_param('forward_gen_beta', beta)

        # compute the (discretized) forward-generation gamma density using appropriate alpha & beta params
        disc_gamma = gen_disc_gamma(MAX_INF_LENGTH, self.forward_gen_alpha, beta)
        self.agents[agent_id].add_param('forward_gen_disc_gamma', disc_gamma)


    def get_infected_agent_ids(self, day):
        return [agent_id for agent_id, agent in self.agents.items() 
                if agent.get_param('day_infection_started') != None and 
                    day - agent.get_param('day_infection_started') <= MAX_INF_LENGTH]


    def get_cum_infected_agent_ids(self):
        return [agent_id for agent_id, agent in self.agents.items()
                if agent.get_param('day_infection_started') != None]
