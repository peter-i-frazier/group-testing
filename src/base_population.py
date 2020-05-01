import random
"""
Implements a generic parent class for a Population object for use in our COVID-testing simulations.

A Population object is responsible for the following:
    * Keeping track of infection-status and quarantine-status across the entire
        population
    * Simulating the day-by-day progression of the disease
    * Making quarantine/social distancing restrictions in response to
        test outcomes

I am currently envisioning the following subclasses:
    * IndividualInteractionPopulation(BasePopulation)
        - implements the population dynamics with pairwise interaction frequency
          between individuals sampled each day from a Poisson distribution
    * HouseholdPopulation(BasePopulation)
        - implements the population dynamics where probability of a new infection
          is a function of current prevalence, with inter-household infection dynamics
    * FluidHouseholdPopulation(HouseholdPopulation)
        - implements the fluid population dynamics, with inter-household infection dynamics
          used to inform the group testing protocol.  i.e. implements the dynamics that
          are currently contained in dynamic_protocol_design.py

BasePopulation has the following attributes:
    * n_agents: how many agents are present in the population. an agent
                may correspond to an individual or a household, depending on
                how the subclass is implemented
    * agents: a set of agent-identifiers
    * quarantine_status: a dictionary mapping agent IDs to a boolean; True indicates
                         the agent is quarantined, False indicates otherwise
    * infection_status: a dict mapping agent IDs to booleans indicating whether or not
                        the agent is infected
    * disease_length: number of days the disease lasts in an infected individual
    * quarantine_length: number of days quarantine lasts
    * days_infected_so_far: dictionary mapping agent IDs to integers indicating number
                            of days the agent has been infected so far
    * days_quarantined_so_far: dict mapping agent IDs to integer, indicating number of
                            days agent has been in quarantine so far
    * days_until_symptomatic: number of days until an individual with the disease presents
                                symptoms

BasePopulation implements the following methods:
    * step(): run one day of disease progression
    * track_infection_status(): increase days_infected_so_far count for each
                                    infected individual, resolve disease for agents who have
                                    reached disease_length
    * track_quarantine_status(): increase days_quarantined_so_far for each infected individual,
                                    resolve quarantine for agents who have reached
                                    quarantine_length
    * resolve_infection(agent_id): implement logic about what happens when disease has 
                                    run its course
    * resolve_quarantine(agent_id): implement logic about what happens when quarantine ends

Probably should add the following methods in the future:
    * quarantine_agent(agent_id)
    * is_agent_symptomatic(agent_id)
    * is_agent_infected(agent_id)
    * is_agent_quarantined(agent_id)
    * get_num_infected()
    * get_num_quarantined()
"""

class BasePopulation:


    def __init__(self, n_agents, 
                        disease_length, 
                        quarantine_length, 
                        days_until_symptomatic, 
                        initial_prevalence = 0):
        self.n_agents = n_agents
        self.agents = set(range(n_agents))
        self.quarantine_status = {agent_id:False for agent_id in self.agents}
        self.infection_status = {agent_id: random.random() < initial_prevalence 
                                    for agent_id in self.agents}
        self.days_infected_so_far = {agent_id: 0 for agent_id in self.agents}
        self.days_quarantined_so_far = {agent_id: 0 for agent_id in self.agents}
        self.days_until_symptomatic = days_until_symptomatic
        self.disease_length = disease_length
        self.quarantine_length = quarantine_length
        self.cumulative_infected_agents = set()


    def get_summary_data(self):
        raise(Exception("todo"))


    def get_summary_data(self):
        return {'num_quarantined':self.get_num_quarantined(),
                'num_infected': self.get_num_infected(),
                'cumulative_num_infected': self.get_cumulative_num_infected()}


    def track_infection_status(self):
        for agent_id in self.agents:
            if self.infection_status[agent_id]:
                self.cumulative_infected_agents.add(agent_id)
                self.days_infected_so_far[agent_id] += 1
                if self.days_infected_so_far[agent_id] >= self.disease_length:
                    self.resolve_infection(agent_id)


    def track_quarantine_status(self):
        for agent_id in self.agents:
            if self.quarantine_status[agent_id]:
                self.days_quarantined_so_far[agent_id] += 1
                if self.days_quarantined_so_far[agent_id] >= self.quarantine_length:
                    self.resolve_quarantine(agent_id)


    def resolve_infection(self, agent_id):
        self.infection_status[agent_id] = False
        self.days_infected_so_far[agent_id] = 0
    

    def resolve_quarantine(self, agent_id):
        self.quarantine_status[agent_id] = False
        self.days_quarantined_so_far[agent_id] = 0
    

    def step(self):
        raise(Exception("step() logic must be implemented by child class"))

    
    def infect_agent(self, agent_id):
        if agent_id not in self.cumulative_infected_agents:
            self.infection_status[agent_id] = True
            self.cumulative_infected_agents.add(agent_id)

   
    def quarantine_agent(self, agent_id):
        self.quarantine_status[agent_id] = True


    def is_agent_infected(self, agent_id):
        # might need to optimize this for speed in the future
        return self.infection_status[agent_id]


    def is_agent_quarantined(self, agent_id):
        return self.quarantine_status[agent_id]


    def is_agent_recovered(self, agent_id):
        return (not self.infection_status[agent_id]) and \
                (agent_id in self.cumulative_infected_agents)


    def get_num_quarantined(self):
        return sum(self.quarantine_status.values())
    

    def get_num_infected(self): 
        return sum(self.infection_status.values())


    def is_agent_symptomatic(self, agent_id):
        return self.days_infected_so_far[agent_id] >= self.days_until_symptomatic


    def get_cumulative_num_infected(self):
        return len(self.cumulative_infected_agents)
