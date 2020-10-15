class QuarantineManager:
    def __init__(self, agents, max_iso_length=21):
        self.agents = agents
        self.max_iso_length=max_iso_length


    def is_agent_isolated(self, agent_id):
        return self.agents[agent_id].get_param('is_in_isolation', False)


    def isolate_agent(self, agent_id, curr_day):
        self.agents[agent_id].add_param('is_in_isolation', True)
        self.agents[agent_id].add_param('day_isolated', curr_day)
        self.agents[agent_id].add_param('day_freed', None)


    def free_agent(self, agent_id, curr_day):
        self.agents[agent_id].add_param('is_in_isolation', False)
        self.agents[agent_id].add_param('day_free', curr_day)


    def step_isolation_removals(self, curr_day):
        for agent_id in self.agents:
            if self.agents[agent_id].get_param('is_in_isolation', False):
                day_isolated = self.agents[agent_id].get_param('day_isolated')
                if curr_day - day_isolated >= self.max_iso_length:
                    self.free_agent(agent_id, curr_day)

