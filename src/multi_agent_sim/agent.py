class Agent:

    def __init__(self, agent_id):
        self.id = agent_id
        self.params = {}

    def add_param(self, key, val):
        self.params[key] = val

    def get_param(self, key, default=None):
        return self.params.get(key, default)
