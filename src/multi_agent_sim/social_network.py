import numpy as np

NB_R_DEFAULT=2
NB_ALPHA_DEFAULT=3
NB_BETA_DEFAULT=3


class SocialNetwork:

    def __init__(self, agents, params):
        self.agents = agents
        self.agent_ids = list(range(len(agents)))

        # read social network parameters
        location_vec_dim = params.get('location_vec_dim', 2)
        network_gamma = params.get('network_gamma', 5)

        self.daily_contacts_distn_type = params.get('daily_contacts_distn_type', 'negative_binomial')
        if self.daily_contacts_distn_type != 'negative_binomial':
            raise(Exception("Only Negative Binomial currently supported for daily contact distribution"))
        nb_r = params.get('neg_bin_r', NB_R_DEFAULT)
        nb_alpha, nb_beta = params.get('neg_bin_p_hyperparams', (NB_ALPHA_DEFAULT, NB_BETA_DEFAULT))
        
        # generate agent parameters
        for agent in self.agents.values():
            location_vec = np.array([np.random.uniform() for _ in range(location_vec_dim)])
            agent.add_param('location_vec', location_vec)
            nb_p = np.random.beta(nb_alpha, nb_beta)
            agent.add_param('nb_r', nb_r)
            agent.add_param('nb_p', nb_p)


        # precompute agent similarities
        self.agent_similarities = np.matrix(
            [[1 / (np.linalg.norm(self.agents[i].get_param('location_vec') - 
                    self.agents[j].get_param('location_vec'))**network_gamma) if i != j else 0
                        for j in self.agents] 
                        for i in self.agents])


    def get_mean_contacts(self, agent_id):
        agent = self.agents[agent_id]
        r = agent.get_param('nb_r')
        p = agent.get_param('nb_p')
        return p * r / (1 - p)


    def get_n_closest_agents(self, agent_id, n):
        if n == 0:
            return []
        similarities = np.array(self.agent_similarities[agent_id,:].T).flatten()
        nth_largest = sorted(similarities, reverse=True)[n-1]
        closest_agents = []
        for idx, similarity in enumerate(similarities):
            if similarity >= nth_largest:
                closest_agents.append(idx)
        return closest_agents


    """
    sample a subset of the population with whom a specific agent comes into contact in any time period
    pays no attention to agent isolation status or infection status
    """
    def sample_contacts(self, agent_idx):
        agent = self.agents[agent_idx]
        nb_r = agent.get_param('nb_r')
        nb_p = agent.get_param('nb_p')

        num_contacts = min(np.random.negative_binomial(nb_r, nb_p), len(self.agent_ids)-1)

        similarity_vec = np.array(self.agent_similarities[agent_idx,:].T).flatten()
        similarity_vec[agent_idx] = 0

        normalizer = sum(similarity_vec)
        probability_vec = similarity_vec / normalizer

        contact_ids = np.random.choice(self.agent_ids, num_contacts, 
                                        replace=False, p=probability_vec)
        return contact_ids


