import numpy as np
from scipy.stats import gamma
from math import exp

GAMMA_A_PARAM = 8 # semi-arbitrarily chosen based on values in kaplan paper
MAX_INFECTION_LENGTH = 30 #also arbitrary

NB_R_DEFAULT = 2 # with beta-params (3,3) the expected value of contacts/day is 3/2 * R = 3

class Infection:
    """
    implements infection dynamics for a single agent
    """

    def __init__(self, day_infected=-1):
        self.init_day_infected = day_infected

        # use gamma distribution to track the "forward generation time density f(a)"
        # reference: https://www.medrxiv.org/content/medrxiv/early/2020/07/16/2020.07.09.20149351.full.pdf 
        # (see e.g. section 5 and Figure 6)
        self.forward_gen_time_density = np.array([gamma.pdf(x,GAMMA_A_PARAM) for x in range(MAX_INFECTION_LENGTH)])

        normalizer = sum(self.forward_gen_time_density)
        self.forward_gen_time_density = self.forward_gen_time_density / normalizer


    def get_infectivity(self, t):
        assert(self.has_infection())
        days_into_infection = t - self.init_day_infected
        if days_into_infection >= MAX_INFECTION_LENGTH:
            return 0
        else:
            return self.forward_gen_time_density[days_into_infection]


    def has_infection(self): # really this returns: is susceptible OR has been infected at some point
        if self.init_day_infected == -1:
            return False
        else:
            return True


    def start_infection(self, day_infected):
        assert(day_infected >= 0)
        self.init_day_infected = day_infected
        



class Agent:

    def __init__(self, 
            contact_vec_dim = 20,
            normalize_agent_vector=True,
            nb_r_multiplier = 1,
            nb_conjugate_beta_params=(3,3),
            contact_recall_window=4,
            use_pessimistic_detectability_curve=False):
        self.use_pessimistic_detectability_curve=use_pessimistic_detectability_curve
        # ensure first coordinate is 1 so no chance of the zero vector
        self.contact_vec = np.array([np.random.uniform() for _ in range(contact_vec_dim)])

        if normalize_agent_vector:
            norm = np.linalg.norm(self.contact_vec)
            self.contact_vec = self.contact_vec / norm

        self.contact_recall_window = contact_recall_window
        self.previous_contacts = [set([]) for _ in range(contact_recall_window)]
        
        alpha, beta = nb_conjugate_beta_params
        self.nb_p = np.random.beta(alpha, beta)

        self.is_in_isolation = False
        self.is_isolated_for_followup = False
        self.time_in_isolation = -1
        
        self.nb_r = NB_R_DEFAULT * nb_r_multiplier

        self.infection = Infection()

        self.past_three_results = [False, False, False]

        self.removed_from_isolation=False


    def sample_num_contacts(self):
        # higher contact magnitude => lower # contacts ... need to improve terminology here
        return np.random.negative_binomial(self.nb_r, self.nb_p)


    def record_contacts(self, contacts):
        self.previous_contacts.pop(0)
        self.previous_contacts.append(contacts)


    def get_avg_contacts(self):
        return self.nb_p * self.nb_r / (1 - self.nb_p)
    

    def is_free_and_susceptible(self):
        return (not self.infection.has_infection()) and (not self.is_in_isolation)


    def get_infectivity(self, t):
        return self.infection.get_infectivity(t)
    

    def get_detectability(self, t):
        if not self.has_infection():
            return 0
        else:
            viral_load = self.infection.get_infectivity(t)
            # sigmoid function with hand-tuned parameters... see design doc for now
            # for a brief justification of these parameters... future work: standardize these params
            if self.use_pessimistic_detectability_curve:
                return 1 / (1 + 1.5 * exp(-30 * (viral_load - 0.04)))
            else:
                return 1 / (1 + 1.25 * exp(-100 * (viral_load - 0.02)))


    def get_non_compliance(self):
        if self.is_in_isolation:
            return 0
        else:
            return self.non_compliance

    def configure_test_params(self, test_sched, non_compliance):
        self.test_sched = test_sched
        self.non_compliance = non_compliance
    

    def has_infection(self):
        return self.infection.has_infection()


    def start_infection(self, t):
        self.infection.start_infection(t)


    def record_test_result(self, result):
        self.past_three_results.pop(0)
        self.past_three_results.append(result)


    def isolate(self):
        self.is_in_isolation=True
        self.is_isolated_for_followup = False


    def isolate_for_followup(self):
        self.is_isolated_for_followup = True
        self.is_in_isolation = True


    def remove_from_followup_isolation(self):
        if not any(self.past_three_results) and self.is_isolated_for_followup:
            self.remove_from_isolation()


    def remove_from_isolation(self):
        assert(not any(self.past_three_results))
        self.is_in_isolation=False
        self.is_isolated_for_followup = False
        self.removed_from_isolation=True


