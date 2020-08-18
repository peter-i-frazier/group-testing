from multi_group_simulation import MultiGroupSimulation
from load_params import load_params
from math import ceil
from copy import copy
import numpy as np

# Create a Multi-Group Dorm simulation where each dorm floor is its own group

class MultiGroupDormSimulation(MultiGroupSimulation):

    def __init__(self, num_dorms, # how many dorm sub-communities to include
            dorm_population, # how many students are split across the dorm sub-communities
            non_dorm_population, # how many individuals are there in the community but not in a dorm
            intra_dorm_contacts, # how many contacts/day does a student have inside their dorm
            inter_dorm_contacts, # how many contacts/day does a student have w/ other dorms
            intra_non_dorm_contacts, # how many contacts/day happens between non-dorm & dorm individuals
            inter_non_dorm_contacts, # how many contacts/day does a non-dorm person have with non-dorm people
            dorm_test_rate, # % of each dorm tested each day
            non_dorm_test_rate, # % of non-dorm tested each day
            quarantine_contacts_multiplier=0.1, # what to multiply the rate-of-contacts by after a dorm goes into quarantine
            quarantine_test_fraction=1.0, #testing fraction after dorm goes into quarantine
            initial_dorm_prevalence=None,
            dorm_outside_infection_p=0,
            safe_days_until_unquarantine=3,
            base_config="/home/jmc678/covid_data/group-testing/src/simulations_v2/params/june8params/nominal.yaml"):
        
        self.quarantined_dorms = set()
        self.num_dorms = num_dorms
        self.quarantine_test_fraction = quarantine_test_fraction

        self.quarantine_contacts_multiplier = quarantine_contacts_multiplier

        self.dorm_test_rate = dorm_test_rate

        self.safe_days_until_unquarantine = safe_days_until_unquarantine
        # this dict records the number of days since a dorm has seen a positive test
        self.days_since_last_positive = {i:0 for i in range(num_dorms)}

        _, base_params = load_params(base_config)

        grp_params_list = []
        grp_names_list = []
        num_grps = num_dorms + 1
        interactions_mtx = np.zeros((num_grps, num_grps))

        dorm_pop_size = ceil(dorm_population / num_dorms)
        

        for i in range(num_dorms):
            # create & save params for the ith dorm population
            grp_params = copy(base_params)
            grp_params['expected_contacts_per_day'] = intra_dorm_contacts
            grp_params['population_size'] = dorm_pop_size
            grp_params['test_population_fraction'] = dorm_test_rate
            if dorm_outside_infection_p != None:
                grp_params['daily_outside_infection_p'] = dorm_outside_infection_p
            if initial_dorm_prevalence != None:
                grp_params['initial_ID_prevalence'] = initial_dorm_prevalence
            grp_params_list.append(grp_params)

            # create next row of interactions matrix
            if num_dorms > 1:
                dorm_interactions = [inter_dorm_contacts / (num_dorms - 1)] * num_dorms + [inter_non_dorm_contacts]
            else:
                dorm_interactions = [0, inter_non_dorm_contacts]
            dorm_interactions[i] = intra_dorm_contacts
            interactions_mtx[i,:] = np.array(dorm_interactions)
            
            # record informative group name
            grp_names_list.append("dorm_group_{}".format(i))

        # instantiate non-dorm group
        grp_params = copy(base_params)
        if non_dorm_population > 0:
            grp_params['expected_contacts_per_day'] = intra_non_dorm_contacts
            grp_params['population_size'] = non_dorm_population
            grp_params['test_population_fraction'] = non_dorm_test_rate
            grp_params_list.append(grp_params)

        # solve for the non-dorm -> dorm contacts/day based on the corresponding
        # population sizes and the dorm -> non-dorm rate (assuming symmetric contacts)
        if non_dorm_population > 0:
            non_dorm_to_dorm_contacts = (dorm_pop_size / non_dorm_population) * inter_non_dorm_contacts
        else:
            non_dorm_to_dorm_contacts = 0
        non_dorm_interactions = [non_dorm_to_dorm_contacts] * num_dorms + [intra_non_dorm_contacts]
        interactions_mtx[num_grps - 1, :] = np.array(non_dorm_interactions)

        if non_dorm_population > 0:
            grp_names_list.append("non_dorm_group")

        super().__init__(grp_params_list, interactions_mtx, grp_names_list)

    
    def run_new_trajectory(self, T):
        quarantined_dorm_counts = []
        for _ in range(T):
            self.step()
            self.update_quarantine_status()
            quarantined_dorm_counts.append(len(self.quarantined_dorms))
        return quarantined_dorm_counts

    
    def update_quarantine_status(self):
        self.update_days_since_last_positive()
        for i in range(self.num_dorms):
            if self.days_since_last_positive[i] == 0 and i not in self.quarantined_dorms:
                self.quarantine_dorm(i)
            elif i in self.quarantined_dorms and self.days_since_last_positive[i] >= self.safe_days_until_unquarantine:
                self.unquarantine_dorm(i)
    

    def update_days_since_last_positive(self):
        for i in range(self.num_dorms):
            if self.sims[i].new_QI_from_self_reports > 0 or self.sims[i].new_QI_from_last_test > 0:
                self.days_since_last_positive[i] = 0
            else:
                self.days_since_last_positive[i] += 1


    def quarantine_dorm(self, i):
        assert(i not in self.quarantined_dorms)
        interaction_mtx = self.get_interaction_mtx()
        for j in range(self.num_dorms+1):
            interaction_mtx[i,j] = self.quarantine_contacts_multiplier * interaction_mtx[i,j]
            if i != j:
                interaction_mtx[j,i] = self.quarantine_contacts_multiplier * interaction_mtx[j,i]

        self.sims[i].test_pop_fraction = self.quarantine_test_fraction

        self.quarantined_dorms.add(i) 
        self.set_interaction_mtx(interaction_mtx)


    def unquarantine_dorm(self, i):
        assert(i in self.quarantined_dorms)
        interaction_mtx = self.get_interaction_mtx()
        for j in range(self.num_dorms+1):
            interaction_mtx[i,j] = interaction_mtx[i,j] / self.quarantine_contacts_multiplier 
            if i != j:
                interaction_mtx[j,i] = interaction_mtx[j,i] / self.quarantine_contacts_multiplier 

        self.sims[i].test_pop_fraction = self.dorm_test_rate

        self.quarantined_dorms.remove(i)
        self.set_interaction_mtx(interaction_mtx)


