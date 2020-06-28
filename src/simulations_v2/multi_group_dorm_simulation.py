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
            quarantine_leakage_contacts_dorm, # dorm -> dorm contacts/day after a dorm is quarantined
            quarantine_leakage_contacts_non_dorm, # dorm -> non-dorm contacts/day after a dorm is quarantined
            quarantine_test_fraction=1.0, #testing fraction after dorm goes into quarantine
            initial_dorm_prevalence=None,
            dorm_outside_infection_p=0,
            base_config="/home/jmc678/covid_data/group-testing/src/simulations_v2/params/june8params/nominal.yaml"):
        
        self.quarantine_leakage_contacts_dorm = quarantine_leakage_contacts_dorm
        self.quarantine_leakage_contacts_non_dorm = quarantine_leakage_contacts_non_dorm
        self.quarantined_dorms = set()
        self.num_dorms = num_dorms
        self.quarantine_test_fraction = quarantine_test_fraction

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
            grp_params['daily_outside_infection_p'] = dorm_outside_infection_p
            if initial_dorm_prevalence != None:
                grp_params['initial_ID_prevalence'] = initial_dorm_prevalence
            grp_params_list.append(grp_params)

            # create next row of interactions matrix
            dorm_interactions = [inter_dorm_contacts / (num_dorms - 1)] * num_dorms + [inter_non_dorm_contacts]
            dorm_interactions[i] = intra_dorm_contacts
            interactions_mtx[i,:] = np.array(dorm_interactions)
            
            # record informative group name
            grp_names_list.append("dorm_group_{}".format(i))

        # instantiate non-dorm group
        grp_params = copy(base_params)
        grp_params['expected_contacts_per_day'] = intra_non_dorm_contacts
        grp_params['population_size'] = non_dorm_population
        grp_params['test_population_fraction'] = non_dorm_test_rate
        grp_params_list.append(grp_params)

        # solve for the non-dorm -> dorm contacts/day based on the corresponding
        # population sizes and the dorm -> non-dorm rate (assuming symmetric contacts)
        non_dorm_to_dorm_contacts = (dorm_pop_size / non_dorm_population) * inter_non_dorm_contacts
        non_dorm_interactions = [non_dorm_to_dorm_contacts] * num_dorms + [intra_non_dorm_contacts]
        interactions_mtx[num_grps - 1, :] = np.array(non_dorm_interactions)

        grp_names_list.append("non_dorm_group")

        super().__init__(grp_params_list, interactions_mtx, grp_names_list)

    
    def run_new_trajectory(self, T):
        quarantined_dorm_counts = []
        for _ in range(T):
            self.step()
            self.quarantine_infected_dorms()
            quarantined_dorm_counts.append(len(self.quarantined_dorms))
        return quarantined_dorm_counts


    def quarantine_infected_dorms(self):
        interaction_mtx = self.get_interaction_mtx()
        for i in range(self.num_dorms):
            if i in self.quarantined_dorms:
                continue
            if self.sims[i].QI > 0:
                interaction_mtx[i,:] = np.array(
                        [self.quarantine_leakage_contacts_dorm / (self.num_dorms - 1)] * self.num_dorms + \
                                               [self.quarantine_leakage_contacts_non_dorm])
                self.quarantined_dorms.add(i)
                self.sims[i].test_pop_fraction = self.quarantine_test_fraction
        self.set_interaction_mtx(interaction_mtx)













