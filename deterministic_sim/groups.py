from typing import List
import numpy as np
'''
This module contains code for working with groups and meta-groups
'''

class meta_group:
    '''
    A "meta-group" is a collection of groups.
    Each group within the meta-group is associated with a level of contact.

    Given that a contact exposes a member of a meta-group, the probability that the person exposed is within group i
    is proportional to the population of group i times group i's level of contact.

    Examples of meta-groups: UG, Grad-Professional, Staff/Faculty, Grad-Research
    '''
    def __init__(self, name, pop, contact_units):
        assert (len(pop) == len(contact_units))
        self.name = name # Name of the meta-group, e.g., "UG"
        self.contact_units = contact_units # Amount of contact that each group has
        self.K = len(contact_units) # Number of groups
        self.pop = pop # Number of people in each group

    # used to be called well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact)
    def infection_matrix(self, infections_per_contact_unit):
        '''
        Returns an infection rate matrix that corresponds to well-mixed interactions within the meta-group.

        The matrix returned has entry [i,j] equal to the total number of secondary infections expected in
        group j that result from a single infection in group i, assuming that each contact_unit results in
        infections_per_contact_unit total secondary infections.

        A person in group i has an amount of contact during their infectious period (summed over all exposed groups)
        equal to self.contact_units[i].

        The units in "contact units" can be total contacts over the course of someone's infectious period,
        in which case infections_per_contact is the same as the probability of transmission given contact.
        It could also be the *rate* at which people have contact per day, in which case infections_per_contact
        should be the probability of transmission given contact times the expected length of a person's
        infectious period.
        '''

        # An incoming contact lands in a population with a probability proportional to its number of outgoing contacts,
        # which is pop[i]*self.contact_units[i].
        q = self.pop * self.contact_units / np.sum(self.pop * self.contact_units)

        # The total number of people infected by a positive in group i is infections_per_contact * marginal_contacts[i]
        r = infections_per_contact_unit * self.contact_units

        # When a person in group i is infected, the number of people infected in group j is then r[i]*q[j]
        infection_rates = np.outer(r, q)
        return infection_rates


    def get_init_SIR(self, initial_infectious: float, initial_recovered: float,
                     weight: str = "population"):
        """Return initial SIR vectors.

        Specify a weight used to distribute the initial infections and initial
        recovered across the groups within this metagroup.

        Use "population" if each person is equally likely to be infected.
        Use "population x contacts" if a person's probability of being infected
        is proportional to their amount of contact.
        Use "most social" if we know that the initial infections were among the most social.

        Args:
            initial_infectious (float): Initial infectious count.
            initial_recovered (float): Initial recovered count.
            weight (str): {population, population x contacts, most_social}
        """
        if weight == "population":
            w = self.pop
        elif weight == "population x contacts":
            w = self.contact_units * self.pop
        elif weight == "most_social":
            w = np.zeros(self.K)
            w[-1] = 1
        else:
            raise ValueError("The provided weight is not supported.")
        w = w / np.sum(w)  # normalize
        R0 = initial_recovered * w
        I0 = initial_infectious * w
        S0 = np.maximum(self.pop - R0 - I0, 0)
        return S0, I0, R0


class population:
    '''
    A population is a collection of meta-groups
    '''

    def __init__(self, meta_group_list: List[meta_group], meta_group_contact_matrix: np.ndarray):
        '''
        meta_group_matrix[i,j] is a matrix that gives the number of contacts that occurred where the source case was
        in meta-group i and the exposed group was in meta-group j.

        We normalize each row so that it sums to 1, so that meta_group[i,j] is the conditional probability that
        the exposed is in meta-group j, given that the source is in meta-group i.  Note that these conditional
        probabilities may be influenced by the population sizes -- if meta-group j is bigger, then it may be more
        likely to be the exposed meta-group.
        '''
        self.meta_group_contact_matrix = meta_group_contact_matrix
        self.meta_group_list = meta_group_list

        n, m = len(meta_group_contact_matrix), len(meta_group_contact_matrix[0])
        assert (n == m)

        k = len(meta_group_list)
        assert (k == n)

        self.idx2groupname = []
        for meta_group in meta_group_list:
            for i in range(meta_group.K):
                self.idx2groupname.append(meta_group.name + " " + str(i))

        self.groupname2idx = {name: i for i, name in enumerate(self.idx2groupname)}


    def infection_matrix(self, infections_per_contact_unit):
        '''
        Returns an infection rate matrix that corresponds to well-mixed interactions within each meta-group,
        and interactions across meta-groups given by self.meta_group_contact_matrix.

        The matrix returned has entry [i,j] equal to the total number of secondary infections expected in
        group j that result from a single infection in group i, assuming that each contact_unit results in
        infections_per_contact_unit total secondary infections.


        Here, a "group" is an integer that corresponds to a meta-group and a group within that meta-group.
        '''
        dim_tot = 0 #total number of meta-group-groups
        cum_tot = [] #help keep track of location in res matrix
        for i in self.meta_group_list:
            cum_tot.append(dim_tot)
            dim_tot += i.K

        res = np.zeros((dim_tot, dim_tot))
        for i in range(len(self.meta_group_list)): #source meta group
            for j in range(self.meta_group_list[i].K): #source meta-group-group
                for k in range(len(self.meta_group_list)): #exposed meta group
                    for l in range(self.meta_group_list[k].K): #expoed meta-group-group
                        q = self.meta_group_list[k].pop[l] * self.meta_group_list[k].contact_units[l] / np.sum(self.meta_group_list[k].pop * self.meta_group_list[k].contact_units)
                        res[cum_tot[i]+j, cum_tot[k]+l] = \
                        infections_per_contact_unit[i] * self.meta_group_list[i].contact_units[j] * self.meta_group_contact_matrix[i,k] * q

        return res

    def metagroup_names(self):
        '''
        Returns a list of the names of the metagroups in this population
        '''
        return [mg.name for mg in self.meta_group_list]


    def idx_to_groupname(self, i):
        '''
        Returns a string naming the group indexed by i, e.g., "UG 6".

        The name is a concatenation of the group's meta-group name and the group's number of contact units
        '''
        return self.idx2groupname[i]

    def metagroup_indices(self, metagroup_names):
        '''
        Returns a list of group indices corresponding to the passed metagroup_names
        '''
        res = []
        for metagroup_name in metagroup_names:
            tmp = []
            for i in range(len(self.idx2groupname)):
                if self.idx2groupname[i][0: len(metagroup_name)] == metagroup_name:
                    tmp.append(i)
            res.append(tmp)
        return res


    def groupname_to_idx(self, groupname):
        '''
        Returns the group index associated with a groupname
        '''
        return self.groupname2idx[groupname]

    def flatten(self, input):
        '''
        Returns a flattened version of inputted array
        amenable to usage in sim(), which only takes 1D SIR etc.
        e.g. marginal contacts, or population counts per meta-group-group
        '''
        tmp = []
        for i in range(len(input)):
            tmp += list(input[i])
        return np.array(tmp)


    def get_init_SIR(self, initial_infectious, initial_recovered):
        """Return initial SIR vectors.

        Assuming initial infections and recovered are distributed proportional
        to the population.

        Args:
            initial_infectious (float): Initial infectious count.
            initial_recovered (float): Initial recovered count.
        """
        meta_group_pops = np.array([sum(g.pop) for g in self.meta_group_list])
        meta_group_prop = meta_group_pops / sum(meta_group_pops)
        initial_infectious = initial_infectious * meta_group_prop
        initial_recovered = initial_recovered * meta_group_prop

        return self.get_init_SIR_vec(initial_infectious, initial_recovered)


    def get_init_SIR_vec(self, initial_infectious : np.ndarray,
                         initial_recovered : np.ndarray,
                         weight: str = "population"):
        """Return initial SIR vectors.

        Args:
            initial_infectious (np.ndarray): Initial infectious count (per meta group).
            initial_recovered (np.ndarray): Initial recovered count (per meta group).
            weight (str): {population, population x contacts, most_social}
        """
        SIR = []
        for i in range(len(self.meta_group_list)):
            group = self.meta_group_list[i]
            S0, I0, R0 = group.get_init_SIR(initial_infectious[i],
                                            initial_recovered[i],
                                            weight=weight)
            SIR.append([S0, I0, R0])
        SIR = np.array(SIR)

        S0 = self.flatten(SIR[:,0])
        I0 = self.flatten(SIR[:,1])
        R0 = self.flatten(SIR[:,2])
        return S0, I0, R0

    def get_outside_rate(self, outside_rates):
        '''
        Return flattened outside_rate array for inputting into sim().

        Args:
            outside_rates: list of #outside infections/day for each metagroup
        '''
        dim_tot = 0 #keep track of 2d structure within flattened array
        cum_tot = []
        for i in self.meta_group_list:
            cum_tot.append(dim_tot)
            dim_tot += i.K
        outside_rate = np.zeros(dim_tot)
        for i in range(len(self.meta_group_list)):
            for j in range(self.meta_group_list[i].K):
                outside_rate[cum_tot[i]+j] =  outside_rates[i]*self.meta_group_list[i].pop[j]/np.sum(self.meta_group_list[i].pop)
        return np.array(outside_rate)

    def get_num_groups(self):
        ''' Returns the number of groups '''
        return len(self.idx2groupname)

    def metagroup2group(self, v : np.array):
        '''
        Given a list of floats v giving values corresponding to metagroups,
        returns a longer vector that has value v[mg] for all groups in metagroup mg
        '''
        out = np.zeros(self.get_num_groups())
        assert(len(v) == len(self.meta_group_list))

        # The output vector is arranged such that the first meta-group has entries
        # corresponding to its groups first, then the second meta-group han entries
        # corresponding to its groups, etc.

        j = 0 # This will be our output index
        for i in range(len(self.meta_group_list)): # iterate over the meta-groups
            K = self.meta_group_list[i].K # number of groups in the meta-group
            out[j:K+j] = v[i] # write K copies of the number v[i]
            j += K
        assert(j == self.get_num_groups())

        return out