import numpy as np

class sim:
    '''
    Simulates COVID in a population using a SIR-style model.
    Infections last a single generation, after which people recover and are immune for the remainder of the simulation.
    Infections may be discovered immediately (e.g., through symptomatic self-reporting or surveillance) or may be
    discovered later (e.g., because of a 1-time asymptomatic test applied to the whole population, or because the
    student had symptoms due to another cause and got a COVID test).

    Supports multiple groups. Takes as input a matrix (infection_matrix) that determines the number of secondary
    infections created in each group j over 1 generation by a single infection in group i.

    Appropriate choices for this matrix can be determined using the micro class.  This takes as input testing
    interventions as well as social and biological parameters to produce this matrix.

    During execution, S, I, and R are TxK matrices that track the number susceptible, infectious, and recovered in each
    of the K groups over T time periods.
    '''
    def __init__(self, max_T: int, init_susceptible: np.ndarray,
                 init_infected: np.ndarray, init_recovered: np.ndarray,
                 infection_rate: np.ndarray, infection_discovery_frac: float = 1,
                 recovered_discovery_frac: float = 1, generation_time: float = 1):
        """Initialize an SIR-style simulation of COVID spread.

        Group dynamics can be captured by providing vectors of
        initial counts. If there are K groups, the length of init_susceptible,
        init_infected, and init_recovered should all be K.

        To model the discovery of people who
        have or had COVID (I or R), we model a fraction infection_discovery_frac
        of new infections as being discovered in the generation that they start.
        If they are not discovered then, they become "hidden recovered".  A
        fraction recovered_discovery_frac of the undiscovered recovered infections
        are discovered in each generation.

        Args:
            max_T (int): Maximum number of time periods to simulate.
            init_susceptible (np.ndarray): Vector of the initial number of \
                people in each group that are susceptible.
            init_infected (np.ndarray): Vector of the initial number of \
                people in each group that are infected.
            init_recovered (np.ndarray): Vector of the initial number of \
                people in each group that are recovered.
            infection_rate (np.ndarray): Matrix where infection_rate[i,j] is \
                the number of new infections that an infected person in \
                group i creates in group j
            infection_discovery_frac (float): Fraction of infections being \
                discovered in the generation that they start. Defaults to 1.
            recovered_discovery_frac (float): Fraction of recovered being \
                discovered in each generation. Defaults to 1.
            generation_time (float): Length of the generation interval in \
                whatever units the user would like to use, e.g., days or weeks. \
                This is not used, except to support plotting. Defaults to 1.
        """
        assert (max_T > 0)

        self.max_T = max_T # Maximum number of periods we can simulate
        self.t = 0 # current time period
        self.generation_time = generation_time
        assert (infection_discovery_frac >= 0) and (infection_discovery_frac <= 1)
        assert (recovered_discovery_frac >= 0) and (recovered_discovery_frac <= 1)
        self.infection_discovery_frac = infection_discovery_frac
        self.recovered_discovery_frac = recovered_discovery_frac

        self.K = len(init_susceptible) # Number of groups
        assert (len(init_infected) == self.K)
        assert (len(init_recovered) == self.K)

        self.S = np.zeros((self.max_T, self.K)) # susceptible
        self.I = np.zeros((self.max_T, self.K)) # infected
        self.R = np.zeros((self.max_T, self.K)) # recovered
        self.D = np.zeros((self.max_T, self.K)) # discovered infections
        self.H = np.zeros((self.max_T, self.K)) # hidden recovered infections

        self.I[0] = init_infected
        self.S[0] = init_susceptible
        self.R[0] = init_recovered

        self.D[0] = self.I[0]*self.infection_discovery_frac
        self.H[0] = self.I[0]*(1-self.infection_discovery_frac)

        assert ((init_susceptible >= 0).all())
        assert ((init_infected >= 0).all())
        assert ((init_recovered >= 0).all())

        assert (infection_rate.shape == (self.K, self.K))
        self.infection_rate = infection_rate


    def step(self, nsteps: int = 1, infection_rate: np.ndarray = None,
             infection_discovery_frac: float = None,
             recovered_discovery_frac: float = None):
        """Take n steps forward in the simulation.

        As parameters like infection_rate can change over time, this function
        can take those parameters as input.

        Args:
            infection_rate (np.ndarray): Matrix where infection_rate[i,j] is \
                the number of new infections that an infected person in \
                group i creates in group j. Defaults to self.infection_rate.
            infection_discovery_frac (float): Fraction of infections being \
                discovered in the generation that they start. \
                Defaults to self.infection_discovery_frac.
            recovered_discovery_frac (float): Fraction of recovered being \
                discovered in each generation. \
                Defaults to self.recovered_discovery_frac.
        """
        assert(nsteps >= 1)
        # take multiple steps if necessary
        if nsteps > 1:
            for _ in range(nsteps):
                self.step(infection_rate=infection_rate,
                          infection_discovery_frac=infection_discovery_frac,
                          recovered_discovery_frac=recovered_discovery_frac)
            return

        if infection_rate is None:
            infection_rate = self.infection_rate
        if infection_discovery_frac is None:
            infection_discovery_frac = self.infection_discovery_frac
        if recovered_discovery_frac is None:
            recovered_discovery_frac = self.recovered_discovery_frac

        t = self.t

        assert(t+1 < self.max_T) # enforce max generation

        # vector giving the fraction susceptible in each group
        frac_susceptible = self.S[t] / (self.S[t] + self.I[t] + self.R[t])

        # The number of new infections in each group that would result, if everyone were susceptible
        # A = I[t-1] is a vector containing the number of infections in each source group
        # np.matmul(A, infection_rate) has a value at entry j of sum(A[k], infection_rate[k,j])
        self.I[t+1] = np.matmul(self.I[t], infection_rate)

        # Adjust this for the fact that not everyone is susceptible. This is an elementwise product.
        self.I[t+1] = self.I[t+1] * frac_susceptible

        # We can't infect more than the number of susceptible people.
        # np.minimum applied to two arrays returns the elementwise minimum.
        self.I[t+1] = np.minimum(self.I[t+1],self.S[t])

        # Move the infected people out of susceptible and in to recovered
        self.S[t+1] = self.S[t]-self.I[t+1]
        self.R[t+1] = self.R[t]+self.I[t]

        # The old hidden recoveries are either discovered (with probability recovered_discovery_frac) or move forward
        # into the next time period as hidden recoveries
        self.D[t+1] = self.H[t]*recovered_discovery_frac # discovery of old hidden recoveries
        self.H[t+1] = self.H[t]*(1-recovered_discovery_frac)

        # New infections are either discovered immediately (with probability infection_discovery_frac) or
        # become hidden recoveries
        self.D[t+1] = self.D[t+1] + self.I[t+1] * infection_discovery_frac
        self.H[t+1] = self.H[t+1] + self.I[t+1] * (1-infection_discovery_frac)

        self.t = self.t + 1 # Move time forward by one step

        return self

    def __get_metric__(self,metric):
        if metric == 'S':
            return self.S
        elif metric == 'I':
            return self.I
        elif metric == 'R':
            return self.R
        elif metric == 'D':
            return self.D
        elif metric == 'H':
            return self.H
        else:
            raise ValueError('metric argument must be one of S,I,R,D,H')

    def get_metric_for_group(self, metric, group, normalize = False, cumulative = False):
        '''
        Returns a vector where component t contains the number of people of a particular type (S,I,R,D,H)
        in generation t in a specific group. If normalize is true, then this is normalized to the group's population
        size.  For example, to get a vector containing the number of people with an active infection in group 0 and
        each generation, call get_metric_in_group('I').
        '''
        assert(group>=0)
        assert(group<self.K)

        y = self.__get_metric__(metric)[:,group]

        if normalize:
            pop = self.S[:, group] + self.I[:, group] + self.R[:, group]  # total population by time
            y = y / pop

        if cumulative:
            return np.cumsum(y, axis=0)
        else:
            return y


    def get_metric(self, metric, aggregate=True, normalize=False, cumulative=False):
        '''
        Returns a vector where component t contains the total number of people of a particular type (S,I,R,D,H)
        in generation t.
        If aggregate is true, then this is summed across the groups.
        If normalize is true, then this is normalized to the population size.
        When normalize = True and aggregate = False, then the fraction returned is the number of people in the group
        infected divided by the size *of that group*
        For example, to get a vector containing the total number of people that were discovered each week,
        call get_metric('D'). To make it cumulative, call get_metric('D',cumulative=True)
        '''
        y = self.__get_metric__(metric)
        if normalize:
            pop = self.S + self.I + self.R  # total population in each group
            y = y / pop

        if aggregate:
            y = np.sum(y, axis=1)

        if cumulative:
            return np.cumsum(y,axis=0)
        else:
            return y

    def get_infected(self, aggregate=True, normalize=False, cumulative=False):
        return self.get_metric('I', aggregate, normalize, cumulative)

    def get_infected_for_group(self, metric, group, normalize = False, cumulative = False):
        return self.get_metric_for_group('I', group, normalize, cumulative)

    def get_discovered(self, aggregate=True, normalize=False, cumulative=False):
        return self.get_metric('D', aggregate, normalize, cumulative)

    def get_discovered_for_group(self, group, normalize=False, cumulative=False):
        return self.get_metric_for_group('D', group, normalize, cumulative)

    def get_isolated(self, group = False, isolation_len = 2):
        '''
        Returns the number of people in isolation during the generation.
        isolation_len is the number of generations that isolation lasts.
        group is the group to look at.  If group is false, then aggregate across everyone.
        '''
        if group == False:
            discovered = self.get_discovered()
        else:
            discovered = self.get_discovered_for_group(group)

        isolated = np.zeros(self.max_T)
        for t in range(self.max_T):
            for i in range(isolation_len):
                if t-i >= 0:
                    # Add in the people who were discovered i generations ago
                    isolated[t] = isolated[t] + discovered[t-i]

        return isolated

    def get_generation_time(self):
        return self.generation_time


'''
Returns an infection rate matrix that corresponds to well-mixed interactions between groups, where each group has a
amount of contact during their infectious period (summed over all exposed groups) given by the vector
marginal_contacts and a population size by the vector pop.  The number of infections per unit of contact in
marginal_contact (infections_per_contact) is a scalar and applies to all contacts. The units of "marginal_contacts"
could be total contacts over the course of someone's infectious period, in which case infections_per_contact is the
same as the probability of transmission given contact.  It could also be the *rate* at which people have contact per
day, in which case infections_per_contact should be the probability of transmission given contact times the expected
length of a person's infectious period.
'''
def well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact):
    assert (len(pop) == len(marginal_contacts))
    frac_pop = pop / np.sum(pop)  # convert to a fraction
    K = len(marginal_contacts)

    # An incoming contact lands in a population with a probability proportional to its number of outgoing contacts,
    # which is q[i] = frac_pop[i]*marginal_contacts[i].
    q = frac_pop * marginal_contacts / np.sum(frac_pop * marginal_contacts)

    # The total number of people infected by a positive in group i is infections_per_contact * marginal_contacts[i]
    r = infections_per_contact * marginal_contacts

    # When a person in group i is infected, the number of people infected in group j is then r[i]*q[j]

    infection_rates = np.outer(r, q)
    return infection_rates
