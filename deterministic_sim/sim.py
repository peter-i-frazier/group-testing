import numpy as np
import micro
import matplotlib.pyplot as plt
import pytest

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
    def __init__(self, max_T, init_susceptible, init_infected, init_recovered, infection_rate, generation_time,
                 infection_discovery_frac=1, recovered_discovery_frac=1):
        '''
        max_T is the maximum number of time periods we will simulate
        init_susceptible, init_infected, init_recovered are vectors containing the number of people in each group that are susceptible, infectious, and recovered.
        infection_rate is a matrix where infection_rate[i,j] is the number of new infections that an infected person in group i creates in group j
        Letting K be the number of groups, the three init_ vectors should be length K.  infection_rate should be KxK.

        generation_time is the length of the generation interval in whatever units the user would like to use, e.g.,
        days or weeks. This is not used, except to support plotting.

        To model the discovery of people who have or had COVID (people in I or R), we model a fraction
        infection_discovery_frac of new infections as being discovered in the generation that they start.
        If they are not discovered then, they become "hidden recovered".  A fraction recovered_discovery_frac of the
        undiscovered recovered infections are discovered in each generation. By default, everything is discovered,
        making these two fractions 1.
        '''
        assert (max_T > 0)

        self.max_T = max_T # Maximum number of periods we can simulate
        self.t = 0 # current time period
        self.generation_time = generation_time
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

    def step(self, nsteps = 1):
        """Move the simulation forward by nsteps"""

        assert(nsteps >= 1)
        if nsteps > 1: # If the user is asking for more than 1 step, run this function the requested number of times
            for i in range(nsteps):
                self.step() # Simulate forward 1 step
            return

        t = self.t

        assert(t+1 < self.max_T) # Otherwise, we will run out of space in the matrices S, I, R when we try to add at t+1

        # vector giving the fraction susceptible in each group
        frac_susceptible = self.S[t] / (self.S[t] + self.I[t] + self.R[t])

        # The number of new infections in each group that would result, if everyone were susceptible
        # A = I[t-1] is a vector containing the number of infections in each source group
        # np.matmul(A, infection_rate) has a value at entry j of sum(A[k], infection_rate[k,j])
        self.I[t+1] = np.matmul(self.I[t], self.infection_rate)

        # Adjust this for the fact that not everyone is susceptible. This is an elementwise product.
        self.I[t+1] = self.I[t+1] * frac_susceptible

        # We can't infect more than the number of susceptible people.
        # np.minimum applied to two arrays returns the elementwise minimum.
        self.I[t+1] = np.minimum(self.I[t+1],self.S[t])

        # Move the infected people out of susceptible and in to recovered
        self.S[t+1] = self.S[t]-self.I[t+1]
        self.R[t+1] = self.R[t]+self.I[t+1] # XXX is this a bug?  I think I[t+1] here should be t.  Write a test case that ensures that the

        # The old hidden recoveries are either discovered (with probability recovered_discovery_frac) or move forward
        # into the next time period as hidden recoveries
        self.D[t+1] = self.H[t]*self.recovered_discovery_frac # discovery of old hidden recoveries
        self.H[t+1] = self.H[t]*(1-self.recovered_discovery_frac)

        # New infections are either discovered immediately (with probability infection_discovery_frac) or
        # become hidden recoveries
        self.D[t+1] = self.D[t+1] + self.I[t+1] * self.infection_discovery_frac
        self.H[t+1] = self.H[t+1] + self.I[t+1] * (1-self.infection_discovery_frac)

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


def test_sim1():
    K = 3
    T = 20

    # The populations are symmetric
    pop = 100 * np.ones(3)/3
    R0 = 100 * np.array([.1, .1, .1])
    I0 = 100 * np.array([.01, .01, .01])
    S0 = pop - R0 - I0

    # Here, each group has an overall rate at which they infect others.
    # These infections are distributed proportionally to the overall population.
    # infection_rate[i,j] is the number of new infections that an infected person in i creates in j.
    contact_rates = np.array([0, 1, 2])
    infection_rate = np.outer(contact_rates,pop/100)

    generation_time = 4/7 # in units of weeks

    s = sim(T,S0,I0,R0,infection_rate,generation_time)
    s.step(T-1)

    y = s.get_metric('I', aggregate=False, normalize=True)
    for i in range(K):
        plt.plot(np.arange(T), y[:,i], label=i)
    plt.legend(title='Group')
    plt.xlabel('Generation')
    plt.ylabel('Fraction Infected')
    plt.savefig('test_sim1.png', facecolor='w')
    plt.close()

    # Because the S,I,R amounts for each population were initially symmetric and the new infections are distributed
    # proportionally across the populations, the number infected should stay symmetric over time.
    # That is, for each value of t, y[t,0], y[t,1], and y[t,2] should be nearly identical.
    assert(np.isclose(y[:,0],y[:,1]).all())
    assert(np.isclose(y[:,0],y[:,2]).all())

def test_sim4():
    K = 3
    T = 20

    # The populations are symmetric
    pop = 100 * np.ones(3)/3
    R0 = 100 * np.array([.1, .1, .1])
    I0 = 100 * np.array([0, .01, .01])
    S0 = pop - R0 - I0

    infections_per_contact = 1
    marginal_contacts = np.array([0,1,2])
    infection_rate = well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact)
    generation_time = 4/7 # in units of weeks

    s = sim(T,S0,I0,R0,infection_rate,generation_time,0,0)
    s.step(T-1)

    # Since no one is discovered (infection_discovery_frac and recovered_discovery_frac are 0 above),
    # the number discovered should be 0
    discovered = s.get_discovered(cumulative=True)
    infected = s.get_infected(cumulative=True)

    plt.plot(np.arange(T), infected, label='Infected')
    plt.plot(np.arange(T), discovered, label='Discovered')
    plt.legend()
    plt.xlabel('Generation')
    plt.ylabel('Fraction Infected')
    plt.savefig('test_sim4.png', facecolor='w')
    plt.close()

    assert(np.isclose(discovered, np.zeros(T)).all())

def test_sim3():
    total_pop = 16000 #total UG population
    K = 12 #number of distinct contact groups
    T = 20 #num generations

    # Based on Xiangyu's adjusted moment match code, this is the fraction of the population (UG-only?  or all students)
    # broken out by the amount of contact that they have, starting from low contact to high contact.
    pop_frac = np.array(
        [0.46666204134859246,
         0.21110377424326393,
         0.13427835918082992,
         0.040993148549700376,
         0.05561449714374835,
         0.03772148571886847,
         0.020557396664413034,
         0.01829685176380435,
         0.003308325172067398,
         0.006056046991853723,
         0.0027991704708900224,
         0.002608902751968122])
    pop = total_pop * pop_frac

    # Calculate the probability of infection.
    # A ballpark estimate is that R0 with 1x / wk testing and a 2-day delay from sampling to isolation is 5,
    # with a 4 day generation time.  We think that this outbreak was driven by the people with 6 contacts per period
    # above because the sum of their proportions of the UG population is about 1900 people. To achieve this, we set
    # infections_per_contact = 5/6 so that the number of secondary infections from someone with 6 contacts is 5.
    # We then adjust this by multiplying by the number of days infectious under our testing strategy, divided by the
    # number under our December Omicron outbreak.
    # infections_per_contact = 5/6 * micro.days_infectious(3.5,7,1) / micro.days_infectious(7,7,2)
    infections_per_contact = 5/6 * micro.days_infectious(np.inf,7,1) / micro.days_infectious(7,7,2)
    generation_time = 4/7 # weeks

    marginal_contacts = np.arange(1,K+1)
    infection_rate = well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact)

    # There were roughly 1900 UG infected during the Omicron outbreak in December 2021.
    # Assume that another 2000 will be infected during winter break
    infected_before_semester = 1900 + 2000

    # Assume 100 active infections to start the semester
    initial_infections = 100

    # Assume a group's previous and new infections are divided proportionally to the amount of contact it has as a
    # group. This is its contact rate * population size
    b = marginal_contacts * pop_frac
    b =  b / np.sum(b)
    R0 = infected_before_semester * b
    I0 = initial_infections * b

    S0 = np.maximum(pop - R0 - I0, 0) # Need to take the max with 0 because R0[i]+S0[i] can be bigger than pop[i]

    s = sim(T,S0,I0,R0,infection_rate,generation_time)
    s.step(T-1)

    infected = s.get_metric('I', aggregate=True)
    cum_infected = s.get_metric('I', aggregate=True,cumulative=True)
    plt.plot(np.arange(T)*generation_time, infected, label='New Infections')
    plt.plot(np.arange(T)*generation_time, cum_infected, label='Cumulative Infections')

    plt.legend()
    plt.xlabel('Weeks')
    plt.ylabel('Number Infected')
    plt.savefig('test_sim3.png', facecolor='w')
    plt.close()

def test_sim2():
    K = 3
    T = 20

    # The populations are symmetric
    pop = 100 * np.ones(3)/3
    R0 = 100 * np.array([.1, .1, .1])
    I0 = 100 * np.array([0, .01, .01])
    S0 = pop - R0 - I0

    infections_per_contact = 1
    marginal_contacts = np.array([0,1,2])
    infection_rate = well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact)
    generation_time = 4/7 # in units of weeks

    s = sim(T,S0,I0,R0,infection_rate,generation_time)
    s.step(T-1)

    # The group with 0 contacts should not have any infections
    assert(np.isclose(s.get_metric_for_group('I', 0),np.zeros(T)).all())

    for i in range(K):
        plt.plot(np.arange(T), s.get_metric_for_group('I', i,normalize=True), label=i)
    plt.legend(title='Group')
    plt.xlabel('Generation')
    plt.ylabel('Fraction Infected')
    plt.savefig('test_sim2.png', facecolor='w')
    plt.close()

def test_sim5():
    total_pop = 16000 #total UG population
    K = 12 #number of distinct contact groups
    T = 20 #num generations

    # Based on Xiangyu's adjusted moment match code, this is the fraction of the population (UG-only?  or all students)
    # broken out by the amount of contact that they have, starting from low contact to high contact.
    pop_frac = np.array(
        [0.46666204134859246,
         0.21110377424326393,
         0.13427835918082992,
         0.040993148549700376,
         0.05561449714374835,
         0.03772148571886847,
         0.020557396664413034,
         0.01829685176380435,
         0.003308325172067398,
         0.006056046991853723,
         0.0027991704708900224,
         0.002608902751968122])
    pop = total_pop * pop_frac

    marginal_contacts = np.arange(1,K+1)

    # There were roughly 1900 UG infected during the Omicron outbreak in December 2021.
    # Assume that another 2000 will be infected during winter break
    infected_before_semester = 1900 + 2000

    # Assume 100 active infections to start the semester
    initial_infections = 100

    # Assume a group's previous and new infections are divided proportionally to the amount of contact it has as a
    # group. This is its contact rate * population size
    b = marginal_contacts * pop_frac
    b =  b / np.sum(b)
    R0 = infected_before_semester * b
    I0 = initial_infections * b

    S0 = np.maximum(pop - R0 - I0, 0) # Need to take the max with 0 because R0[i]+S0[i] can be bigger than pop[i]

    generation_time = 4/7 # weeks
    symptomatic_rate = .3 # used for the fraction of new infections will be discovered without surevillance

    # We calibrate our probability of infection to the December Omicron outbreak
    # We use a ballpark estimate of the effective R0 during that period, and an estimate of which group was the most
    # important one in driving that outbreak.  We then assume that the probability of infection is such that this
    # group's effective R0 (under the testing intervention at the time) was equal to the ballpark estimate.
    # Assuming a generation time of 4 days, and 50% day-over-day growth, we get a ballpark estimate of 1.5^4 = 5.06
    # for the december effective R0
    dec_effective_R0 = 10
    dec_contacts_of_key_group = 6
    dec_infections_per_contact = dec_effective_R0 / dec_contacts_of_key_group
    dec_days_infectious = micro.days_infectious(7,7,2)

    # We believe that boosters reduce the transmission of virus by a factor of 2
    booster_effectiveness = 0.5


    # A ballpark estimate is that R0 with 1x / wk testing and a 2-day delay from sampling to isolation is 5,
    # with a 4 day generation time.  We think that this outbreak was driven by the people with 6 contacts per period
    # above because the sum of their proportions of the UG population is about 1900 people. To achieve this, we set
    # infections_per_contact = 5/6 so that the number of secondary infections from someone with 6 contacts is 5.
    # We then adjust this by multiplying by the number of days infectious under our testing strategy, divided by the
    # number under our December Omicron outbreak.

    # 1x / week testing and 2 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(7,7,2) / dec_days_infectious
    infection_rate = well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact)
    s = sim(T, S0, I0, R0, infection_rate, generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='1x/wk, 2d delay')
    plt.subplot(212)
    plt.plot(np.arange(T)*generation_time, s.get_isolated(), label='1x/wk, 2d delay')

    # 2x / week testing and 1.5 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(3.5,7,1.5) / dec_days_infectious
    infection_rate = well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact)
    s = sim(T, S0, I0, R0, infection_rate, generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='2x/wk, 1.5d delay')
    plt.subplot(212)
    plt.plot(np.arange(T)*generation_time, s.get_isolated(), label='2x/wk, 1.5d delay')

    # 2x / week testing and 2 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(3.5,7,2) / dec_days_infectious
    infection_rate = well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact)
    s = sim(T, S0, I0, R0, infection_rate, generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='2x/wk, 2d delay')
    plt.subplot(212)
    plt.plot(np.arange(T)*generation_time, s.get_isolated(), label='2x/wk, 2d delay')

    # 2x / week testing and 1 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(3.5,7,1) / dec_days_infectious
    infection_rate = well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact)
    s = sim(T, S0, I0, R0, infection_rate, generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='2x/wk, 1d delay')
    plt.subplot(212)
    plt.plot(np.arange(T)*generation_time, s.get_isolated(), label='2x/wk, 1d delay')

    # 1x / week testing and 1 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(7,7,1) / dec_days_infectious
    infection_rate = well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact)
    s = sim(T,S0,I0,R0,infection_rate,generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='1x/wk, 1d delay')
    plt.subplot(212)
    plt.plot(np.arange(T) * generation_time, s.get_isolated(), label='1x/wk, 1d delay')

    # 1x / week testing and 1.5 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(7,7,1.5) / dec_days_infectious
    infection_rate = well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact)
    s = sim(T,S0,I0,R0,infection_rate,generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='1x/wk, 1.5d delay')
    plt.subplot(212)
    plt.plot(np.arange(T) * generation_time, s.get_isolated(), label='1x/wk, 1.5d delay')

    # No surveillance
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(np.inf,7,1) / dec_days_infectious
    infection_rate = well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact)
    infection_discovery_frac = symptomatic_rate # 30% are asymptomatic
    recovered_discovery_frac = .01 # 1% of the population is tested for any reason in a given generation
    s = sim(T,S0,I0,R0,infection_rate,generation_time,infection_discovery_frac,recovered_discovery_frac)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), 'k-', label='No surveillance, Discovered')
    plt.plot(np.arange(T)*generation_time, s.get_infected(aggregate=True,cumulative=True), 'k--', label='No surveillance, Infected')
    plt.subplot(212)
    plt.plot(np.arange(T) * generation_time, s.get_isolated(), 'k', label='No surveillance')


    plt.subplot(211)
    plt.title('Dec Effective R0 = {}, Symptomatic Rate = {}'.format(dec_effective_R0, symptomatic_rate))
    plt.rcParams.update({'font.size': 8})
    plt.legend()
    #plt.xlabel('Weeks')
    plt.ylabel('UG Infected')

    plt.subplot(212)
    plt.rcParams.update({'font.size': 8})
    plt.legend()
    plt.xlabel('Weeks')
    plt.ylabel('UG in Isolation')

    plt.savefig('test_sim5.png', facecolor='w')
    plt.close()
