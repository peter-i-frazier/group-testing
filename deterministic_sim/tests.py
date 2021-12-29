import numpy as np
import micro
from sim import sim
from groups import meta_group, population
from micro import __days_infectious_perfect_sensitivity__, days_infectious
import matplotlib.pyplot as plt


def is_constant_population_size(s):
    """Checks that the population size of each group is contant across sim."""
    pop = s.get_metric("S") + s.get_metric("I") + s.get_metric("R")
    return all(np.isclose(pop, np.max(pop)))


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

    assert is_constant_population_size(s)

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
    mg = meta_group('Test', pop, marginal_contacts)
    infection_rate = mg.infection_matrix(infections_per_contact)
    generation_time = 4/7 # in units of weeks

    s = sim(T,S0,I0,R0,infection_rate,0,0,generation_time)
    s.step(T-1)

    assert is_constant_population_size(s)

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
    # infections_per_contact = 5/6 * micro.days_infectious(3.5,1) / micro.days_infectious(7,2)
    infections_per_contact = 5/6 * micro.days_infectious(np.inf,1) / micro.days_infectious(7,2)
    generation_time = 4/7 # weeks

    marginal_contacts = np.arange(1,K+1)
    mg = meta_group('Test', pop, marginal_contacts)
    infection_rate = mg.infection_matrix(infections_per_contact)

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

    assert is_constant_population_size(s)

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
    mg = meta_group('Test', pop, marginal_contacts)
    infection_rate = mg.infection_matrix(infections_per_contact)
    generation_time = 4/7 # in units of weeks

    s = sim(T,S0,I0,R0,infection_rate,generation_time)
    s.step(T-1)

    # The group with 0 contacts should not have any infections
    assert(np.isclose(s.get_metric_for_group('I', 0),np.zeros(T)).all())
    assert is_constant_population_size(s)

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
    dec_days_infectious = micro.days_infectious(7,2)

    # We believe that boosters reduce the transmission of virus by a factor of 2
    booster_effectiveness = 0.5


    # A ballpark estimate is that R0 with 1x / wk testing and a 2-day delay from sampling to isolation is 5,
    # with a 4 day generation time.  We think that this outbreak was driven by the people with 6 contacts per period
    # above because the sum of their proportions of the UG population is about 1900 people. To achieve this, we set
    # infections_per_contact = 5/6 so that the number of secondary infections from someone with 6 contacts is 5.
    # We then adjust this by multiplying by the number of days infectious under our testing strategy, divided by the
    # number under our December Omicron outbreak.

    # 1x / week testing and 2 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(7,2) / dec_days_infectious
    infection_rate = meta_group('UG', pop, marginal_contacts).infection_matrix(infections_per_contact)
    s = sim(T, S0, I0, R0, infection_rate, generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='1x/wk, 2d delay')
    plt.subplot(212)
    plt.plot(np.arange(T)*generation_time, s.get_isolated(), label='1x/wk, 2d delay')

    # 2x / week testing and 1.5 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(3.5,1.5) / dec_days_infectious
    infection_rate = meta_group('UG', pop, marginal_contacts).infection_matrix(infections_per_contact)
    s = sim(T, S0, I0, R0, infection_rate, generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='2x/wk, 1.5d delay')
    plt.subplot(212)
    plt.plot(np.arange(T)*generation_time, s.get_isolated(), label='2x/wk, 1.5d delay')

    # 2x / week testing and 2 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(3.5,2) / dec_days_infectious
    infection_rate = meta_group('UG', pop, marginal_contacts).infection_matrix(infections_per_contact)
    s = sim(T, S0, I0, R0, infection_rate, generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='2x/wk, 2d delay')
    plt.subplot(212)
    plt.plot(np.arange(T)*generation_time, s.get_isolated(), label='2x/wk, 2d delay')

    # 2x / week testing and 1 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(3.5,1) / dec_days_infectious
    infection_rate = meta_group('UG', pop, marginal_contacts).infection_matrix(infections_per_contact)
    s = sim(T, S0, I0, R0, infection_rate, generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='2x/wk, 1d delay')
    plt.subplot(212)
    plt.plot(np.arange(T)*generation_time, s.get_isolated(), label='2x/wk, 1d delay')

    # 1x / week testing and 1 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(7,1) / dec_days_infectious
    infection_rate = meta_group('UG', pop, marginal_contacts).infection_matrix(infections_per_contact)
    s = sim(T,S0,I0,R0,infection_rate,generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='1x/wk, 1d delay')
    plt.subplot(212)
    plt.plot(np.arange(T) * generation_time, s.get_isolated(), label='1x/wk, 1d delay')

    # 1x / week testing and 1.5 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(7,1.5) / dec_days_infectious
    infection_rate = meta_group('UG', pop, marginal_contacts).infection_matrix(infections_per_contact)
    s = sim(T,S0,I0,R0,infection_rate,generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='1x/wk, 1.5d delay')
    plt.subplot(212)
    plt.plot(np.arange(T) * generation_time, s.get_isolated(), label='1x/wk, 1.5d delay')

    # No surveillance
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(np.inf,1) / dec_days_infectious
    infection_rate = meta_group('UG', pop, marginal_contacts).infection_matrix(infections_per_contact)
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


def test_days_infectious():

    # If there is infinite time between tests, the time infectious in the presence of testing should be whatever
    # the length of the maximum infectious period is. This should be true regardless of the sensitivity.
    assert __days_infectious_perfect_sensitivity__(np.inf,1,max_infectious_days = 5) == 5
    assert days_infectious(np.inf,1,max_infectious_days = 5) == 5

    # Time infectious with testing should be less than the maximum
    assert days_infectious(2,1,max_infectious_days = 5) < 5

    # Testing with more delay should result in a longer time infectious
    assert days_infectious(7,1) > days_infectious(2, 1)

    # A shorter max infectioun time should result in fewer days infectious
    assert days_infectious(5, 1, max_infectious_days=5) < days_infectious(5, 1, max_infectious_days=7)

    # Codes should agree when sensitivity is perfect
    assert days_infectious(7,1,sensitivity=1) == __days_infectious_perfect_sensitivity__(7,1)
    assert days_infectious(4,2,sensitivity=1) == __days_infectious_perfect_sensitivity__(4,2)

    # A lower sensitivity should result in more days infectious
    assert days_infectious(5, 1, sensitivity=.5) > days_infectious(5, 1, sensitivity=.7)

    # A sensitivity of 0 should be like not testing
    assert days_infectious(2,1,max_infectious_days = 5, sensitivity = 0) == 5

    # days infectious should be larger when the days between tests is larger
    for d in range(20):
        assert days_infectious(d+2, 1) > days_infectious(d+1, 1)

    days_between_tests = np.arange(1, 20, .1)
    y70 = [days_infectious(d, 1, sensitivity = 0.7) for d in days_between_tests]
    y90 = [days_infectious(d, 1, sensitivity = 0.9) for d in days_between_tests]
    y100 = [days_infectious(d, 1, sensitivity = 1) for d in days_between_tests]
    plt.plot(days_between_tests, y70, label = 'Sensitivity = 70%')
    plt.plot(days_between_tests, y90, label = 'Sensitivity = 90%')
    plt.plot(days_between_tests, y100, label = 'Sensitivity = 100%')
    plt.legend()
    plt.xlabel('Number of days between tests')
    plt.ylabel('Expected Days infectious')
    plt.savefig('test_days_infectious1.png', facecolor='w')
    plt.close()

    # Compare what we are using here to what we used before
    days_between_tests = np.arange(1, 20, .1)
    y_old = [days_infectious(d, 1, sensitivity=1, max_infectious_days=7) for d in days_between_tests]
    y_new = [days_infectious(d, 1) for d in days_between_tests]
    plt.plot(days_between_tests, y_old, 'k--', label='Sensitivity = 100%, max_infectious = 7 days, 1 day delay')
    plt.plot(days_between_tests, y_new, 'k-', label='Sensitivity = 60%, max_infectious = default, 1 day delay')

    y_old = [days_infectious(d, 2, sensitivity=1, max_infectious_days=7) for d in days_between_tests]
    y_new = [days_infectious(d, 2) for d in days_between_tests]
    plt.plot(days_between_tests, y_old, 'b--', label='Sensitivity = 100%, max_infectious = 7 days, 2 day delay')
    plt.plot(days_between_tests, y_new, 'b-', label='Sensitivity = 60%, max_infectious = default, 2 day delay')

    plt.legend()
    plt.xlabel('Number of days between tests')
    plt.ylabel('Expected Days infectious')
    plt.savefig('test_days_infectious2.png', facecolor='w')
    plt.close()

def test_sim6():
    '''order of meta-groups: UG, GR, VM/GM/LA, fac/staff'''
    total_pops = [16000, 8000, 3000, 10000]
    T = 20 #num generations

    # Based on Xiangyu's adjusted moment match code, this is the fraction of the population (UG-only?  or all students)
    # broken out by the amount of contact that they have, starting from low contact to high contact.
    pop_fracs = np.array([
            [0.3771077769,
            0.1744692313,
            0.2130738259,
            0.07227578963,
            0.05147885598,
            0.04364552029,
            0.02718384893,
            0.01451678525,
            0.02624836593],
            [0.5790258783,
            0.1683859957,
            0.05712348893,
            0.04359733917,
            0.03548852695,
            0.0451325654,
            0.03935404056,
            0.01167552257,
            0.01055548401,
            0.009661158416],
            [0.472724139,
            0.1443461065,
            0.03264544331,
            0.1245769902,
            0.04056260275,
            0.03439035051,
            0.04498079947,
            0.01334486452,
            0.04825874066,
            0.04416996305],
            [0.5559842653,
            0.1914068446,
            0.09034172148,
            0.07325926546,
            0.03157069457,
            0.02081854892,
            0.01296646907,
            0.01615689928,
            0,
            0.005729729349,
            0.00176556195]])
    pops = total_pops[:]
    for i in range(len(total_pops)):
        pops[i] = total_pops[i]*np.array(pop_fracs[i])


    marginal_contacts = np.array([np.arange(1,len(pops[0])+1),
                    np.arange(1,len(pops[1])+1),
                    np.arange(1,len(pops[2])+1),
                    np.arange(1,len(pops[3])+1)])

    # flatten doesn't work when rows have different lengths :(
    tmp = []
    for i in range(len(marginal_contacts)):
        tmp+=list(marginal_contacts[i]) 
    marginal_contacts_flat = np.array(tmp)

    tmp = []
    for i in range(len(pops)):
        tmp+=list(pops[i])
    pops_flat = np.array(tmp)
    # There were roughly 1900 UG infected during the Omicron outbreak in December 2021.
    # Assume that another 2000 will be infected during winter break
    infected_before_semester = 1900 + 2000

    # Assume 100 active infections to start the semester
    initial_infections = 100

    # Assume a group's previous and new infections are divided proportionally to the amount of contact it has as a
    # group. This is its contact rate * population size
    b = marginal_contacts_flat * pops_flat
    b =  b / np.sum(b)
    R0 = infected_before_semester * b
    I0 = initial_infections * b

    S0 = np.maximum(pops_flat - R0 - I0, 0) # Need to take the max with 0 because R0[i]+S0[i] can be bigger than pop[i]

    generation_time = 4/7 # weeks
    symptomatic_rate = .3 # used for the fraction of new infections will be discovered without surevillance

    # We calibrate our probability of infection to the December Omicron outbreak
    # We use a ballpark estimate of the effective R0 during that period, and an estimate of which group was the most
    # important one in driving that outbreak.  We then assume that the probability of infection is such that this
    # group's effective R0 (under the testing intervention at the time) was equal to the ballpark estimate.
    # Assuming a generation time of 4 days, and 50% day-over-day growth, we get a ballpark estimate of 1.5^4 = 5.06
    # for the december effective R0
    dec_effective_R0 = 5
    dec_contacts_of_key_group = 6
    dec_infections_per_contact = dec_effective_R0 / dec_contacts_of_key_group
    dec_days_infectious = micro.days_infectious(7,2)

    # We believe that boosters reduce the transmission of virus by a factor of 2
    booster_effectiveness = 0.5

    meta_matrix = np.array([[0.8351648352, 0.03296703297, 0.05054945055, 0.08131868132],
                [0.02033898305, 0.9593220339, 0.02033898305, 0],
                [0.1634615385, 0.1057692308, 0.6923076923, 0.03846153846],
                [0.06875, 0.0125, 0.05, 0.86875]])
    # A ballpark estimate is that R0 with 1x / wk testing and a 2-day delay from sampling to isolation is 5,
    # with a 4 day generation time.  We think that this outbreak was driven by the people with 6 contacts per period
    # above because the sum of their proportions of the UG population is about 1900 people. To achieve this, we set
    # infections_per_contact = 5/6 so that the number of secondary infections from someone with 6 contacts is 5.
    # We then adjust this by multiplying by the number of days infectious under our testing strategy, divided by the
    # number under our December Omicron outbreak.

    # 1x / week testing and 2 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(7,2) / dec_days_infectious
    ug = meta_group('UG', pops[0], marginal_contacts[0])
    gr = meta_group('GR', pops[1], marginal_contacts[1])
    pro = meta_group('PR', pops[2], marginal_contacts[2])
    fac = meta_group('FS', pops[3], marginal_contacts[3])
    popul = population([ug, gr, pro, fac], meta_matrix)
    # print(popul.infection_matrix(infections_per_contact))
    s = sim(T, S0, I0, R0, popul.infection_matrix(infections_per_contact), generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='1x/wk, 2d delay')
    plt.subplot(212)
    plt.plot(np.arange(T)*generation_time, s.get_isolated(), label='1x/wk, 2d delay')

        # 2x / week testing and 1.5 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(3.5,1.5) / dec_days_infectious
    ug = meta_group('UG', pops[0], marginal_contacts[0])
    gr = meta_group('GR', pops[1], marginal_contacts[1])
    pro = meta_group('PR', pops[2], marginal_contacts[2])
    fac = meta_group('FS', pops[3], marginal_contacts[3])
    popul = population([ug, gr, pro, fac], meta_matrix)
    # print(popul.infection_matrix(infections_per_contact))
    s = sim(T, S0, I0, R0, popul.infection_matrix(infections_per_contact), generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='2x/wk, 1.5d delay')
    plt.subplot(212)
    plt.plot(np.arange(T)*generation_time, s.get_isolated(), label='2x/wk, 1.5d delay')

    # 2x / week testing and 2 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(3.5,2) / dec_days_infectious
    ug = meta_group('UG', pops[0], marginal_contacts[0])
    gr = meta_group('GR', pops[1], marginal_contacts[1])
    pro = meta_group('PR', pops[2], marginal_contacts[2])
    fac = meta_group('FS', pops[3], marginal_contacts[3])
    popul = population([ug, gr, pro, fac], meta_matrix)
    # print(popul.infection_matrix(infections_per_contact))
    s = sim(T, S0, I0, R0, popul.infection_matrix(infections_per_contact), generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='2x/wk, 2d delay')
    plt.subplot(212)
    plt.plot(np.arange(T)*generation_time, s.get_isolated(), label='2x/wk, 2d delay')

    # 2x / week testing and 1 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(3.5,1) / dec_days_infectious
    ug = meta_group('UG', pops[0], marginal_contacts[0])
    gr = meta_group('GR', pops[1], marginal_contacts[1])
    pro = meta_group('PR', pops[2], marginal_contacts[2])
    fac = meta_group('FS', pops[3], marginal_contacts[3])
    popul = population([ug, gr, pro, fac], meta_matrix)
    # print(popul.infection_matrix(infections_per_contact))
    s = sim(T, S0, I0, R0, popul.infection_matrix(infections_per_contact), generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='2x/wk, 1d delay')
    plt.subplot(212)
    plt.plot(np.arange(T)*generation_time, s.get_isolated(), label='2x/wk, 1d delay')

    # 1x / week testing and 1 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(7,1) / dec_days_infectious
    ug = meta_group('UG', pops[0], marginal_contacts[0])
    gr = meta_group('GR', pops[1], marginal_contacts[1])
    pro = meta_group('PR', pops[2], marginal_contacts[2])
    fac = meta_group('FS', pops[3], marginal_contacts[3])
    popul = population([ug, gr, pro, fac], meta_matrix)
    # print(popul.infection_matrix(infections_per_contact))
    s = sim(T, S0, I0, R0, popul.infection_matrix(infections_per_contact), generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='1x/wk, 1d delay')
    plt.subplot(212)
    plt.plot(np.arange(T) * generation_time, s.get_isolated(), label='1x/wk, 1d delay')

    # 1x / week testing and 1.5 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(7,1.5) / dec_days_infectious
    ug = meta_group('UG', pops[0], marginal_contacts[0])
    gr = meta_group('GR', pops[1], marginal_contacts[1])
    pro = meta_group('PR', pops[2], marginal_contacts[2])
    fac = meta_group('FS', pops[3], marginal_contacts[3])
    popul = population([ug, gr, pro, fac], meta_matrix)
    # print(popul.infection_matrix(infections_per_contact))
    s = sim(T, S0, I0, R0, popul.infection_matrix(infections_per_contact), generation_time)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label='1x/wk, 1.5d delay')
    plt.subplot(212)
    plt.plot(np.arange(T) * generation_time, s.get_isolated(), label='1x/wk, 1.5d delay')

    # No surveillance
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(np.inf,1) / dec_days_infectious
    ug = meta_group('UG', pops[0], marginal_contacts[0])
    gr = meta_group('GR', pops[1], marginal_contacts[1])
    pro = meta_group('PR', pops[2], marginal_contacts[2])
    fac = meta_group('FS', pops[3], marginal_contacts[3])
    popul = population([ug, gr, pro, fac], meta_matrix)
    infection_discovery_frac = symptomatic_rate # 30% are asymptomatic
    recovered_discovery_frac = .01 # 1% of the population is tested for any reason in a given generation
    s = sim(T,S0,I0,R0,popul.infection_matrix(infections_per_contact),generation_time,infection_discovery_frac,recovered_discovery_frac)
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
    plt.ylabel('Total Infected')

    plt.subplot(212)
    plt.rcParams.update({'font.size': 8})
    plt.legend()
    plt.xlabel('Weeks')
    plt.ylabel('Total in Isolation')

    plt.savefig('test_sim6.png', facecolor='w')
    plt.close()
