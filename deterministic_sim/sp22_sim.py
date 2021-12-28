import numpy as np
import pandas as pd
import micro
from sim import sim, well_mixed_infection_rate
import matplotlib.pyplot as plt


def main():
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
    # Assume that another 175 will be infected during winter break
    # According to @Brian Liu's analysis
    infected_before_semester = 1900 + 175
    # infected_before_semester = 1900 + 2000

    # Assume 25 active infections to start the semester
    # According to @Brian Liu's analysis
    initial_infections = 25
    # initial_infections = 100

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
    dec_effective_R0 = 7
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

    def sim_test_regime(tests_per_week, delay):
        """Simulate a testing regime"""
        days_between_tests = 7 / tests_per_week
        infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(days_between_tests,delay) / dec_days_infectious
        infection_rate = well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact)
        s = sim(T, S0, I0, R0, infection_rate=infection_rate, generation_time=generation_time)
        s.step(4)
        infection_rate = well_mixed_infection_rate(pop, 0.85*marginal_contacts, infections_per_contact)
        s.step(T-1-4, infection_rate=infection_rate)

        label = "%dx/wk, %.1fd delay" % (tests_per_week, delay)
        plt.subplot(211)
        plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label=label)
        plt.subplot(212)
        plt.plot(np.arange(T)*generation_time, s.get_isolated(), label=label)

    sim_test_regime(1,2)
    sim_test_regime(1,1.5)
    sim_test_regime(1,1)
    sim_test_regime(2,2)
    sim_test_regime(2,1.5)
    sim_test_regime(2,1)

    # 2x / week testing and 1 day delay
    infections_per_contact = booster_effectiveness * dec_infections_per_contact * micro.days_infectious(3.5,1) / dec_days_infectious
    infection_rate = well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact)
    infection_discovery_frac = symptomatic_rate # 30% are asymptomatic
    recovered_discovery_frac = .01 # 1% of the population is tested for any reason in a given generation
    s = sim(T, S0, I0, R0, infection_rate=infection_rate, generation_time=generation_time,
            infection_discovery_frac=infection_discovery_frac,
            recovered_discovery_frac=recovered_discovery_frac)
    s.step(4)
    infection_rate = well_mixed_infection_rate(pop, 0.85*marginal_contacts, infections_per_contact)
    s.step(T-1-4, infection_rate=infection_rate)

    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), \
             color = 'powderblue', label='2x/wk, 1d delay')
    plt.subplot(212)
    plt.plot(np.arange(T)*generation_time, s.get_isolated(), color = 'powderblue', label='2x/wk, 1d delay')




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


if __name__ == "main":
    main()
