import numpy as np
import pandas as pd
import micro
from sim import sim, well_mixed_infection_rate
import matplotlib.pyplot as plt

def test_sim_calibration():
    total_pop = 16000 # total UG population
    K = 12 # number of distinct contact groups
    T = 5 # num generations
    initial_group = 11

    # Based on Xiangyu's adjusted moment match code, this is the fraction of the
    # population (UG-only?  or all students) broken out by the amount of contact
    # that they have, starting from low contact to high contact.
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

    # Assume no recovered and 6 infected with Omicron
    R0 = np.zeros(K)
    I0 = np.zeros(K)
    I0[initial_group] = 6

    # Need to take the max with 0 because R0[i]+S0[i] can be bigger than pop[i]
    S0 = np.maximum(pop - R0 - I0, 0)

    generation_time = 4 # days
    symptomatic_rate = .3 # fraction of new infections will be discovered without surevillance

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

    # [12/1 to 12/9] 1x / week testing with 36hr delay
    infections_per_contact = dec_infections_per_contact * micro.days_infectious(7,1.5) / dec_days_infectious
    infection_rate = well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact)
    s = sim(T, S0, I0, R0, infection_rate=infection_rate, generation_time=generation_time)
    s.step(2)

    # [12/10 to 12/16] 1x / week testing with 3 day delay
    infections_per_contact = dec_infections_per_contact * micro.days_infectious(7,3) / dec_days_infectious
    infection_rate = well_mixed_infection_rate(pop, 0.5*marginal_contacts, infections_per_contact)
    s.step(1, infection_rate=infection_rate)

    infections_per_contact = dec_infections_per_contact * micro.days_infectious(7,3) / dec_days_infectious
    infection_rate = well_mixed_infection_rate(pop, 0.33*marginal_contacts, infections_per_contact)
    s.step(1, infection_rate=infection_rate)

    I = s.get_metric('I', aggregate=True, cumulative=True)
    plt.plot(np.arange(T)*generation_time, I, label="model I")

    # Plot actual counts
    dec_daily_positives = list(pd.read_csv("data/dec_infections.csv")['positives'])
    dec_positives = np.cumsum(dec_daily_positives)
    plt.plot(np.arange(15), dec_positives[:15], label="actual")

    plt.title(f"Dec Effective R0 = {dec_effective_R0}, Symptomatic Rate = {symptomatic_rate}")
    plt.rcParams.update({'font.size': 8})
    plt.xlabel('Days Since Dec. 1')
    plt.ylabel('Cumulative Infected')
    plt.legend()
    plt.savefig('test_sim6.png', facecolor='w')
    plt.close()