import numpy as np
import pandas as pd
import micro
from sim import sim
from groups import meta_group
import matplotlib.pyplot as plt


TOTAL_POP = 16000       # total UG population
K = 12                  # number of distinct contact groups
T = 5                   # num generations
INITIAL_GROUP = 11      # index of group containing initial infections
INITIAL = 6             # number of initial infections
SYMPTOMATIC_RATE = .3   # fraction of new infections will be discovered
                        # without surevillance

# We aim to calibrate the effective R0 and generation_time of the simulation
# to the observed December 2021 outbreak. We assume that the most social of
# the population of gets infected first. Based on the 1900 observed infections
# during this period, we arrive at the contacts of the key group being 6.
# TODO (hwr): I still think I could understand this a little better..

dec_effective_R0 = 7    # effective R0 (december)
generation_time = 4     # generation time (days)
dec_contacts_of_key_group = 6
DEC_INFECTIONS_PER_CONTACT = dec_effective_R0 / dec_contacts_of_key_group
DEC_DAYS_INFECTIOUS = micro.days_infectious(7,2)

# Based on Xiangyu's adjusted moment match code, this is the fraction of the
# population (UG-only?  or all students) broken out by the amount of contact
# that they have, starting from low contact to high contact.
POP_FRAC = np.array(
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
MARGINAL_CONTACTS = np.arange(1,K+1)


def main():
    # ===================================================================
    # [Initialize] Assume no recovered and set initial Omicron infections
    # ===================================================================

    pop = TOTAL_POP * POP_FRAC
    R0 = np.zeros(K)
    I0 = np.zeros(K)
    I0[INITIAL_GROUP] = INITIAL
    S0 = np.maximum(pop - R0 - I0, 0)

    # ========================================================================
    # [Run] Increase testing delay and reduce interaction over duration of sim
    # ========================================================================

    # [12/1 to 12/9] 1x / week testing with 36hr delay
    infections_per_contact = DEC_INFECTIONS_PER_CONTACT * micro.days_infectious(7,1.5) / DEC_DAYS_INFECTIOUS
    infection_rate = meta_group("UG", pop, MARGINAL_CONTACTS).infection_matrix(infections_per_contact)
    s = sim(T, S0, I0, R0, infection_rate=infection_rate, generation_time=generation_time)
    s.step(2)

    # [12/10 to 12/16] 1x / week testing with 3 day delay
    infections_per_contact = DEC_INFECTIONS_PER_CONTACT * micro.days_infectious(7,3) / DEC_DAYS_INFECTIOUS
    infection_rate = meta_group("UG", pop, 0.5 * MARGINAL_CONTACTS).infection_matrix(infections_per_contact)
    s.step(1, infection_rate=infection_rate)

    infections_per_contact = DEC_INFECTIONS_PER_CONTACT * micro.days_infectious(7,3) / DEC_DAYS_INFECTIOUS
    infection_rate = meta_group("UG", pop, 0.33 * MARGINAL_CONTACTS).infection_matrix(infections_per_contact)
    s.step(1, infection_rate=infection_rate)

    # ==================================================
    # [Plot] Actual Infections vs. Simulation Infections
    # ==================================================

    I = s.get_metric('I', aggregate=True, cumulative=True)
    plt.plot(np.arange(T)*generation_time, I, label="simulated")

    # plot actual counts
    dec_daily_positives = list(pd.read_csv("data/dec_infections.csv")['positives'])
    dec_positives = np.cumsum(dec_daily_positives)
    plt.plot(np.arange(15), dec_positives[:15], label="actual")

    plt.title(f"Actual vs. Simulated Infection Trajectories\nDec Effective R0 = {dec_effective_R0}, Symptomatic Rate = {SYMPTOMATIC_RATE}")
    plt.rcParams.update({'font.size': 8})
    plt.xlabel('Days Since Dec. 1')
    plt.ylabel('Cumulative Infected')
    plt.legend()
    plt.savefig('calibration.png', facecolor='w')
    plt.close()


if __name__ == "__main__":
    main()
