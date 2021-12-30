import numpy as np
import pandas as pd
import micro
from sim import sim
from groups import meta_group
import matplotlib.pyplot as plt

TOTAL_POP = 10000       # total fac/staff population
K = 11                  # number of distinct contact groups
T = 5                   # num generations
INITIAL_GROUP = 10      # index of group containing initial infections
INITIAL = 1             # number of initial infections
SYMPTOMATIC_RATE = .3   # fraction of new infections will be discovered
                        # without surevillance

# We aim to calibrate the effective R0 and generation_time of the simulation
# to the observed December 2021 outbreak. We assume that the most social of
# the population of gets infected first.

# This is the parameter we aim to calibrate.  An employee in a group with a given number of contacts units
# generates dec_fs_infected_per_day_unit * contact_units infections per day generated during study week in the December
# Omicron surge
dec_fs_infected_per_day_unit = 0.16

generation_time = 4     # generation time (days)

# Based on Xiangyu's adjusted moment match code, this is the fraction of the
# population (UG-only?  or all students) broken out by the amount of contact
# that they have, starting from low contact to high contact.
POP_FRAC = np.array(
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
    0.00176556195])
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
    infections_per_contact = dec_fs_infected_per_day_unit * micro.days_infectious(7,1.5)
    infection_rate = meta_group("FS", pop, MARGINAL_CONTACTS).infection_matrix(infections_per_contact)
    s = sim(T, S0, I0, R0, infection_rate=infection_rate, generation_time=generation_time)
    s.step(2)

    # [12/10 to 12/16] 1x / week testing with 3 day delay
    infections_per_contact = dec_fs_infected_per_day_unit * micro.days_infectious(7,3)
    infection_rate = meta_group("FS", pop, 0.5 * MARGINAL_CONTACTS).infection_matrix(infections_per_contact)
    s.step(1, infection_rate=infection_rate)

    infections_per_contact = dec_fs_infected_per_day_unit * micro.days_infectious(7,3)
    infection_rate = meta_group("FS", pop, 0.33 * MARGINAL_CONTACTS).infection_matrix(infections_per_contact)
    s.step(1, infection_rate=infection_rate)

    # ==================================================
    # [Plot] Actual Infections vs. Simulation Infections
    # ==================================================

    I = s.get_metric('I', aggregate=True, cumulative=True)
    plt.plot(np.arange(T)*generation_time, I, label="simulated")

    # plot actual counts
    dec_daily_positives = list(pd.read_csv("data/dec_infections_fs.csv")['positives'])
    dec_positives = np.cumsum(dec_daily_positives)
    plt.plot(np.arange(15), dec_positives[:15], label="actual")


    plt.title(f"Actual vs. Simulated Infection Trajectories\n"
              f"dec_fs_infected_per_day_unit={dec_fs_infected_per_day_unit}, Symptomatic Rate = {SYMPTOMATIC_RATE}")
    plt.rcParams.update({'font.size': 8})
    plt.xlabel('Days Since Dec. 1')
    plt.ylabel('Cumulative Infected')
    plt.legend()
    plt.savefig('calibration_fs.png', facecolor='w')
    plt.close()


if __name__ == "__main__":
    main()
