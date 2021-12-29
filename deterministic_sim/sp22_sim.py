import numpy as np
import micro
from sim import sim, well_mixed_infection_rate
import matplotlib
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)


TOTAL_POP = 16000               # total UG population
K = 12                          # number of distinct contact groups
T = 20                          # num generations
SYMPTOMATIC_RATE = .3           # fraction of new infections will be discovered
                                # without surevillance
BOOSTER_EFFECTIVENESS = 0.5     # reduction in transmission from boosters

infected_from_outbreak = 1900
infected_over_break = 175       # According to @Brian Liu's analysis
INFECTED_BEFORE_SEMESTER = infected_from_outbreak + infected_over_break
INITIAL_INFECTIONS = 25         # According to @Brian Liu's analysis

dec_effective_R0 = 7            # calibrated from calibration.py by @hwr26
R0_REDUCTION = 0.85             # reduction in R0 once classes start
generation_time = 4             # generation time (days)
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

    # =====================================================================
    # [Initialize] Assume a group's previous and new infections are divided
    # proportionally to the amount of contact it has as a group.
    # =====================================================================

    pop = TOTAL_POP * POP_FRAC
    b = MARGINAL_CONTACTS * POP_FRAC
    b =  b / np.sum(b)
    R0 = INFECTED_BEFORE_SEMESTER * b
    I0 = INITIAL_INFECTIONS * b
    S0 = np.maximum(pop - R0 - I0, 0)

    # ========================================
    # [Run] Reduce R0 once the semester begins
    # ========================================

    def sim_test_regime(tests_per_week, delay, color):
        """Simulate a testing regime"""
        days_between_tests = 7 / tests_per_week
        infections_per_contact = BOOSTER_EFFECTIVENESS * DEC_INFECTIONS_PER_CONTACT * micro.days_infectious(days_between_tests,delay) / DEC_DAYS_INFECTIOUS
        infection_rate = well_mixed_infection_rate(pop, MARGINAL_CONTACTS, infections_per_contact)
        s = sim(T, S0, I0, R0, infection_rate=infection_rate, generation_time=generation_time)
        s.step(4)
        infection_rate = well_mixed_infection_rate(pop, R0_REDUCTION * MARGINAL_CONTACTS, infections_per_contact)
        s.step(T-1-4, infection_rate=infection_rate)

        label = "%dx/wk, %.1fd delay" % (tests_per_week, delay)
        plt.subplot(211)
        plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), label=label, color=color)
        plt.subplot(212)
        plt.plot(np.arange(T)*generation_time, s.get_isolated(), label=label, color=color)

    sim_test_regime(1,2,"crimson")
    sim_test_regime(1,1.5,"orangered")
    sim_test_regime(1,1,"coral")
    sim_test_regime(2,2,"navy")
    sim_test_regime(2,1.5,"royalblue")
    sim_test_regime(2,1,"powderblue")

    # No surveillance
    infections_per_contact = BOOSTER_EFFECTIVENESS * DEC_INFECTIONS_PER_CONTACT * micro.days_infectious(np.inf,1) / DEC_DAYS_INFECTIOUS
    infection_rate = well_mixed_infection_rate(pop, MARGINAL_CONTACTS, infections_per_contact)
    infection_discovery_frac = SYMPTOMATIC_RATE
    recovered_discovery_frac = .01 # 1% of the population is tested for any reason in a given generation
    s = sim(T, S0, I0, R0, infection_rate=infection_rate, generation_time=generation_time,
        infection_discovery_frac=infection_discovery_frac,
        recovered_discovery_frac=recovered_discovery_frac)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*generation_time, s.get_infected(aggregate=True,cumulative=True), 'k--', label='No surveillance, Infected')
    plt.plot(np.arange(T)*generation_time, s.get_discovered(aggregate=True,cumulative=True), 'k-', label='No surveillance, Discovered')
    plt.subplot(212)
    plt.plot(np.arange(T) * generation_time, s.get_isolated(), 'k', label='No surveillance')

    # ====================================
    # [Plot] Comparison of testing regimes
    # ====================================

    plt.subplot(211)
    plt.title(f"Dec Effective R0 = {dec_effective_R0}, Symptomatic Rate = {SYMPTOMATIC_RATE}")
    plt.rcParams.update({'font.size': 8})
    plt.legend()
    plt.ylabel('UG Infected')

    plt.subplot(212)
    plt.rcParams.update({'font.size': 8})
    plt.legend()
    plt.xlabel('Weeks')
    plt.ylabel('UG in Isolation')

    plt.savefig('sp22_sim.png', facecolor='w')
    plt.close()


if __name__ == "__main__":
    main()
