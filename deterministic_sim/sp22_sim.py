import sys
import numpy as np
import yaml
import micro
from sim import sim
from groups import well_mixed_infection_rate_one_meta_group
import matplotlib
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

NOMINAL = yaml.safe_load(open("nominal.yaml", "r"))

# TODO (hwr26): auto genrerate this
# These parameters are passed to get_isolated to get isolation days for different populations of students
# under different isolation policies
ISOLATION_LEN = 3
ONCAMPUS_FRAC = 0.5 # fraction of students that live on campus
ISOLATION_FRAC_ALL_5DAY = np.array([1, .4, .1]) # 80% 5 day / 20% 10 day
ISOLATION_FRAC_ALL_10DAY = np.array([1, 1, .5]) # 100% 10 day
ISOLATION_FRAC_ONCAMPUS_5DAY = ONCAMPUS_FRAC * ISOLATION_FRAC_ALL_5DAY
ISOLATION_FRAC_ONCAMPUS_10DAY = ONCAMPUS_FRAC * ISOLATION_FRAC_ALL_10DAY

def main(**kwargs):

    params = NOMINAL
    params.update(kwargs)

    T = params['T']
    GENERATION_TIME = params['generation_time']
    SYMPTOMATIC_RATE = params['symptomatic_rate']
    BOOSTER_EFFECTIVENESS = params['booster_effectiveness']
    R0_REDUCTION = params['R0_reduction']
    DEC_INFECTIONS_PER_CONTACT = params['dec_effective_R0'] / \
                                 params['dec_contacts_of_key_group']
    DEC_DAYS_INFECTIOUS = micro.days_infectious(params['dec_days_between_tests'],
                                                params['dec_isolation_delay'])

    # =====================================================================
    # [Initialize] Assume a group's previous and new infections are divided
    # proportionally to the amount of contact it has as a group.
    # =====================================================================

    pop = params['total_pop'] * np.array(params['pop_frac'])
    MARGINAL_CONTACTS = np.arange(1, params['K']+1)
    b = MARGINAL_CONTACTS * np.array(params['pop_frac'])
    b =  b / np.sum(b)
    total_initial = params['infected_from_outbreak'] + params['infected_over_break']
    R0 = total_initial * b
    I0 = params['initial_infections'] * b
    S0 = np.maximum(pop - R0 - I0, 0)

    # ========================================
    # [Run] Reduce R0 once the semester begins
    # ========================================

    def sim_test_regime(tests_per_week, delay, color):
        """Simulate a testing regime"""
        days_between_tests = 7 / tests_per_week
        infections_per_contact = BOOSTER_EFFECTIVENESS * DEC_INFECTIONS_PER_CONTACT * micro.days_infectious(days_between_tests,delay) / DEC_DAYS_INFECTIOUS
        infection_rate = well_mixed_infection_rate_one_meta_group(pop, MARGINAL_CONTACTS, infections_per_contact)
        s = sim(T, S0, I0, R0, infection_rate=infection_rate, generation_time=GENERATION_TIME)
        s.step(4)
        infection_rate = well_mixed_infection_rate_one_meta_group(pop, R0_REDUCTION * MARGINAL_CONTACTS, infections_per_contact)
        s.step(T-1-4, infection_rate=infection_rate)

        label = "%dx/wk, %.1fd delay" % (tests_per_week, delay)
        plt.subplot(211)
        plt.plot(np.arange(T)*GENERATION_TIME/7, s.get_discovered(aggregate=True,cumulative=True), label=label, color=color)
        plt.subplot(212)
        isolated = s.get_isolated(isolation_len=ISOLATION_LEN, isolation_frac=ISOLATION_FRAC_ONCAMPUS_5DAY)
        plt.plot(np.arange(T)*GENERATION_TIME/7, isolated, label=label, color=color)

    sim_test_regime(1,2,"crimson")
    sim_test_regime(1,1.5,"orangered")
    sim_test_regime(1,1,"coral")
    sim_test_regime(2,2,"navy")
    sim_test_regime(2,1.5,"royalblue")
    sim_test_regime(2,1,"powderblue")

    # No surveillance
    infections_per_contact = BOOSTER_EFFECTIVENESS * DEC_INFECTIONS_PER_CONTACT * micro.days_infectious(np.inf,1) / DEC_DAYS_INFECTIOUS
    infection_rate = well_mixed_infection_rate_one_meta_group(pop, MARGINAL_CONTACTS, infections_per_contact)
    infection_discovery_frac = SYMPTOMATIC_RATE
    recovered_discovery_frac = .01 # 1% of the population is tested for any reason in a given generation
    s = sim(T, S0, I0, R0, infection_rate=infection_rate, generation_time=GENERATION_TIME,
        infection_discovery_frac=infection_discovery_frac,
        recovered_discovery_frac=recovered_discovery_frac)
    s.step(T-1)
    plt.subplot(211)
    plt.plot(np.arange(T)*GENERATION_TIME/7, s.get_infected(aggregate=True,cumulative=True), 'k--', label='No surveillance, Infected')
    plt.plot(np.arange(T)*GENERATION_TIME/7, s.get_discovered(aggregate=True,cumulative=True), 'k-', label='No surveillance, Discovered')
    plt.subplot(212)
    isolated = s.get_isolated(isolation_len=ISOLATION_LEN, isolation_frac=ISOLATION_FRAC_ONCAMPUS_5DAY)
    plt.plot(np.arange(T) * GENERATION_TIME/7, isolated, 'k', label='No surveillance')

    # ====================================
    # [Plot] Comparison of testing regimes
    # ====================================

    plt.subplot(211)
    plt.title(f"Dec Effective R0 = {params['dec_effective_R0']}, Symptomatic Rate = {SYMPTOMATIC_RATE}")
    plt.rcParams.update({'font.size': 8})
    plt.legend()
    plt.ylabel('UG Infected')

    plt.subplot(212)
    plt.rcParams.update({'font.size': 8})
    plt.legend()
    plt.xlabel('Weeks')
    plt.ylabel('UG in Isolation (on-campus 5day)')

    plt.savefig('sp22_sim.png', facecolor='w')
    plt.close()


if __name__ == "__main__":
    override_params = dict(arg.split('=') for arg in sys.argv[1:])
    override_params = {k:float(v) for k,v in override_params.items()}
    main(**override_params)
