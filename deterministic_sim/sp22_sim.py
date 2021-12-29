import sys
import numpy as np
import json
import yaml
import micro
from sim import sim
from groups import meta_group, population
import matplotlib
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

NOMINAL = yaml.safe_load(open("nominal.yaml", "r"))

def main(**kwargs):

    params = NOMINAL
    params.update(kwargs)

    # Include the parameters defined in the JSON file
    json_params = json.load(open(params["json_path"]))
    params.update(json_params)

    T = params['T']
    GENERATION_TIME = params['generation_time']
    SYMPTOMATIC_RATE = params['symptomatic_rate']
    BOOSTER_EFFECTIVENESS = params['booster_effectiveness']
    R0_REDUCTION = params['R0_reduction']
    DEC_INFECTIONS_PER_CONTACT = params['dec_effective_R0'] / \
                                 params['dec_contacts_of_key_group']
    DEC_DAYS_INFECTIOUS = micro.days_infectious(params['dec_days_between_tests'],
                                                params['dec_isolation_delay'])

    # TODO (hwr26): Update with Sam's auto-generation of these arrays
    on_campus_frac = params['on_campus_frac']
    isolation_frac_5_day = np.array([1, .4, .1]) # 80% 5 day / 20% 10 day
    isolation_frac_10_day = np.array([1, 1, .5]) # 100% 10 day
    ISOLATION_LEN = 3
    ISOLATION_FRAC_ON_CAMPUS_5DAY = on_campus_frac * isolation_frac_5_day
    ISOLATION_FRAC_ON_CAMPUS_10DAY = on_campus_frac * isolation_frac_10_day

    # =====================================================================
    # [Initialize] Assume a group's previous and new infections are divided
    # proportionally to the amount of contact it has as a group.
    # =====================================================================

    # TODO (hwr26): To be replaced with code Sam is working on
    pops = params["total_pops"][:]
    for i in range(len(params["total_pops"])):
        pops[i] = params["total_pops"][i]*np.array(params["pop_fracs"][i])

    marginal_contacts = np.array([np.arange(1,len(pops[0])+1),
                                  np.arange(1,len(pops[1])+1),
                                  np.arange(1,len(pops[2])+1),
                                  np.arange(1,len(pops[3])+1)])


    group_names = ['UG', 'GR', 'PR', 'FS']
    meta_groups = [meta_group(group_names[i], pops[i], marginal_contacts[i]) \
                    for i in range(4)]
    popul = population(meta_groups, np.array(params['meta_matrix']))

    marginal_contacts_flat = popul.flatten(marginal_contacts)
    pops_flat = popul.flatten(pops)
    b =  marginal_contacts_flat * pops_flat
    b =  b / np.sum(b)
    total_initial = params['infected_from_outbreak'] + params['infected_over_break']
    R0 = total_initial * b
    I0 = params['initial_infections'] * b
    S0 = np.maximum(pops_flat - R0 - I0, 0)

    # ========================================
    # [Run] Reduce R0 once the semester begins
    # ========================================

    def sim_test_regime(tests_per_week, delay, color):
        days_between_tests = 7 / tests_per_week
        infections_per_contact = BOOSTER_EFFECTIVENESS * DEC_INFECTIONS_PER_CONTACT * micro.days_infectious(days_between_tests, delay) / DEC_DAYS_INFECTIOUS
        infection_matrix = popul.infection_matrix(infections_per_contact)
        s = sim(T, S0, I0, R0, infection_rate=infection_matrix,
                generation_time=GENERATION_TIME)
        s.step(T-1)

        label = "%dx/wk, %.1fd delay" % (tests_per_week, delay)
        plt.subplot(211)
        plt.plot(np.arange(T)*GENERATION_TIME, s.get_discovered(aggregate=True,cumulative=True), label=label, color=color)
        plt.subplot(212)
        plt.plot(np.arange(T)*GENERATION_TIME, s.get_isolated(), label=label, color=color)

    sim_test_regime(1,2,"crimson")
    sim_test_regime(1,1.5,"orangered")
    sim_test_regime(1,1,"coral")
    sim_test_regime(2,2,"navy")
    sim_test_regime(2,1.5,"royalblue")
    sim_test_regime(2,1,"powderblue")

    # No surveillance
    infections_per_contact = BOOSTER_EFFECTIVENESS * DEC_INFECTIONS_PER_CONTACT * micro.days_infectious(np.inf,1) / DEC_DAYS_INFECTIOUS
    infection_discovery_frac = SYMPTOMATIC_RATE
    recovered_discovery_frac = .01 # 1% of the population is tested for any reason in a given generation
    infection_rate=popul.infection_matrix(infections_per_contact)
    s = sim(T, S0, I0, R0, infection_rate=infection_rate,
            infection_discovery_frac=infection_discovery_frac,
            recovered_discovery_frac=recovered_discovery_frac,
            generation_time=GENERATION_TIME)
    s.step(T-1)

    plt.subplot(211)
    plt.plot(np.arange(T)*GENERATION_TIME, s.get_discovered(aggregate=True,cumulative=True), 'k-', label='No surveillance, Discovered')
    plt.plot(np.arange(T)*GENERATION_TIME, s.get_infected(aggregate=True,cumulative=True), 'k--', label='No surveillance, Infected')
    plt.subplot(212)
    plt.plot(np.arange(T) * GENERATION_TIME, s.get_isolated(), 'k', label='No surveillance')


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
