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
    DEC_INFECTIONS_PER_CONTACT = params['dec_effective_R0'] / \
                                 params['dec_contacts_of_key_group']
    DEC_DAYS_INFECTIOUS = micro.days_infectious(params['dec_days_between_tests'],
                                                params['dec_isolation_delay'])

    # =====================================================================
    # [Initialize] Assume a group's previous and new infections are divided
    # proportionally to the amount of contact it has as a group.
    # =====================================================================

    population_count = params["population_count"]
    population_names = params["population_names"]
    meta_groups = []
    for i in range(len(population_count)):
        name = population_names[i]
        pop = population_count[i] * np.array(params['pop_fracs'][i])
        contact_units = np.arange(1, len(pop) + 1)
        meta_groups.append(meta_group(name, pop, contact_units))

    initial_infectious = params['initial_infections']
    initial_recovered = params['infected_from_outbreak'] + params['infected_over_break']

    popul = population(meta_groups, np.array(params['meta_matrix']))
    S0, I0, R0 = popul.get_init_SIR(initial_infectious, initial_recovered)

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
        isolated = s.get_isolated(iso_lengths=params["isolation_durations"],
                                  iso_props=params["isolation_fracs"],
                                  on_campus_frac=params["on_campus_frac"])
        plt.plot(np.arange(T)*GENERATION_TIME, isolated, label=label, color=color)

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
    isolated = s.get_isolated(iso_lengths=params["isolation_durations"],
                              iso_props=params["isolation_fracs"],
                              on_campus_frac=params["on_campus_frac"])
    plt.plot(np.arange(T) * GENERATION_TIME, isolated, 'k', label='No surveillance')


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
