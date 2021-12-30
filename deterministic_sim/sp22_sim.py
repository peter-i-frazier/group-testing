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
import plotting
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

NOMINAL = yaml.safe_load(open("nominal.yaml", "r"))

def main(**kwargs):
    if "yaml" in kwargs:
        params = yaml.safe_load(open(kwargs["yaml"], "r"))
    else:
        params = yaml.safe_load(open("nominal.yaml", "r"))

    params.update(kwargs)

    # Include the parameters defined in the JSON file
    json_params = json.load(open(params["json_path"]))
    params.update(json_params)

    T = params['T']
    GENERATION_TIME = params['generation_time']
    SYMPTOMATIC_RATE = params['symptomatic_rate']
    NO_SURVEILLANCE_TEST_RATE = params["no_surveillance_test_rate"]
    R0_REDUCTION = params['R0_reduction']
    BOOSTER_EFFECTIVENESS = params['booster_effectiveness']
    INFECTIONS_PER_DAY_PER_CONTACT_UNIT = \
        np.array(params['infections_per_day_per_contact_unit'])

    # =====================================================================
    # [Initialize] Assume a group's previous and new infections are divided
    # proportionally to the amount of contact it has as a group.
    # =====================================================================

    population_count = params["population_count"]
    population_names = params["population_names"]
    initial_infections = params['initial_infections']
    past_infections = params['past_infections']
    meta_groups = []
    for i in range(len(population_count)):
        name = population_names[i]
        pop = population_count[i] * np.array(params['pop_fracs'][i])
        contact_units = np.arange(1, len(pop) + 1)
        meta_groups.append(meta_group(name, pop, contact_units))

    popul = population(meta_groups, np.array(params['meta_matrix']))
    S0, I0, R0 = popul.get_init_SIR_vec(initial_infections, past_infections)

    # ========================================
    # [Run] Reduce R0 once the semester begins
    # ========================================

    test_regime_names = []
    test_regime_sims = []
    test_regime_colors = []
    def sim_test_regime(tests_per_week, delay, color):
        days_between_tests = 7 / tests_per_week
        infections_per_contact_unit = BOOSTER_EFFECTIVENESS * INFECTIONS_PER_DAY_PER_CONTACT_UNIT * micro.days_infectious(days_between_tests, delay)
        infection_rate = popul.infection_matrix(infections_per_contact_unit)
        s = sim(T, S0, I0, R0, infection_rate=infection_rate,
                generation_time=GENERATION_TIME)
        s.step(4)
        s.step(T-1-4, infection_rate=R0_REDUCTION * infection_rate)

        label = "%dx/wk, %.1fd delay" % (tests_per_week, delay)
        test_regime_names.append(label)
        test_regime_sims.append(s)
        test_regime_colors.append(color)

    sim_test_regime(1,2,"crimson")
    sim_test_regime(1,1.5,"orangered")
    sim_test_regime(1,1,"coral")
    sim_test_regime(2,2,"navy")
    sim_test_regime(2,1.5,"royalblue")
    sim_test_regime(2,1,"powderblue")

    # No surveillance
    infections_per_contact_unit = BOOSTER_EFFECTIVENESS * INFECTIONS_PER_DAY_PER_CONTACT_UNIT * micro.days_infectious(np.inf,1)
    infection_discovery_frac = SYMPTOMATIC_RATE
    recovered_discovery_frac = NO_SURVEILLANCE_TEST_RATE
    infection_rate = popul.infection_matrix(infections_per_contact_unit)
    s = sim(T, S0, I0, R0, infection_rate=infection_rate,
            infection_discovery_frac=infection_discovery_frac,
            recovered_discovery_frac=recovered_discovery_frac,
            generation_time=GENERATION_TIME)
    s.step(4)
    s.step(T-1-4, infection_rate=R0_REDUCTION * infection_rate)

    # ====================================
    # [Plot] Comparison of testing regimes
    # ====================================

    def old_plot():
        plotting.plot_sm_test_regime_comparison(test_regime_names,
            test_regime_sims, test_regime_colors, s, params)

    def new_plot():
        plotting.plot_comprehensive_summary(s, popul, params)

    old_plot()


if __name__ == "__main__":
    _override_params = dict(arg.split('=') for arg in sys.argv[1:])
    override_params = {}
    for k, v in _override_params.items():
        if k == "yaml":
            override_params[k] = v
        elif k == "T":
            override_params[k] = int(v)
        else:
            override_params[k] = float(v)
    main(**override_params)
