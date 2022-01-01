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

def main(yaml_file='nominal.yaml', simple_plot=False, out_file='sp22_sim.png', **kwargs):

    params = yaml.safe_load(open(yaml_file, "r"))
    params.update(kwargs)

    # Include the parameters defined in the JSON file
    json_params = json.load(open(params["json_path"]))
    params.update(json_params)

    T = params['T']
    GENERATION_TIME = params['generation_time']
    SYMPTOMATIC_RATE = params['symptomatic_rate']
    NO_SURVEILLANCE_TEST_RATE = params["no_surveillance_test_rate"]
    CLASSWORK_TRANSMISSION_MULTIPLIER = params['classwork_transmission_multiplier']
    BOOSTER_EFFECTIVENESS = params['booster_effectiveness']
    INFECTIONS_PER_DAY_PER_CONTACT_UNIT = \
        np.array(params['infections_per_day_per_contact_unit'])
    MAX_INFECTIOUS_DAYS = params['max_infectious_days']
    PCR_SENSITIVITY = params['pcr_sensitivity']

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
    S0, I0, R0 = popul.get_init_SIR_vec(initial_infections, past_infections,
                                        weight="population x contacts")
    outside_rates = params['outside_rates']
    outside_rate = popul.get_outside_rate(outside_rates)


    # ==================================================
    # [Run] Reduce transmission once the semester begins
    # ==================================================

    test_regime_names = []
    test_regime_sims = []
    test_regime_colors = []
    def sim_test_regime(tests_per_week, delay, color):
        if tests_per_week == 0: # No surveillance
            days_between_tests = np.inf
            infection_discovery_frac = SYMPTOMATIC_RATE
            recovered_discovery_frac = NO_SURVEILLANCE_TEST_RATE
            label = "No surveillance"
        else:
            days_between_tests = 7 / tests_per_week
            infection_discovery_frac = 1
            recovered_discovery_frac = 1
            label = "%dx/wk, %.1fd delay" % (tests_per_week, delay)

        days_infectious = micro.days_infectious(days_between_tests, delay, sensitivity = PCR_SENSITIVITY, \
                                                max_infectious_days=MAX_INFECTIOUS_DAYS)
        infections_per_contact_unit = BOOSTER_EFFECTIVENESS * INFECTIONS_PER_DAY_PER_CONTACT_UNIT * days_infectious

        infection_rate = popul.infection_matrix(infections_per_contact_unit)
        s = sim(T, S0, I0, R0, infection_rate=infection_rate,
                infection_discovery_frac=infection_discovery_frac,
                recovered_discovery_frac=recovered_discovery_frac,
                generation_time=GENERATION_TIME,
                outside_rate=outside_rate)
        # 3 generations = 12 days, modeling 7 days of no instruction, followed by roughly a week of not much HW in
        # the first week of virtual classes
        s.step(3)
        infections_per_contact_unit = CLASSWORK_TRANSMISSION_MULTIPLIER * infections_per_contact_unit
        infection_rate = popul.infection_matrix(infections_per_contact_unit)
        s.step(T-1-3, infection_rate=infection_rate)

        test_regime_names.append(label)
        test_regime_sims.append(s)
        test_regime_colors.append(color)

    sim_test_regime(1,2,"crimson")
    sim_test_regime(1,1.5,"orangered")
    sim_test_regime(1,1,"coral")
    sim_test_regime(2,2,"navy")
    sim_test_regime(2,1.5,"royalblue")
    sim_test_regime(2,1,"powderblue")
    sim_test_regime(0,1,"black") # No surveillance

    # =================
    # [Plot] Make plots
    # ==================

    if simple_plot:
        plotting.plot_sm_test_regime_comparison(out_file, test_regime_names,
            test_regime_sims, test_regime_colors, params)
    else:
        plotting.plot_comprehensive_summary(out_file, test_regime_names,
            test_regime_sims, test_regime_colors, params, popul)

def usage():
    ''' Print usage message '''
    print('Usage:')
    print('python sp22_sim.py [-h] [--help] [--yaml=file.yaml] [--simple] [--out=file.png] [yaml-overrides]')
    print('--simple makes simple plots, rather than comprehensive ones')
    print('--yaml=file.yaml uses file.yaml instead of nominal.yaml')
    print('--outfile=file.png saves the plot to file.png instead of sp22_sim.png')
    print('--help and -h print this message and exit')
    print('yaml-overrides replace parameters in the yaml file and are of the form parameter_name=value')
    print('Example: python sp22_sim.py --yaml=nominal_ug.yaml --out=nominal_ug.png T=10')

if __name__ == "__main__":

    # Default values for arguments not specified in the yaml
    yaml_file = 'nominal.yaml'
    simple_plot = False # whether to make simple or comprehensive plots
    out_file = 'sp22_sim.png' # The filename for the plots

    # Parameters from the yaml file to override
    override_params = {}

    for arg in sys.argv[1:]:
        # Handle arguments without an equal sign
        if arg == '-h' or arg == '--help':
            usage()
            exit(0)
            continue
        elif arg == '--simple':
            simple_plot = True
            continue

        # If we get to this point, the argument should have the form key=value
        # Arguments that are not part of the things specified in the YAML file
        assert(len(arg.split('='))==2)
        k,v = arg.split('=') # key, value

        # Arguments that are not YAML overrides
        if k == "--yaml":
            yaml_file = str(v)
        elif k == '--out':
            out_file = str(v)

        # YAML overrides
        elif k == "T":
            # This yaml over
            override_params[k] = int(v)
        else:
            override_params[k] = float(v)

    main(yaml_file, simple_plot, out_file, **override_params)
