import sys
import numpy as np
import json
import yaml
from plotting import Trajectory
from strategy import Strategy
import micro
from sim import sim
from groups import meta_group, population
from testing_regime import TestingRegime
import matplotlib
import matplotlib.pyplot as plt
import warnings
import plotting
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


def initialize_population(params):
    """Initialize the population from the simulation params."""

    population_count = params["population_count"]
    population_names = params["population_names"]
    meta_groups = []
    for i in range(len(population_count)):
        name = population_names[i]
        pop = population_count[i] * np.array(params['pop_fracs'][i])
        contact_units = np.arange(1, len(pop) + 1)
        meta_groups.append(meta_group(name, pop, contact_units))
    return population(meta_groups, np.array(params['meta_matrix']))


def main(yaml_file='nominal.yaml', simple_plot=False, out_file='sp22_sim.png', **kwargs):

    # =======================
    # [Initialize Parameters]
    # =======================

    params = yaml.safe_load(open(yaml_file, "r"))
    params.update(kwargs)

    # Include the parameters defined in the JSON file
    json_params = json.load(open(params["json_path"]))
    params.update(json_params)

    T = params['T']
    GENERATION_TIME = params['generation_time']
    CLASSWORK_TRANSMISSION_MULTIPLIER = params['classwork_transmission_multiplier']
    BOOSTER_EFFECTIVENESS = params['booster_effectiveness']
    INFECTIONS_PER_DAY_PER_CONTACT_UNIT = \
        np.array(params['infections_per_day_per_contact_unit'])

    # If set, this replaces the detailed description of parameters in the plot with a simple summary
    if 'simple_param_summary' in params:
        SIMPLE_PARAM_SUMMARY = params['simple_param_summary']
    else:
        SIMPLE_PARAM_SUMMARY = None

    # ========================================================================
    # [Initialize Strategies]
    # TODO (hwr26): At a future point, would be nice to store this information
    # in multiple YAML files.
    # ========================================================================

    popul = initialize_population(params)

    # some common testing regimes used in these strategies
    no_testing_testing_regime=TestingRegime(popul=popul, tests_per_week=0,
                                            test_delay=1, params=params)
    ug_prof_2x_week_testing_regime= \
        TestingRegime(popul=popul, tests_per_week={ 'UG':2, 'GR':0, 'PR':2, 'FS':0},
                      test_delay=1, params=params)

    # No extra testing (note all strategies assume testing for cause)
    no_testing_strategy = \
        Strategy(name="No Testing",
             pct_discovered_in_pre_departure=0,
             pct_discovered_in_arrival_test=0,
             testing_regimes=[no_testing_testing_regime, no_testing_testing_regime],
             transmission_multipliers=[1, CLASSWORK_TRANSMISSION_MULTIPLIER],
             period_lengths=[3,T-3-1])

    # SCENARIO 1:
    # all students do a pre-departure antigen test with sensitivity 50%
    # all students again an arrival PCR test with sensitivity 67% and no isolation delay
    #     0.5 caught in pre-departure
    #     (1 - 0.5) * 0.67 = 0.335 caught in arrival testing
    #     1 - 0.5 - 0.335 = 0.165 not caught

    # SCENARIO 2:
    # half of students do a pre-departure test with sensitivity 50%
    # 75% of students (independently chosen) do an arrival PCR 67% sensitivity and no test delay
    #     0.5 * 0.5 = 0.25 caught in pre-departure
    #     (1 - 0.25) * (0.75 * 0.67) = 0.38 caught in arrival testing
    #     1 - 0.25 - 0.38 = 0.37 not caught

    # Pre-departure + arrival testing. No surveillance at any point
    arrival_testing_strategy = \
        Strategy(name="Only Pre-Departure + Arrival Testing",
             pct_discovered_in_pre_departure=0.25,   # Using Scenario 2
             pct_discovered_in_arrival_test=0.38,  # Using Scenario 2
             testing_regimes=[no_testing_testing_regime, no_testing_testing_regime],
             transmission_multipliers=[1, CLASSWORK_TRANSMISSION_MULTIPLIER],
             period_lengths=[3,T-3-1])

    # Pre-departure + arrival testing. Surveillance of UG and professional
    # students before classes and during virtual instruction at 2x/wk.
    # It does not surveil GR or FS.
    surge_testing_strategy = \
        Strategy(name="UG+Prof. 2x/wk in Virtual Instr. Only",
             pct_discovered_in_pre_departure=0.25,   # Using Scenario 2
             pct_discovered_in_arrival_test=0.38,  # Using Scenario 2
             testing_regimes=[ug_prof_2x_week_testing_regime,
                              ug_prof_2x_week_testing_regime,
                              no_testing_testing_regime],
             transmission_multipliers=[1, CLASSWORK_TRANSMISSION_MULTIPLIER,
                                       CLASSWORK_TRANSMISSION_MULTIPLIER],
             period_lengths=[3,3,T-6-1])

    # ==================================
    # [Run] Compare a list of strategies
    # ==================================


    def sim_test_regime(tests_per_week, delay, color):

        regime = TestingRegime(popul,tests_per_week,delay,params)

        strategy = \
            Strategy(name=regime.name,
                # TODO (hwr26): This matches what we had before but we should
                # choose a scenario to run with
                pct_discovered_in_pre_departure=0.25,   # Using Scenario 2
                pct_discovered_in_arrival_test=0.38,  # Using Scenario 2
                testing_regimes=[regime, regime],
                transmission_multipliers=[1, CLASSWORK_TRANSMISSION_MULTIPLIER],
                period_lengths=[3,T-3-1])

        sim_test_strategy(strategy, color)


    def sim_test_strategy(strategy:Strategy, color:str):
        for i in range(strategy.periods):
            regime = strategy.testing_regimes[i]
            # infections_per_contact_unit =
            #   BOOSTER_EFFECTIVENESS * transmission_multipliers[i] * INFECTIONS_PER_DAY_PER_CONTACT_UNIT * \
            #   regime.get_days_infectious()
            infections_per_contact_unit = BOOSTER_EFFECTIVENESS * \
                                          np.multiply(strategy.transmission_multipliers[i], \
                                                      np.multiply(
                                                          INFECTIONS_PER_DAY_PER_CONTACT_UNIT, \
                                                          regime.get_days_infectious()))
            infection_rate = popul.infection_matrix(infections_per_contact_unit)
            infection_discovery_frac = popul.metagroup2group(regime.get_infection_discovery_frac())
            recovered_discovery_frac = popul.metagroup2group(regime.get_recovered_discovery_frac())

            initial_infections = strategy.get_initial_infections(params)
            past_infections = strategy.get_past_infections(params)
            S0, I0, R0 = popul.get_init_SIR_vec(initial_infections, past_infections,
                                                weight="population x contacts")
            outside_rates = params['outside_rates']
            outside_rate = popul.get_outside_rate(outside_rates)

            if i == 0: # instantiate simulation object
                s = sim(T, S0, I0, R0, infection_rate,
                        infection_discovery_frac=infection_discovery_frac ,
                        recovered_discovery_frac=recovered_discovery_frac,
                        generation_time=GENERATION_TIME, outside_rate=outside_rate)
            s.step(strategy.period_lengths[i], infection_rate=infection_rate,
               infection_discovery_frac = infection_discovery_frac,
               recovered_discovery_frac = recovered_discovery_frac)

        return Trajectory(strategy, s, color)

    """
    TODO (hwr26): I don't think this is necessary anymore..

    I ran these three different simulations and confirmed via inspection of the plots that they
    gave the same output is the same. It would be good to move this into a test case.

    sim_test_regime(2,1,"powderblue")

    sim_test_complex_regime(
        [ { 'UG':2, 'GR':2, 'PR':2, 'FS':2}, 2 ], # testing frequencies.  UG and PR are tested 2x / wk in period 1
        [ 1, 1], #test delay
        [ 1, CLASSWORK_TRANSMISSION_MULTIPLIER], # transmission multipliers
        [ 3, T-3-1 ], # period lengths
        'powderblue', '2x/wk test')

    sim_test_complex_regime(
        [ 2, 2 ], # testing frequencies.  UG and PR are tested 2x / wk in period 1
        [ 1, 1], #test delay
        [ 1, CLASSWORK_TRANSMISSION_MULTIPLIER], # transmission multipliers
        [ 3, T-3-1 ], # period lengths
        'powderblue', '2x/wk test')
    """
    plot = 3

    if plot == 1:
        trajectories = [
            sim_test_strategy(surge_testing_strategy, 'purple'),
            sim_test_regime(0,1,"black") # No surveillance
        ]
    if plot == 2:
        trajectories = [
            sim_test_regime(1,2,"crimson"),
            sim_test_regime(1,1,"orangered"), # used to be coral in the plots with 1.5 day delay
            sim_test_regime(2,2,"navy"),
            sim_test_regime(2,1,"royalblue"), # used to be powderblue in the plots with 1.5 day delay
            sim_test_regime(0,1,"black") # No surveillance
        ]
    if plot == 3:
        trajectories = [
            sim_test_strategy(no_testing_strategy, 'black'),
            sim_test_strategy(arrival_testing_strategy, 'royalblue'),
            sim_test_strategy(surge_testing_strategy, 'purple')
        ]

    # =================
    # [Plot] Make plots
    # ==================

    if simple_plot:
        plotting.plot_sm_test_regime_comparison(out_file, trajectories, params)
    else:
        plotting.plot_comprehensive_summary(out_file, trajectories, params, popul, SIMPLE_PARAM_SUMMARY)

    plotting.plot_hospitalization('sp22_sim_hosp.png', trajectories, params, popul)


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
        elif k == 'simple_param_summary':
            override_params[k] = v # string
        else:
            override_params[k] = float(v)

    main(yaml_file, simple_plot, out_file, **override_params)
