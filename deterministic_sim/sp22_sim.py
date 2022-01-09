import sys
import numpy as np
import json
import yaml
from strategy import Strategy
from sim_helper import sim_test_regime, sim_test_strategy
from groups import population
from testing_regime import TestingRegime
import matplotlib
import matplotlib.pyplot as plt
import warnings
import plotting
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


def main(yaml_file='nominal.yaml', simple_plot=False, out_file='sp22_sim.png', **kwargs):

    # =======================
    # [Initialize Parameters]
    # =======================

    params = yaml.safe_load(open(yaml_file, "r"))
    params.update(kwargs)

    # Include the parameters defined in the JSON file
    json_params = json.load(open(params["json_path"]))
    params.update(json_params)

    # convert to numpy matrix
    params["meta_matrix"] = \
        np.array([list(row.values()) for row in params["meta_matrix"].values()])

    T = params['T']
    CLASSWORK_TRANSMISSION_MULTIPLIER = \
        list(params['classwork_transmission_multiplier'].values())

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

    popul = population.from_scenario(params)

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
             # pct_discovered_in_pre_departure=0.5,   # Using Scenario 1
             # pct_discovered_in_arrival_test=0.335,  # Using Scenario 1
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
             # pct_discovered_in_pre_departure=0.5,   # Using Scenario 1
             # pct_discovered_in_arrival_test=0.335,  # Using Scenario 1
             pct_discovered_in_pre_departure=0.25,  # Using Scenario 2
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

    plot = 1

    if plot == 1:
        trajectories = [
            sim_test_strategy(params, popul, surge_testing_strategy, 'purple'),
            sim_test_strategy(params, popul, arrival_testing_strategy, 'black'),
        ]
    if plot == 2:
        trajectories = [
            sim_test_regime(params,popul,1,2,"crimson"),
            sim_test_regime(params,popul,1,1,"orangered"), # used to be coral in the plots with 1.5 day delay
            sim_test_regime(params,popul,2,2,"navy"),
            sim_test_regime(params,popul,2,1,"royalblue"), # used to be powderblue in the plots with 1.5 day delay
            sim_test_regime(params,popul,0,1,"black") # No surveillance
        ]
    if plot == 3:
        trajectories = [
            sim_test_strategy(params, popul, no_testing_strategy, 'royalblue'),
            sim_test_strategy(params, popul, arrival_testing_strategy, 'black'),
            sim_test_strategy(params, popul, surge_testing_strategy, 'purple')
        ]

    # =================
    # [Plot] Make plots
    # ==================

    if simple_plot:
        plotting.plot_sm_test_regime_comparison(out_file, trajectories, params)
    else:
        plotting.plot_comprehensive_summary(out_file, trajectories, params, popul, SIMPLE_PARAM_SUMMARY)

    plotting.plot_hospitalization('sp22_sim_hosp.png', trajectories, params, popul)
    plotting.summary_statistics('sp22_sim_summary_stats.csv', trajectories, params, popul)

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
