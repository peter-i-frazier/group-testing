import sys
import numpy as np
import json
import yaml
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

    # If set, this replaces the detailed description of parameters in the plot with a simple summary
    if 'simple_param_summary' in params:
        SIMPLE_PARAM_SUMMARY = params['simple_param_summary']
    else:
        SIMPLE_PARAM_SUMMARY = None

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

    def sim_test_regime(tests_per_week, delay, color, label = None):

        regime = TestingRegime(popul,tests_per_week,delay,params)

        if label is None:
            label = regime.get_name()

        infections_per_contact_unit = BOOSTER_EFFECTIVENESS * \
                                      np.multiply(INFECTIONS_PER_DAY_PER_CONTACT_UNIT, \
                                                  regime.get_days_infectious())

        # This infection rate is before we have applied testing, boosters,
        # and any period-specific changes to transmission rates. It is a matrix where entry [i,j]
        # is the number of
        infection_rate = popul.infection_matrix(infections_per_contact_unit)

        s = sim(T, S0, I0, R0, infection_rate=infection_rate,
                infection_discovery_frac=popul.metagroup2group(regime.get_infection_discovery_frac()),
                recovered_discovery_frac=popul.metagroup2group(regime.get_recovered_discovery_frac()),
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

    def sim_test_complex_regime(tests_per_week, delay, transmission_multipliers, \
                                period_lengths, color, label):
        '''
        tests_per_week, delay, transmission_multipliers, and period_lengths should all be lists of the same length.
        These correspond to different periods of time that we wish to simulate.
        Period i is described by the i-th element in each of these lists.

        tests_per_week[i] should be either a scalar giving a value to apply across the whole population,
        or should be a dictionary whose keys are the meta-group names (as specified in the population)
        and whose value is the number of tests to do for that meta-group in the period.

        delay[i] should be structured similarly and gives the test delay associated with this meta-group
        in this period.

        transmission_multipliers[i] is either a float or a list of floats (one for each meta-group).
        It multplies transmission by for each meta-group in this period.

        period_lengths[i] is the number of generations that this period lasts

        color is the name of a color to use in plotting for this testing regime

        label is a string to use in legends when referring to this testing regime
        '''
        for i in range(len(period_lengths)):
            regime = TestingRegime(popul, tests_per_week[i], delay[i], params)
            # infections_per_contact_unit =
            #   BOOSTER_EFFECTIVENESS * transmission_multipliers[i] * INFECTIONS_PER_DAY_PER_CONTACT_UNIT * \
            #   regime.get_days_infectious()
            infections_per_contact_unit = BOOSTER_EFFECTIVENESS * \
                                          np.multiply(transmission_multipliers[i], \
                                                      np.multiply(
                                                          INFECTIONS_PER_DAY_PER_CONTACT_UNIT, \
                                                          regime.get_days_infectious()))
            infection_rate = popul.infection_matrix(infections_per_contact_unit)
            infection_discovery_frac = popul.metagroup2group(regime.get_infection_discovery_frac())
            recovered_discovery_frac = popul.metagroup2group(regime.get_recovered_discovery_frac())
            if i == 0: # instantiate simulation object
                s = sim(T, S0, I0, R0, infection_rate,
                        infection_discovery_frac=infection_discovery_frac ,
                        recovered_discovery_frac=recovered_discovery_frac,
                        generation_time=GENERATION_TIME, outside_rate=outside_rate)
            s.step(period_lengths[i], infection_rate=infection_rate,
               infection_discovery_frac = infection_discovery_frac,
               recovered_discovery_frac = recovered_discovery_frac)

        test_regime_names.append(label)
        test_regime_sims.append(s)
        test_regime_colors.append(color)

    # This tests UG and Professional students 2x week during virtual instruction (first 6 generations),
    # and then stops surveilling them. It does not surveil GR or FS.
    # TODO pf98 rather than hard-coding the test regimes we want to plot, let's put them into a separate yaml file
    plot = 1
    if plot == 1:
        sim_test_complex_regime(
            [ { 'UG':2, 'GR':0, 'PR':2, 'FS':0},
              {'UG': 2, 'GR': 0, 'PR': 2, 'FS': 0},
              0 ], # testing frequencies.  UG and PR are tested 2x / wk in period 1
            [ 1, 1, 1], #test delay
            [ 1, CLASSWORK_TRANSMISSION_MULTIPLIER, CLASSWORK_TRANSMISSION_MULTIPLIER], # transmission multipliers
            [ 3, 3, T-6-1 ], # period lengths
            'purple', 'UG+Prof. 2x/wk in virtual instr., no other surveillance ')

    """
    TODO pf98

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

    if plot == 2:
        sim_test_regime(1,2,"crimson")
        sim_test_regime(1,1,"orangered") # used to be coral in the plots with 1.5 day delay
        sim_test_regime(2,2,"navy")
        sim_test_regime(2,1,"royalblue") # used to be powderblue in the plots with 1.5 day delay

    # plots 1 and 2
    sim_test_regime(0,1,"black") # No surveillance

    # =================
    # [Plot] Make plots
    # ==================
    if simple_plot:
        plotting.plot_sm_test_regime_comparison(out_file, test_regime_names,
            test_regime_sims, test_regime_colors, params)
    else:
        plotting.plot_comprehensive_summary(out_file, test_regime_names,
            test_regime_sims, test_regime_colors, params, popul, SIMPLE_PARAM_SUMMARY)

    plotting.plot_hospitalization('sp22_sim_hosp.png', test_regime_names, test_regime_sims, test_regime_colors, params, popul)
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
