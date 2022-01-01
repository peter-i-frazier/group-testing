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

def main(yaml_file='nominal.yaml', out_file='days_infectious.png', **kwargs):

    params = yaml.safe_load(open(yaml_file, "r"))
    params.update(kwargs)

    # Include the parameters defined in the JSON file
    json_params = json.load(open(params["json_path"]))
    params.update(json_params)

    T = params['T']
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

    GENERATION_TIME = params['generation_time']


    # ==================================================
    # [Run]
    # ==================================================

    days_between_tests = np.arange(1,20)
    delay = 1
    secondary_infections = {}
    for mgname in popul.metagroup_names():
        secondary_infections[mgname] = np.zeros(len(days_between_tests))

    for i in range(len(days_between_tests)):
        d = days_between_tests[i]
        days_infectious = micro.days_infectious(d, delay, sensitivity = PCR_SENSITIVITY, \
                                                max_infectious_days=MAX_INFECTIOUS_DAYS)
        infections_per_contact_unit = BOOSTER_EFFECTIVENESS * INFECTIONS_PER_DAY_PER_CONTACT_UNIT * days_infectious
        infection_rate = popul.infection_matrix(infections_per_contact_unit)
        # infection_rate is a matrix where infection_rate[i,j] is the number of people in group j exposed by a positive in group i

        # Instantiate a simulation object so we can get the fraction susceptible.
        # These options don't matter for this quantity:
        # infection_discovery_frac, recovered_discovery_frac, generation_time, outside_rate
        s = sim(T, S0, I0, R0, infection_rate=infection_rate)

        # This is a column vector where susc[j] is the fraction of group j that is susceptible
        frac_susc = s.get_frac_susceptible()

        # This is a column vector where entry i gives the number of secondary infections resulting from a source
        # infection in group i
        secondary_infections_by_group = np.matmul(infection_rate, frac_susc)

        # This is a column vector whose entry i is the fraction of this group that is infected at the start
        infected_by_group = s.get_frac_infected()

        population_by_group = s.get_pop_count()


        for mgname in popul.metagroup_names():
            # Group indices for those groups in the metagroup
            idx = popul.metagroup_indices(mgname)[0]

            # Vectors giving fraction infected, population size, and number of secondary infections
            # for those groups within the meta-group
            mg_frac_infected = np.array([infected_by_group[i] for i in idx])
            mg_pop = np.array([population_by_group[i] for i in idx])
            mg_secondary_infections = np.array([secondary_infections_by_group[i] for i in idx])

            # Calculate the conditional probability that a person is from group i, given that they are
            # infected and from this meta-group
            mg_num_infected = mg_pop * mg_frac_infected # Number of people infected in each group
            mg_prob = mg_num_infected / np.sum(mg_num_infected)

            secondary_infections[mgname][i] = np.dot(mg_prob, mg_secondary_infections)

    for mgname in popul.metagroup_names():
        plt.plot(days_between_tests,secondary_infections[mgname],label=mgname)

    plt.xlabel('Days between tests')
    plt.ylabel('Effective R0 at the start of the simulation')
    plt.legend()
    plt.axhline(y=1,color='black',linestyle='dashed')

    plt.savefig(out_file)

def usage():
    ''' Print usage message '''
    print('Usage:')
    print('python days_infectious.py [-h] [--help] [--yaml=file.yaml] [--out=file.png] [yaml-overrides]')
    print('--yaml=file.yaml uses file.yaml instead of nominal.yaml')
    print('--outfile=file.png saves the plot to file.png instead of sp22_sim.png')
    print('--help and -h print this message and exit')
    print('yaml-overrides replace parameters in the yaml file and are of the form parameter_name=value')
    print('Example: python days_infectious.py --yaml=nominal_ug.yaml --out=nominal_ug.png pcr_sensitivity=0.8')

if __name__ == "__main__":

    # Default values for arguments not specified in the yaml
    yaml_file = 'nominal.yaml'
    out_file = 'days_infectious.png' # The filename for the plots

    # Parameters from the yaml file to override
    override_params = {}

    for arg in sys.argv[1:]:
        # Handle arguments without an equal sign
        if arg == '-h' or arg == '--help':
            usage()
            exit(0)
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

    main(yaml_file, out_file, **override_params)
