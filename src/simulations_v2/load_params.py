import yaml
from subdivide_severity import subdivide_severity
import os
import numpy as np
from analysis_helpers import poisson_waiting_function, binomial_exit_function

# upper bound on how far the recursion can go in the yaml-depency tree
MAX_DEPTH=5

# simulation parameters which can be included as yaml-keys but are not required
# they are set to a default value of 0 if not included in the yaml-config
DEFAULT_ZERO_PARAMS = ['initial_E_count',
                        'initial_pre_ID_count',
                        'initial_ID_count',
                        'initial_SyID_mild_count',
                        'initial_SyID_severe_count']


# yaml-keys which share the same key as the simulation parameter, and
# can be copied over one-to-one
COPY_DIRECTLY_YAML_KEYS = ['exposed_infection_p', 'expected_contacts_per_day', 
                            'perform_contact_tracing', 'contact_tracing_delay', 
                            'cases_isolated_per_contact', 'cases_quarantined_per_contact',
            'use_asymptomatic_testing', 'contact_trace_testing_frac', 'days_between_tests',
            'test_population_fraction','test_protocol_QFNR','test_protocol_QFPR',
            'initial_ID_prevalence', 'population_size', 'daily_outside_infection_p'] + \
            DEFAULT_ZERO_PARAMS

def update_sev_prevalence(curr_prevalence_dist, new_asymptomatic_pct):
    new_dist = [new_asymptomatic_pct]
    remaining_mass = sum(curr_prevalence_dist[1:])

    # need to scale so that param_val + x * remaning_mass == 1
    scale = (1 - new_asymptomatic_pct) / remaining_mass
    idx = 1
    while idx < len(curr_prevalence_dist):
        new_dist.append(curr_prevalence_dist[idx] * scale)
        idx += 1
    assert(np.isclose(sum(new_dist), 1))
    return np.array(new_dist)


def load_age_sev_params(param_file):
    with open(param_file) as f:
        age_sev_params = yaml.load(f)

    subparams = age_sev_params['prob_severity_given_age']
    prob_severity_given_age = np.array([
        subparams['agegroup1'],
        subparams['agegroup2'],
        subparams['agegroup3'],
        subparams['agegroup4'],
        subparams['agegroup5'],
    ])

    prob_infection = np.array(age_sev_params['prob_infection_by_age'])
    prob_age = np.array(age_sev_params['prob_age'])
    return subdivide_severity(prob_severity_given_age, prob_infection, prob_age)


# reads multigroup stochastic sim parameters from a yaml config file
def load_multigroup_params(param_file):
    with open(param_file) as f:
        params = yaml.load(f)

    # go through params that point to other directories: start by changing
    # the current working directory so that relative file paths can be parsed
    cwd = os.getcwd()

    nwd = os.path.dirname(os.path.realpath(param_file))
    os.chdir(nwd)

    # read mandatory arguments
    assert('_num_groups' in params)
    num_groups = params['_num_groups']
    assert('_scenario_name' in params)
    scenario_name = params['_scenario_name']


    # load multi-group inputs
    group_params = []
    group_names = []

    # first get individual-group param dicts
    if '_group_configs' not in params:
        raise(Exception("file {} missing _group_configs argument".format(param_file)))

    group_configs = params['_group_configs']

    for i in range(num_groups):
        group_key = '_group_{}'.format(i)
        if group_key not in group_configs:
            raise(Exception("file {} _group_configs argument missing {} specifier".format(param_file, group_key)))
        group_config = group_configs[group_key]

        group_name, group_config_parsed = load_params(additional_params = group_config)
        group_params.append(group_config_parsed)
        group_names.append(group_name)

    # change working-directory back
    os.chdir(cwd)

    # next get inter-group interaction constants
    if '_inter_group_expected_contacts' not in params:
        raise(Exception("file {} missing _inter_group_expected_contacts argument".format(param_file)))

    intergroup_contacts = params['_inter_group_expected_contacts']

    interactions_mtx = np.zeros((num_groups, num_groups))
    for i in range(num_groups):
        interactions_mtx[i,i] = group_params[i]['expected_contacts_per_day']

        key_i = '_group_{}'.format(i)
        if key_i in intergroup_contacts:
            for j in range(num_groups):
                key_j = '_group_{}'.format(j)
                if key_j in intergroup_contacts[key_i]:
                    interactions_mtx[i, j] = intergroup_contacts[key_i][key_j]
            
    return group_params, group_names, interactions_mtx


    




# reads stochastic-simulation parameters from a yaml config file
# supports depence between config files, so that one param file
# can point to another file and params from the pointed-to-file
# are loaded first
def load_params(param_file=None, param_file_stack=[], additional_params = {}):
    if param_file != None:
        assert(len(additional_params) == 0)
        with open(param_file) as f:
            params = yaml.load(f)
        # go through params that point to other directories: start by changing
        # the current working directory so that relative file paths can be parsed
        cwd = os.getcwd()

        nwd = os.path.dirname(os.path.realpath(param_file))
        os.chdir(nwd)

    else:
        params = additional_params
    
    
    if '_inherit_config' in params:
        if len(param_file_stack) >= MAX_DEPTH:
            raise(Exception("yaml config dependency depth exceeded max depth"))
        new_param_file = params['_inherit_config']
        scenario_name, base_params = load_params(new_param_file, param_file_stack + [param_file])
    else:
        scenario_name = None
        base_params = {}

    if '_age_severity_config' in params:
        age_sev_file = params['_age_severity_config']
        severity_dist = load_age_sev_params(age_sev_file)
        base_params['severity_prevalence'] = severity_dist
    else:
        severity_dist = None

    if '_scenario_name' in params:
        scenario_name = params['_scenario_name']
    else:
        # the top-level param-file needs a name
        if len(param_file_stack) == 0:
            raise(Exception("need to specify a _scenario_name value"))
    
    if param_file != None:
        # change working-directory back
        os.chdir(cwd)

    # process the main params loaded from yaml, as well as the additional_params
    # optionally passed as an argument, and store them in base_params 
    for yaml_key, val in params.items():
        # skip the meta-params
        if yaml_key[0] == '_': 
            continue

        if yaml_key == 'ID_time_params':
            assert(len(val)==2)
            mean_time_ID = val[0]
            max_time_ID = val[1]
            base_params['max_time_ID'] = max_time_ID
            base_params['ID_time_function_vals'] = (max_time_ID, mean_time_ID)#poisson_waiting_function(max_time_ID, mean_time_ID)

        elif yaml_key == 'E_time_params':
            assert(len(val) == 2)
            base_params['max_time_exposed'] = val[1]
            base_params['exposed_time_function_vals'] = val#poisson_waiting_function(val[1], val[0])

        elif yaml_key == 'Sy_time_params':
            assert(len(val) == 2)
            base_params['max_time_SyID_mild'] = val[1]
            base_params['SyID_mild_time_function_vals'] = val# = poisson_waiting_function(val[1], val[0])
            base_params['max_time_SyID_severe'] = val[1]
            base_params['SyID_severe_time_function_vals'] = val# poisson_waiting_function(val[1], val[0])

        elif yaml_key == 'asymptomatic_daily_self_report_p':
            base_params['mild_symptoms_daily_self_report_p'] = val

        elif yaml_key == 'symptomatic_daily_self_report_p':
            base_params['severe_symptoms_daily_self_report_p'] = val

        elif yaml_key == 'daily_leave_QI_p':
            base_params['sample_QI_exit_function_val'] = val#binomial_exit_function(val)#(lambda n: np.random.binomial(n, val))

        elif yaml_key == 'daily_leave_QS_p':
            base_params['sample_QS_exit_function_val'] = val#binomial_exit_function(val)#(lambda n: np.random.binomial(n, val))

        elif yaml_key == 'asymptomatic_pct_mult':
            if 'severity_prevalence' not in base_params:
                raise(Exception("encountered asymptomatic_pct_mult with no corresponding severity_dist to modify"))
            new_asymptomatic_p = val * base_params['severity_prevalence'][0]
            base_params['severity_prevalence'] = update_sev_prevalence(base_params['severity_prevalence'], 
                                                                            new_asymptomatic_p)

        elif yaml_key in COPY_DIRECTLY_YAML_KEYS:
            base_params[yaml_key] = val

        else:
            raise(Exception("encountered unknown parameter {}".format(yaml_key)))



    # the pre-ID state is not being used atm so fill it in with some default params here
    if 'max_time_pre_ID' not in base_params:
        base_params['max_time_pre_ID'] = 4
        base_params['pre_ID_time_function'] = poisson_waiting_function(max_time=4, mean_time=0)

    # the following 'initial_count' variables are all defaulted to 0
    for paramname in DEFAULT_ZERO_PARAMS:
        if paramname not in base_params:
            base_params[paramname] = 0

    if 'pre_ID_state' not in base_params:
        base_params['pre_ID_state'] = 'detectable'

    if 'mild_severity_levels' not in base_params:
        base_params['mild_severity_levels'] = 1

    return scenario_name, base_params

