import numpy as np
from analysis_helpers import poisson_waiting_function
from subdivide_severity import subdivide_severity

# in reality I think this value will vary a lot in the first few days, 
# and then reach some kind of steady-state, and I'm not sure what makes the most
# sense to use here.  I am setting it to the very pessimistic value of 100% of
# self-reporters are severe, which yields the smallest infectious window size
pct_self_reports_severe = 1

daily_self_report_severe = 0.2
daily_self_report_mild = 0

# avg_infectious_window = (avg time in ID state) + (avg time in Sy state prior to self-reporting)
avg_infectious_window = 4 + pct_self_reports_severe * (1 / daily_self_report_severe) 
if daily_self_report_mild != 0:
    avg_infectious_window += (1 - pct_self_reports_severe) * (1 / daily_self_report_mild)
population_size = 34310

daily_contacts = 25

num_isolations = avg_infectious_window * daily_contacts * 0.4 * 0.026
num_quarantines = max(7 - num_isolations, 0)

prob_severity_given_age = np.array([[0.15, 0.84, 0.01, 0],\
                                    [0.1, 0.77, 0.10, 0.03],\
                                    [0.1, 0.73, 0.10, 0.07],\
                                    [0.1, 0.67, 0.13, 0.10],\
                                    [0.05, 0.55, 0.2, 0.2]])

prob_infection = np.array([0.018, 0.022, 0.029, 0.042, 0.042])
prob_age = np.array([0, 0.85808522, 0.13170574, 0.00878566, 0.00142338])

base_params = {
    'max_time_exposed': 4,
    'exposed_time_function': poisson_waiting_function(max_time=4, mean_time=2),
    
    'max_time_pre_ID': 4,
    'pre_ID_time_function': poisson_waiting_function(max_time=4, mean_time=0),
    
    'max_time_ID': 8,
    'ID_time_function': poisson_waiting_function(max_time=8, mean_time=4),
    
    'max_time_SyID_mild': 18,
    'SyID_mild_time_function': poisson_waiting_function(max_time=18, mean_time=14),
    
    'max_time_SyID_severe': 18,
    'SyID_severe_time_function': poisson_waiting_function(max_time=18, mean_time=14),
    
    'sample_QI_exit_function': (lambda n: np.random.binomial(n, 0.05)),
    'sample_QS_exit_function': (lambda n: np.random.binomial(n, 0.3)),
    
    'exposed_infection_p': 0.026,
    'expected_contacts_per_day': daily_contacts,
    
    'mild_severity_levels': 1,
    'severity_prevalence': subdivide_severity(prob_severity_given_age, prob_infection, prob_age),
    'mild_symptoms_daily_self_report_p': daily_self_report_mild,
    'severe_symptoms_daily_self_report_p': daily_self_report_severe,
    
    'days_between_tests': 300,
    'test_population_fraction': 0,
    
    'test_protocol_QFNR': 0.19,
    'test_protocol_QFPR': 0.005,
    
    'perform_contact_tracing': True,
    'contact_tracing_delay': 2,
    'cases_isolated_per_contact': num_isolations,
    'cases_quarantined_per_contact': num_quarantines,
    'contact_recall': 0.4,

    'pre_ID_state': 'detectable',
    
    'population_size': population_size,
    'initial_E_count': 0,
    'initial_pre_ID_count': 0,
    'initial_ID_count': 0,
    'initial_ID_prevalence': 0.005,
    'initial_SyID_mild_count': 0,
    'initial_SyID_severe_count': 0
}

base_params_testing = base_params.copy()
base_params_testing['days_between_tests'] = 1
base_params_testing['test_population_fraction'] = 0.07
