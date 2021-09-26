import yaml
import numpy as np


with open('nominal.yaml') as f:
    nominal_params = yaml.load(f, Loader=yaml.FullLoader)

"""
parameter update history:
[x] contact matrix updated 9/25/2021
[x] contact tracing effectiveness updated 9/25/2021
[x] initial prevalence updated 9/25/2021
yet to update: testing frequency

Group 1: UG, Greek/athlete
Group 2: UG, non Greek/athlete
Group 3: MBA
Group 4: Grad/prof, non MBA
"""
contact_matrix = np.array([[0.6946, 0.1965, 0.0000, 0.0232],
    [0.0744, 0.2852, 0.0000, 0.0041],
    [0.0000, 0.0000, 0.3939, 0.0227],
    [0.0000, 0.0179, 0.0000, 0.0759]])
interaction_matrix = np.transpose(contact_matrix)
np.savetxt("interaction_matrix.csv", interaction_matrix)
population_sizes = [3329, 9033, 534, 5227]
# arrival_positives = np.array([10, 43, 1, 5])
outside_infections = np.array([7, 3, 3, 1])
outside_infection_rates = np.divide(outside_infections, np.array(population_sizes) * 125)
# num_free_infectious_per_arrival_positive = 1.08
# secondary_infections_from_arrival_positives = np.dot(interaction_matrix, arrival_positives)
total_num_free_infectious = np.array([12.9, 55.5, 2.3, 20.3])
initial_prevalence = np.divide(total_num_free_infectious, np.array(population_sizes))
# print("free and infectious (missed by arrival testing):", arrival_positives * num_free_infectious_per_arrival_positive)
# print("free and infectious (non-compliant)",  free_and_infectious)
print("total free and infectious:", total_num_free_infectious)
print("outside_infection_rate:", outside_infection_rates)


"""
Calibration for students
group 1
"""
params_group_1 = nominal_params.copy()
params_group_1['population_size'] = population_sizes[0]
params_group_1['test_population_fraction'] = 3/7
params_group_1['expected_contacts_per_day'] = float(contact_matrix[0,0])
params_group_1['cases_isolated_per_contact'] = 0.854
params_group_1['cases_quarantined_per_contact'] = 3.083
params_group_1['initial_ID_prevalence'] = float(initial_prevalence[0])
params_group_1['daily_outside_infection_p'] = float(outside_infection_rates[0])
params_group_1['_scenario_name'] = 'Group 1 Students Parameters, Spring 2021'

with open('group_1_students_spring_2021.yaml', 'w') as f:
    yaml.dump(params_group_1, f)


"""
Calibration for students
group 2
"""
params_group_2 = params_group_1.copy()
params_group_2['population_size'] = population_sizes[1]
params_group_2['test_population_fraction'] = 2/7
params_group_2['expected_contacts_per_day'] = float(contact_matrix[1,1])
params_group_2['initial_ID_prevalence'] = float(initial_prevalence[1])
params_group_2['daily_outside_infection_p'] = float(outside_infection_rates[1])
params_group_2['_scenario_name'] = 'Group 2 Students Parameters, Spring 2021'

with open('group_2_students_spring_2021.yaml', 'w') as f:
    yaml.dump(params_group_2, f)


"""
Calibration for students
group 3
"""
params_group_3 = params_group_1.copy()
params_group_3['population_size'] = population_sizes[2]
params_group_3['test_population_fraction'] = 1/7
params_group_3['expected_contacts_per_day'] = float(contact_matrix[2,2])
params_group_3['initial_ID_prevalence'] = float(initial_prevalence[2])
params_group_3['daily_outside_infection_p'] = float(outside_infection_rates[2])
params_group_3['_scenario_name'] = 'Group 3 Students Parameters, Spring 2021'

with open('group_3_students_spring_2021.yaml', 'w') as f:
    yaml.dump(params_group_3, f)


"""
Calibration for students
group 4
"""
params_group_4 = params_group_1.copy()
params_group_4['population_size'] = population_sizes[3]
params_group_4['test_population_fraction'] = 1/7
params_group_4['expected_contacts_per_day'] = float(contact_matrix[3,3])
params_group_4['initial_ID_prevalence'] = float(initial_prevalence[3])
params_group_4['daily_outside_infection_p'] = float(outside_infection_rates[3])
params_group_4['_scenario_name'] = 'Group 4 Students Parameters, Spring 2021'

with open('group_4_students_spring_2021.yaml', 'w') as f:
    yaml.dump(params_group_4, f)




"""
private calibration for faculty + staff
"""
params_faculty_staff = nominal_params.copy()
params_faculty_staff['population_size'] = 10283
params_faculty_staff['test_population_fraction'] = 0.146
params_faculty_staff['initial_ID_prevalence'] = 68 / params_faculty_staff['population_size']
params_faculty_staff['cases_isolated_per_contact'] = 0.035
params_faculty_staff['daily_outside_infection_p'] = 0
params_faculty_staff['_scenario_name'] = 'Faculty + Staff Parameters, Spring 2021'

with open('faculty_staff_spring_2021.yaml', 'w') as f:
    yaml.dump(params_faculty_staff, f)
