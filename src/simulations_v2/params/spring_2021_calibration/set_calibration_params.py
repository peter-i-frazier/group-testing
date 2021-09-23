import yaml
import numpy as np


with open('nominal.yaml') as f:
    nominal_params = yaml.load(f, Loader=yaml.FullLoader)

"""
parameter update history:
[x] last updated 9/22/2021

Group 1: UG, Greek/athlete
Group 2: UG, non Greek/athlete
Group 3: MBA
Group 4: Grad/prof, non MBA
"""
contact_matrix = np.array([[0.7307, 0.2004, 0.0000, 0.0232],
    [0.0881, 0.2995, 0.0000, 0.0041],
    [0.0000, 0.0000, 0.3939, 0.0265],
    [0.0000, 0.0179, 0.0089, 0.0804]]) # to update
interaction_matrix = np.transpose(contact_matrix)
np.savetxt("interaction_matrix.csv", interaction_matrix)
population_sizes = [3329, 9033, 534, 5227]
arrival_positives = np.array([12, 49, 3, 7])
outside_infections = np.array([4, 2, 3, 1])
outside_infection_rates = np.divide(outside_infections, np.array(population_sizes) * 125)
num_free_infectious_per_arrival_positive = 1.08
secondary_infections_from_arrival_positives = np.dot(interaction_matrix, arrival_positives)
total_num_free_infectious = arrival_positives * num_free_infectious_per_arrival_positive + \
    secondary_infections_from_arrival_positives
initial_prevalence = np.divide(total_num_free_infectious, np.array(population_sizes))
print("free and infectious:", arrival_positives * num_free_infectious_per_arrival_positive)
print("secondary_infections:",  secondary_infections_from_arrival_positives)
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
params_group_1['cases_isolated_per_contact'] = 0.944
params_group_1['cases_quarantined_per_contact'] = 3.233
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
# params_faculty_staff_pre_semester_private = nominal_params.copy()
# params_faculty_staff_pre_semester_private['population_size'] = 10283
# params_faculty_staff_pre_semester_private['test_population_fraction'] = 0
# params_faculty_staff_pre_semester_private['expected_contacts_per_day'] = 10
# params_faculty_staff_pre_semester_private['initial_ID_prevalence'] = 0
# params_faculty_staff_pre_semester_private['test_protocol_QFNR'] = 1 - 0.6 # 0.4
# params_faculty_staff_pre_semester_private['cases_isolated_per_contact'] = 0.255
# params_faculty_staff_pre_semester_private['daily_outside_infection_p'] = 0
# params_faculty_staff_pre_semester_private['_scenario_name'] = 'Faculty + Staff (pre-semester) Parameters, Private'
#
# with open('faculty_staff_pre_semester_private.yaml', 'w') as f:
#     yaml.dump(params_faculty_staff_pre_semester_private, f)
#
# params_faculty_staff_post_movein_private = params_faculty_staff_pre_semester_private.copy()
# params_faculty_staff_post_movein_private['population_size'] = 0
# params_faculty_staff_post_movein_private['test_population_fraction'] = 12.74/130
# params_faculty_staff_post_movein_private['_scenario_name'] = 'Faculty + Staff (post move-in) Parameters, Private'
#
# with open('faculty_staff_post_movein_private.yaml', 'w') as f:
#     yaml.dump(params_faculty_staff_post_movein_private, f)
