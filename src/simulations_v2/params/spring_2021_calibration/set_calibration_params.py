import yaml

with open('nominal.yaml') as f:
    nominal_params = yaml.load(f, Loader=yaml.FullLoader)

"""
parameter update history:
[] test frequencies
[] params for employees
[x] others: last updated 9/20/2021

Group 1: UG, Greek/athlete
Group 2: UG, non Greek/athlete
Group 3: MBA
Group 4: Grad/prof, non MBA
"""

"""
private calibration for students
group 1
"""
params_group_1 = nominal_params.copy()
params_group_1['population_size'] = 3329
params_group_1['test_population_fraction'] = 3/7
params_group_1['expected_contacts_per_day'] = 0.9464
params_group_1['cases_isolated_per_contact'] = 0.8729
params_group_1['cases_quarantined_per_contact'] = 3.11
params_group_1['initial_ID_prevalence'] = 14.04 / params_group_1['population_size']
params_group_1['daily_outside_infection_p'] = 1.44E-5
params_group_1['_scenario_name'] = 'Group 1 Students Parameters, Spring 2021'

with open('group_1_students_spring_2021.yaml', 'w') as f:
    yaml.dump(params_group_1, f)


"""
private calibration for students
group 2
"""
params_group_2 = params_group_1.copy()
params_group_2['population_size'] = 9033
params_group_2['test_population_fraction'] = 2/7
params_group_2['expected_contacts_per_day'] = 0.3205
params_group_2['initial_ID_prevalence'] = 56.16 / params_group_2['population_size']
params_group_2['daily_outside_infection_p'] = 1.77E-6
params_group_2['_scenario_name'] = 'Group 2 Students Parameters, Spring 2021'

with open('group_2_students_spring_2021.yaml', 'w') as f:
    yaml.dump(params_group_2, f)


"""
private calibration for students
group 3
"""
params_group_3 = params_group_1.copy()
params_group_3['population_size'] = 534
params_group_3['test_population_fraction'] = 1/7
params_group_3['expected_contacts_per_day'] = 0.5156
params_group_3['initial_ID_prevalence'] = 3.24 / params_group_3['population_size']
params_group_3['daily_outside_infection_p'] = 7.49E-5
params_group_3['_scenario_name'] = 'Group 3 Students Parameters, Spring 2021'

with open('group_3_students_spring_2021.yaml', 'w') as f:
    yaml.dump(params_group_3, f)


"""
private calibration for students
group 4
"""
params_group_4 = params_group_1.copy()
params_group_4['population_size'] = 5227
params_group_4['test_population_fraction'] = 1/7
params_group_4['expected_contacts_per_day'] = 0.0849
params_group_4['initial_ID_prevalence'] = 7.56 / params_group_4['population_size']
params_group_4['daily_outside_infection_p'] = 1.53E-6
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
