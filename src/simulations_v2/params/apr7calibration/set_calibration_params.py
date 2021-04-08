import yaml

with open('nominal.yaml') as f:
    nominal_params = yaml.load(f, Loader=yaml.FullLoader)

"""
parameter update history:
[x] expected_contacts_per_day updated 4/7
[x] cases_isolated_per_contact -- students: updated 1/11
[x] cases_quarantined_per_contact updated 4/7
[x] cases_isolated_per_contact -- faculty+staff: updated 1/12
[x] test sensitivity updated 4/7


private calibration for students
group 1
"""
params_group_1_pre_semester = nominal_params.copy()
params_group_1_pre_semester['population_size'] = 3533
params_group_1_pre_semester['test_population_fraction'] = 0.0255
params_group_1_pre_semester['expected_contacts_per_day'] = 0.736
params_group_1_pre_semester['cases_isolated_per_contact'] = 1.329
params_group_1_pre_semester['cases_quarantined_per_contact'] = 4.785
params_group_1_pre_semester['initial_ID_prevalence'] = 0
params_group_1_pre_semester['test_protocol_QFNR'] = 1 - 0.6 * 0.974 # 0.4156
params_group_1_pre_semester['daily_outside_infection_p'] = 1.42E-5
params_group_1_pre_semester['_scenario_name'] = 'Group 1 Students (pre-semester) Parameters, Private'

with open('group_1_students_pre_semester_private.yaml', 'w') as f:
    yaml.dump(params_group_1_pre_semester, f)

params_group_1_post_movein_private = params_group_1_pre_semester.copy()
params_group_1_post_movein_private['population_size'] = 0
params_group_1_post_movein_private['test_population_fraction'] = 0.285714
params_group_1_post_movein_private['_scenario_name'] = 'Group 1 Students (post move-in) Parameters, Private'

with open('group_1_students_post_movein_private.yaml', 'w') as f:
    yaml.dump(params_group_1_post_movein_private, f)


"""
private calibration for students
group 2
"""
params_group_2_pre_semester = params_group_1_pre_semester.copy()
params_group_2_pre_semester['population_size'] = 8434
params_group_2_pre_semester['test_population_fraction'] = 0.0255
params_group_2_pre_semester['expected_contacts_per_day'] = 0.148
params_group_2_pre_semester['initial_ID_prevalence'] = 0
params_group_2_pre_semester['daily_outside_infection_p'] = 7.11E-6
params_group_2_pre_semester['_scenario_name'] = 'Group 2 Students (pre-semester) Parameters, Private'

with open('group_2_students_pre_semester_private.yaml', 'w') as f:
    yaml.dump(params_group_2_pre_semester, f)

params_group_2_post_movein_private = params_group_2_pre_semester.copy()
params_group_2_post_movein_private['population_size'] = 0
params_group_2_post_movein_private['test_population_fraction'] = 0.285714
params_group_2_post_movein_private['_scenario_name'] = 'Group 2 Students (post move-in) Parameters, Private'

with open('group_2_students_post_movein_private.yaml', 'w') as f:
    yaml.dump(params_group_2_post_movein_private, f)


"""
private calibration for students
group 3
"""
params_group_3_pre_semester = params_group_1_pre_semester.copy()
params_group_3_pre_semester['population_size'] = 6202
params_group_3_pre_semester['test_population_fraction'] = 0
params_group_3_pre_semester['expected_contacts_per_day'] = 0.067
params_group_3_pre_semester['initial_ID_prevalence'] = 0
params_group_3_pre_semester['daily_outside_infection_p'] = 6.45E-6
params_group_3_pre_semester['_scenario_name'] = 'Group 3 Students (pre-semester) Parameters, Private'

with open('group_3_students_pre_semester_private.yaml', 'w') as f:
    yaml.dump(params_group_3_pre_semester, f)

params_group_3_post_movein_private = params_group_3_pre_semester.copy()
params_group_3_post_movein_private['population_size'] = 0
params_group_3_post_movein_private['test_population_fraction'] = 0.142857
params_group_3_post_movein_private['_scenario_name'] = 'Group 3 Students (post move-in) Parameters, Private'

with open('group_3_students_post_movein_private.yaml', 'w') as f:
    yaml.dump(params_group_3_post_movein_private, f)



"""
private calibration for faculty + staff
"""
params_faculty_staff_pre_semester_private = nominal_params.copy()
params_faculty_staff_pre_semester_private['population_size'] = 10283
params_faculty_staff_pre_semester_private['test_population_fraction'] = 0
params_faculty_staff_pre_semester_private['expected_contacts_per_day'] = 10
params_faculty_staff_pre_semester_private['initial_ID_prevalence'] = 0
params_faculty_staff_pre_semester_private['test_protocol_QFNR'] = 1 - 0.6 # 0.4
params_faculty_staff_pre_semester_private['cases_isolated_per_contact'] = 0.255
params_faculty_staff_pre_semester_private['_scenario_name'] = 'Faculty + Staff (pre-semester) Parameters, Private'

with open('faculty_staff_pre_semester_private.yaml', 'w') as f:
    yaml.dump(params_faculty_staff_pre_semester_private, f)

params_faculty_staff_pre_semester_private = params_faculty_staff_pre_semester_private.copy()
params_faculty_staff_pre_semester_private['population_size'] = 0
params_faculty_staff_pre_semester_private['test_population_fraction'] = 12.74/130
params_faculty_staff_pre_semester_private['_scenario_name'] = 'Faculty + Staff (post move-in) Parameters, Private'

with open('faculty_staff_post_movein_private.yaml', 'w') as f:
    yaml.dump(params_faculty_staff_pre_semester_private, f)
