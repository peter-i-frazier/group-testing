import yaml

with open('nominal.yaml') as f:
    nominal_params = yaml.load(f, Loader=yaml.FullLoader)

"""
public calibration for students
"""
params_students_pre_semester_public = nominal_params.copy()
params_students_pre_semester_public['population_size'] = 7053
params_students_pre_semester_public['test_population_fraction'] = 0.0176
params_students_pre_semester_public['expected_contacts_per_day'] = 10
params_students_pre_semester_public['initial_ID_prevalence'] = 0
params_students_pre_semester_public['_scenario_name'] = 'Students (pre-semester) Parameters, Public'

with open('student_pre_semester_public.yaml', 'w') as f:
    yaml.dump(params_students_pre_semester_public, f)

params_students_post_movein_public = nominal_params.copy()
params_students_post_movein_public['population_size'] = 0
params_students_post_movein_public['test_population_fraction'] = 0.233
params_students_post_movein_public['expected_contacts_per_day'] = 10
params_students_post_movein_public['initial_ID_prevalence'] = 0
params_students_post_movein_public['_scenario_name'] = 'Students (post move-in) Parameters, Public'

with open('students_post_movein_public.yaml', 'w') as f:
    yaml.dump(params_students_post_movein_public, f)
