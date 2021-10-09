test_policy = {}
test_policy[1] = {'ug_ga_vax_test_frequency': 0/7, 'ug_ga_unvax_test_frequency': 1/7,
            'ug_other_vax_test_frequency': 0/7, 'ug_other_unvax_test_frequency': 1/7,
            'grad_vax_test_frequency': 0/7, 'grad_unvax_test_frequency': 1/7,
            'employee_vax_test_frequency': 1/7, 'employee_unvax_test_frequency': 2/7}

test_policy[2] = {'ug_ga_vax_test_frequency': 0/7, 'ug_ga_unvax_test_frequency': 2/7,
            'ug_other_vax_test_frequency': 0/7, 'ug_other_unvax_test_frequency': 2/7,
            'grad_vax_test_frequency': 0/7, 'grad_unvax_test_frequency': 2/7,
            'employee_vax_test_frequency': 1/7, 'employee_unvax_test_frequency': 2/7}

test_policy[3] = {'ug_ga_vax_test_frequency': 1/7, 'ug_ga_unvax_test_frequency': 1/7,
            'ug_other_vax_test_frequency': 1/7, 'ug_other_unvax_test_frequency': 1/7,
            'grad_vax_test_frequency': 1/7, 'grad_unvax_test_frequency': 1/7,
            'employee_vax_test_frequency': 1/7, 'employee_unvax_test_frequency': 2/7}

test_policy[4] = {'ug_ga_vax_test_frequency': 1/7, 'ug_ga_unvax_test_frequency': 2/7,
            'ug_other_vax_test_frequency': 1/7, 'ug_other_unvax_test_frequency': 2/7,
            'grad_vax_test_frequency': 1/7, 'grad_unvax_test_frequency': 2/7,
            'employee_vax_test_frequency': 1/7, 'employee_unvax_test_frequency': 2/7}

test_policy[5] = {'ug_ga_vax_test_frequency': 2/7, 'ug_ga_unvax_test_frequency': 2/7,
            'ug_other_vax_test_frequency': 2/7, 'ug_other_unvax_test_frequency': 2/7,
            'grad_vax_test_frequency': 2/7, 'grad_unvax_test_frequency': 2/7,
            'employee_vax_test_frequency': 1/7, 'employee_unvax_test_frequency': 2/7}

test_policy[6] = {'ug_ga_vax_test_frequency': 2/7, 'ug_ga_unvax_test_frequency': 2/7,
            'ug_other_vax_test_frequency': 1/7, 'ug_other_unvax_test_frequency': 2/7,
            'grad_vax_test_frequency': 1/7, 'grad_unvax_test_frequency': 2/7,
            'employee_vax_test_frequency': 1/7, 'employee_unvax_test_frequency': 2/7}

test_policy[7] = {'ug_ga_vax_test_frequency': 3/7, 'ug_ga_unvax_test_frequency': 3/7,
            'ug_other_vax_test_frequency': 1/7, 'ug_other_unvax_test_frequency': 2/7,
            'grad_vax_test_frequency': 1/7, 'grad_unvax_test_frequency': 2/7,
            'employee_vax_test_frequency': 1/7, 'employee_unvax_test_frequency': 2/7}

test_policy[8] = {'ug_ga_vax_test_frequency': 0/7, 'ug_ga_unvax_test_frequency': 0/7,
            'ug_other_vax_test_frequency': 0/7, 'ug_other_unvax_test_frequency': 0/7,
            'grad_vax_test_frequency': 0/7, 'grad_unvax_test_frequency': 0/7,
            'employee_vax_test_frequency': 1/7, 'employee_unvax_test_frequency': 2/7}

#0x / wk vax others, 2x / wk for unvax & greek
#0x / wk vax others, 1x / wk for unvax & greek

test_policy[9] = {'ug_ga_vax_test_frequency': 2/7, 'ug_ga_unvax_test_frequency': 2/7,
            'ug_other_vax_test_frequency': 0/7, 'ug_other_unvax_test_frequency': 2/7,
            'grad_vax_test_frequency': 0/7, 'grad_unvax_test_frequency': 2/7,
            'employee_vax_test_frequency': 0/7, 'employee_unvax_test_frequency': 2/7}

test_policy[10] = {'ug_ga_vax_test_frequency': 1/7, 'ug_ga_unvax_test_frequency': 1/7,
            'ug_other_vax_test_frequency': 0/7, 'ug_other_unvax_test_frequency': 1/7,
            'grad_vax_test_frequency': 0/7, 'grad_unvax_test_frequency': 1/7,
            'employee_vax_test_frequency': 0/7, 'employee_unvax_test_frequency': 1/7}


policy_labels = {1: '0x/week vax, 1x/week unvax',
                2: '0x/week vax, 2x/week unvax',
                3: '1x/week vax, 1x/week unvax',
                4: '1x/week vax, 2x/week unvax',
                5: '2x/week vax, 2x/week unvax', 
                6: '1x/week others, 2x/week greek, unvax',
                7: '3x/week greek, 2x/week unvax, 1x/week vax',
                8: 'No testing',
                9: '0x/week others, 2x/week greek, unvax',
                10: '0x/week others, 1x/week greek, unvax'}
