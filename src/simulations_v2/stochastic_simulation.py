"""
Implements class containing the stochastic simulation logic outlined in the
following doc:
https://docs.google.com/document/d/18wv_2vcH9tKx1OJ0PpJoI8QZMTSutzX9f44mNKgVS1g/edit#
"""

import numpy as np
import pandas as pd
from math import ceil
from scipy.stats import poisson
import pdb
# from analysis_helpers import binomial_exit_function


def binomial_exit_function(p):
    return (lambda n: np.random.binomial(n, p))


import functools

@functools.lru_cache(maxsize=128)
def poisson_pmf(max_time, mean_time):
    pmf = list()
    for i in range(max_time):
        pmf.append(poisson.pmf(i, mean_time))
    pmf.append(1-np.sum(pmf))
    return np.array(pmf)


def poisson_waiting_function(max_time, mean_time):
    return (lambda n: np.random.multinomial(n, poisson_pmf(max_time, mean_time)))


def poisson_waiting_function2(n, max_time, mean_time):
    return np.random.multinomial(n, poisson_pmf(max_time, mean_time))


class StochasticSimulation:
    def __init__(self, params):

        # Meta-parameters governing the maximum number of days an
        # individual spends in each 'infection' state
        if 'max_time_E' in params:
            self.max_time_E = params['max_time_E']
        else:
            self.max_time_E = params['max_time_exposed']
        self.max_time_pre_ID = params['max_time_pre_ID']
        self.max_time_ID = params['max_time_ID']
        self.max_time_SyID_mild = params['max_time_SyID_mild']
        self.max_time_SyID_severe = params['max_time_SyID_severe']

        # parameters governing distribution over time spent in each
        # of these infection states:
        # Assumptions about the sample_X_times variables:
        # sample_X_times(n) returns a numpy array times of length max_time_X+1
        # such that times[k] is the number of people who stay in state X
        # for k time periods.
        # (However, the length of the corresponding queue arrays will just be max_time_X,
        # not max_time_X + 1)
        # We assume that sum(times) == n

        if 'mean_time_E' in params:
            mean_time = params['mean_time_E']
            self.sample_E_times = poisson_waiting_function(max_time=self.max_time_E, mean_time=mean_time)
        else:
            # self.sample_E_times = params['exposed_time_function']
            self.sample_E_times = poisson_waiting_function(max_time=params['max_time_exposed'], mean_time=params['mean_time_exposed'])

        if 'mean_time_pre_ID' in params:
            mean_time = params['mean_time_pre_ID']
            self.sample_pre_ID_times = poisson_waiting_function(max_time=self.max_time_pre_ID, mean_time=mean_time)
        else:
            # self.sample_pre_ID_times = params['pre_ID_time_function']
            self.sample_pre_ID_times = poisson_waiting_function(max_time=4, mean_time=0)  # copied over from load_params.py

        if 'mean_time_ID' in params:
            mean_time = params['mean_time_ID']
            self.sample_ID_times = poisson_waiting_function(max_time=self.max_time_ID, mean_time=mean_time)
        else:
            # self.sample_ID_times = params['ID_time_function']
            self.sample_ID_times = poisson_waiting_function(params['max_time_ID'], params['mean_time_ID'])

        if 'mean_time_SyID_mild' in params:
            mean_time = params['mean_time_SyID_mild']
            self.sample_SyID_mild_times = poisson_waiting_function(max_time=self.max_time_SyID_mild, mean_time = mean_time)
        else:
            # self.sample_SyID_mild_times = params['SyID_mild_time_function']
            self.sample_SyID_mild_times = poisson_waiting_function(max_time=params['max_time_syID_mild'], mean_time=params['mean_time_syID_mild'])

        if 'mean_time_SyID_severe' in params:
            mean_time = params['mean_time_SyID_severe']
            self.sample_SyID_severe_times = poisson_waiting_function(max_time=self.max_time_SyID_severe, mean_time=mean_time)
        else:
            # self.sample_SyID_severe_times = params['SyID_severe_time_function']
            self.sample_SyID_severe_times = poisson_waiting_function(max_time=params['max_time_syID_mild'], mean_time=params['mean_time_syID_mild'])

        # assumption: sample_QI_exit_count(n) returns a number m <= n
        #             indicating the number of people in the state QI
        #             who exit quarantine, given than n people initially
        #             start there
        # self.sample_QI_exit_count = params['sample_QI_exit_function']
        # self.sample_QS_exit_count = params['sample_QS_exit_function']

        # update so reference in params doesn't include a lambda -sw
        self.sample_QI_exit_count = binomial_exit_function(params['sample_QI_exit_function_param'])
        self.sample_QS_exit_count = binomial_exit_function(params['sample_QS_exit_function_param'])

        # parameters governing distribution over transition out of
        # each infection state
        self.exposed_infection_p = params['exposed_infection_p']
        self.daily_contacts_lambda = params['expected_contacts_per_day']

        # probability that a susceptible individual gets infected from the 'outside' on any given day
        self.daily_outside_infection_p = params['daily_outside_infection_p']

        # mild_severity_levels is the number of severity levels that are contained within the mild class.
        # We assume that they are the first entries in teh severity_prevalence array
        # severity_prevalence is an array that has the distribution of severity levels for any infected patient
        self.mild_severity_levels = params['mild_severity_levels']
        self.severity_prevalence = params['severity_prevalence']
        self.mild_symptoms_p = np.sum(self.severity_prevalence[:self.mild_severity_levels])

        # parameters governing symptomatic daily self reporting
        self.mild_self_report_p = params['mild_symptoms_daily_self_report_p']
        self.severe_self_report_p = params['severe_symptoms_daily_self_report_p']

        # parameters governing test protocol
        use_asymptomatic_testing = params['use_asymptomatic_testing']
        if use_asymptomatic_testing:
            self.days_between_tests = params['days_between_tests']
            self.test_pop_fraction = params['test_population_fraction']
            self.test_QFNR = params['test_protocol_QFNR']
            self.test_QFPR = params['test_protocol_QFPR']
            self.contact_trace_testing_frac = params['contact_trace_testing_frac']
        else:
            self.days_between_tests = 300
            self.test_pop_fraction = 0

            self.test_QFNR = 0.19
            self.test_QFPR = 0.005
            self.contact_trace_testing_frac = 1

        self.perform_contact_tracing = params['perform_contact_tracing']
        self.contact_tracing_delay = params['contact_tracing_delay']

        # new parameters governing contact tracing
        # these can be fractions -- we will round them to integers in the code
        self.cases_isolated_per_contact = params['cases_isolated_per_contact']
        self.cases_quarantined_per_contact = params['cases_quarantined_per_contact']

        # flag governing meaning of the pre-ID state
        self.pre_ID_state = params['pre_ID_state']
        assert(self.pre_ID_state in ['infectious','detectable'])

        # parameters governing initial state of simulation
        self.pop_size = params['population_size']
        self.init_E_count = params['initial_E_count']
        self.init_pre_ID_count = params['initial_pre_ID_count']
        self.init_ID_count = params['initial_ID_count']
        self.init_SyID_mild_count = params['initial_SyID_mild_count']
        self.init_SyID_severe_count = params['initial_SyID_severe_count']

        self.init_ID_prevalence = params['initial_ID_prevalence']
        if 'init_ID_prevalence_stochastic' in params:
            self.init_ID_prevalence_stochastic = params['init_ID_prevalence_stochastic']
        else:
            self.init_ID_prevalence_stochastic = False

        if 'arrival_testing_proportion' in params:
            self.arrival_testing_proportion = params['arrival_testing_proportion']
        else:
            self.arrival_testing_proportion = self.test_pop_fraction

        self.init_S_count = self.pop_size - self.init_E_count - \
            self.init_pre_ID_count - self.init_ID_count - \
            self.init_SyID_mild_count - self.init_SyID_severe_count
        assert(self.init_S_count >= 0)

        # instantiate state variables and relevant simulation variables
        self.reset_initial_state()


    def reset_initial_state(self):
        if self.init_ID_prevalence:
            if self.init_ID_prevalence_stochastic:
                init_ID_count = np.random.binomial(self.pop_size, self.init_ID_prevalence)
            else:
                init_ID_count = ceil(self.pop_size * self.init_ID_prevalence)
        else:
            init_ID_count = self.init_ID_count

        self.S = self.init_S_count + self.init_ID_count - init_ID_count

        # all of the following state vectors have the following convention:
        # state[k] is how many people have k days left to go.
        E_sample = self.sample_E_times(self.init_E_count)
        self.E = E_sample[1:]

        pre_ID_sample = self.sample_pre_ID_times(self.init_pre_ID_count + E_sample[0])
        self.pre_ID = pre_ID_sample[1:]

        ID_sample = self.sample_ID_times(init_ID_count + pre_ID_sample[0])
        self.ID = ID_sample[1:]

        additional_mild = np.random.binomial(ID_sample[0], self.mild_symptoms_p)
        additional_severe = ID_sample[0] - additional_mild

        SyID_mild_sample = self.sample_SyID_mild_times(self.init_SyID_mild_count + additional_mild)
        self.SyID_mild = SyID_mild_sample[1:]

        SyID_severe_sample = self.sample_SyID_severe_times(self.init_SyID_severe_count + additional_severe)
        self.SyID_severe = SyID_severe_sample[1:]

        # contact_trace_queue[k] are the number of quarantined individuals who have k
        # days remaining until the results from their contact trace comes in
        self.contact_trace_queue = [0] * (self.contact_tracing_delay + 1)

        self.QS = 0
        self.QI = 0

        self.QI_mild = 0
        self.QI_severe = 0

        self.R = SyID_mild_sample[0] + SyID_severe_sample[0]
        self.R_mild = SyID_mild_sample[0]
        self.R_severe = SyID_severe_sample[0]

        self.cumulative_outside_infections = 0
        var_labels = self.get_state_vector_labels()
        self.sim_df = pd.DataFrame(columns=var_labels)
        self._append_sim_df()
        self.current_day = 0
        self.last_test_day = -1
        self.new_QI_from_last_test = 0
        self.new_QS_from_last_test = 0
        self.new_QI_from_self_reports = 0

    def run_new_trajectory(self, T):
        self.reset_initial_state()
        for _ in range(T):
            self.step()

        for i in range(len(self.severity_prevalence)):
            if i < self.mild_severity_levels:
                # Mild severity
                self.sim_df['severity_'+str(i)] = self.sim_df['cumulative_mild'] * (self.severity_prevalence[i] / self.mild_symptoms_p)
            else:
                # Severe symptoms
                self.sim_df['severity_'+str(i)] = self.sim_df['cumulative_severe'] * (self.severity_prevalence[i] / (1 - self.mild_symptoms_p))

        return self.sim_df

    def step_contact_trace(self, new_QI):
        """ resolve contact traces at the front of the queue and add new QIs to the back
        of the contact trace queue"""

        # update the contact trace queue
        self.contact_trace_queue[self.contact_tracing_delay] += new_QI
        resolve_today_QI = self.contact_trace_queue[0]
        self._shift_contact_queue()

        # compute how many cases we find
        #total_contacts = int(resolve_today_QI * self.contact_trace_infectious_window \
        #                                * self.daily_contacts_lambda)
        #total_contacts_traced = np.random.binomial(total_contacts, self.contact_tracing_c)
        #total_cases_isolated = np.random.binomial(total_contacts_traced, self.exposed_infection_p)
        #total_contacts_quarantined = min(self.S, total_contacts_traced - total_cases_isolated)

        total_contacts_quarantined = min(self.S, int(self.cases_quarantined_per_contact * resolve_today_QI))
        # add susceptible people to the quarantine state
        self.S = self.S - total_contacts_quarantined
        self.QS = self.QS + total_contacts_quarantined

        total_cases_isolated = int(self.cases_isolated_per_contact * resolve_today_QI)

        # trace these cases across E, pre-ID and ID states

        initial_isolations = total_cases_isolated

        leave_E = int(min(sum(self.E), total_cases_isolated))
        self._trace_E_queue(leave_E)
        total_cases_isolated -= leave_E

        leave_pre_ID = min(sum(self.pre_ID), total_cases_isolated)
        self._trace_pre_ID_queue(leave_pre_ID)
        total_cases_isolated -= leave_pre_ID

        leave_ID = min(sum(self.ID), total_cases_isolated)
        self._trace_ID_queue(leave_ID)
        total_cases_isolated -= leave_ID

        leave_SyID_severe = min(sum(self.SyID_severe), total_cases_isolated)
        self._trace_SyID_severe_queue(leave_SyID_severe)
        total_cases_isolated -= leave_SyID_severe

        leave_SyID_mild = min(sum(self.SyID_mild), total_cases_isolated)
        self._trace_SyID_mild_queue(leave_SyID_mild)
        total_cases_isolated -= leave_SyID_mild

        #print("initial isolations: {}, final isolations: {}".format(initial_isolations, total_cases_isolated))


    def _trace_SyID_severe_queue(self, leave_SyID_severe):
        assert(leave_SyID_severe <= sum(self.SyID_severe))
        self.QI = self.QI + leave_SyID_severe
        self.QI_severe += leave_SyID_severe
        idx = self.max_time_SyID_severe - 1
        while leave_SyID_severe > 0:
            leave_SyID_severe_at_idx = min(self.SyID_severe[idx], leave_SyID_severe)
            self.SyID_severe[idx] -= leave_SyID_severe_at_idx
            leave_SyID_severe -= leave_SyID_severe_at_idx
            idx -= 1


    def _trace_SyID_mild_queue(self, leave_SyID_mild):
        assert(leave_SyID_mild <= sum(self.SyID_mild))
        self.QI = self.QI + leave_SyID_mild
        self.QI_mild += leave_SyID_mild
        idx = self.max_time_SyID_mild - 1
        while leave_SyID_mild > 0:
            leave_SyID_mild_at_idx = min(self.SyID_mild[idx], leave_SyID_mild)
            self.SyID_mild[idx] -= leave_SyID_mild_at_idx
            leave_SyID_mild -= leave_SyID_mild_at_idx
            idx -= 1
   

    def _trace_E_queue(self, leave_E):
        assert(leave_E <= sum(self.E))
        self.QI = self.QI + leave_E
        leave_E_mild = np.random.binomial(leave_E, self.mild_symptoms_p)
        leave_E_severe = leave_E - leave_E_mild
        self.QI_mild += leave_E_mild
        self.QI_severe += leave_E_severe
        idx = self.max_time_E - 1
        while leave_E > 0:
            leave_E_at_idx = min(self.E[idx], leave_E)
            self.E[idx] -= leave_E_at_idx
            leave_E -= leave_E_at_idx
            idx -= 1

    def _trace_pre_ID_queue(self, leave_pre_ID):
        assert(leave_pre_ID <= sum(self.pre_ID))
        self.QI = self.QI + leave_pre_ID
        leave_pre_ID_mild = np.random.binomial(leave_pre_ID, self.mild_symptoms_p)
        leave_pre_ID_severe = leave_pre_ID - leave_pre_ID_mild
        self.QI_mild += leave_pre_ID_mild
        self.QI_severe += leave_pre_ID_severe
        idx = self.max_time_pre_ID - 1
        while leave_pre_ID > 0:
            leave_pre_ID_at_idx = min(self.pre_ID[idx], leave_pre_ID)
            self.pre_ID[idx] -= leave_pre_ID_at_idx
            leave_pre_ID -= leave_pre_ID_at_idx
            idx -= 1

    def _trace_ID_queue(self, leave_ID):
        assert(leave_ID <= sum(self.ID))
        self.QI = self.QI + leave_ID
        leave_ID_mild = np.random.binomial(leave_ID, self.mild_symptoms_p)
        leave_ID_severe = leave_ID - leave_ID_mild
        self.QI_mild += leave_ID_mild
        self.QI_severe += leave_ID_severe
        idx = self.max_time_ID - 1
        while leave_ID > 0:
            leave_ID_at_idx = min(self.ID[idx], leave_ID)
            self.ID[idx] -= leave_ID_at_idx
            leave_ID -= leave_ID_at_idx
            idx -= 1

    def _shift_contact_queue(self):
        idx = 0
        while idx <= self.contact_tracing_delay - 1:
            self.contact_trace_queue[idx] = self.contact_trace_queue[idx+1]
            idx += 1
        self.contact_trace_queue[self.contact_tracing_delay] = 0

    def run_test(self):
        """
        Execute one step of the testing logic.
        """

        # infectious_test_pop = free_infectious * self.test_pop_fraction
        # fluid_new_QI = infectious_test_pop * (1 - self.test_QFNR)

        # the probability that a free infected individual is quarantined
        # on this round of testing.

        # If arrival testing is specified, uses that on frist day, otherwise
        # all test_pop_fractions default to self.test_pop_fraction (configured in
        # param reads -SW)
        if self.current_day == 0:
            test_pop_fraction = self.arrival_testing_proportion
        else:
            test_pop_fraction = self.test_pop_fraction

        new_QI_p = test_pop_fraction * (1 - self.test_QFNR)

        # sample the number of free infected people who end up quarantined
        new_QI_from_ID = np.random.binomial(self.ID, new_QI_p)
        new_QI_from_ID_mild = np.random.binomial(new_QI_from_ID, self.mild_symptoms_p)
        new_QI_from_ID_severe = new_QI_from_ID - new_QI_from_ID_mild
        new_QI_from_SyID_mild = np.random.binomial(self.SyID_mild, new_QI_p)
        new_QI_from_SyID_severe = np.random.binomial(self.SyID_severe, new_QI_p)

        # update counts in relevant states
        self.ID = self.ID - new_QI_from_ID
        self.SyID_mild = self.SyID_mild - new_QI_from_SyID_mild
        self.SyID_severe = self.SyID_severe - new_QI_from_SyID_severe

        new_QI = sum(new_QI_from_ID) + sum(new_QI_from_SyID_mild) + sum(new_QI_from_SyID_severe)
        new_QI_mild = sum(new_QI_from_ID_mild) + sum(new_QI_from_SyID_mild)
        new_QI_severe = sum(new_QI_from_ID_severe) + sum(new_QI_from_SyID_severe)

        # do the above for pre-ID state, if it is detectable
        if self.pre_ID_state == 'detectable':
            new_QI_from_pre_ID = np.random.binomial(self.pre_ID, new_QI_p)
            new_QI_from_pre_ID_mild = np.random.binomial(new_QI_from_pre_ID, self.mild_symptoms_p)
            new_QI_from_pre_ID_severe = new_QI_from_pre_ID - new_QI_from_pre_ID_mild
            self.pre_ID = self.pre_ID - new_QI_from_pre_ID
            new_QI += sum(new_QI_from_pre_ID)
            new_QI_mild += sum(new_QI_from_pre_ID_mild)
            new_QI_severe += sum(new_QI_from_pre_ID_severe)

        # add to QI individuals from E, and from pre-ID (if state is 'infectious'), using
        # the false-positive rate for undetectable individuals
        new_QI_undetectable_p = test_pop_fraction * self.test_QFPR

        new_QI_from_E = np.random.binomial(self.E, new_QI_undetectable_p)
        new_QI_from_E_mild = np.random.binomial(new_QI_from_E, self.mild_symptoms_p)
        new_QI_from_E_severe = new_QI_from_E - new_QI_from_E_mild
        self.E = self.E - new_QI_from_E
        new_QI += sum(new_QI_from_E)
        new_QI_mild += sum(new_QI_from_E_mild)
        new_QI_severe += sum(new_QI_from_E_severe)

        if self.pre_ID_state == 'infectious':
            new_QI_from_pre_ID = np.random.binomial(self.pre_ID, new_QI_undetectable_p)
            new_QI_from_pre_ID_mild = np.random.binomial(new_QI_from_pre_ID, self.mild_symptoms_p)
            new_QI_from_pre_ID_severe = new_QI_from_pre_ID - new_QI_from_pre_ID_mild
            self.pre_ID = self.pre_ID - new_QI_from_pre_ID
            new_QI += sum(new_QI_from_pre_ID)
            new_QI_mild += sum(new_QI_from_pre_ID_mild)
            new_QI_severe += sum(new_QI_from_pre_ID_severe)

        # add to QS individuals from S, due to false positives
        new_QS_p = test_pop_fraction * self.test_QFPR
        # sample number of free susceptible people who become quarantined
        new_QS_from_S = np.random.binomial(self.S, new_QS_p)
        self.S = self.S - new_QS_from_S

        # update QS and QI
        self.QS = self.QS + new_QS_from_S
        self.QI = self.QI + new_QI
        self.QI_mild += new_QI_mild
        self.QI_severe += new_QI_severe

        self.new_QI_from_last_test = new_QI
        self.new_QS_from_last_test = new_QS_from_S

        return new_QI


    def isolate_self_reports(self):
        mild_self_reports = np.random.binomial(self.SyID_mild, self.mild_self_report_p)
        self.SyID_mild = self.SyID_mild - mild_self_reports
        new_QI = sum(mild_self_reports)
        new_QI_mild = sum(mild_self_reports)

        severe_self_reports = np.random.binomial(self.SyID_severe, self.severe_self_report_p)
        self.SyID_severe = self.SyID_severe - severe_self_reports
        new_QI += sum(severe_self_reports)
        new_QI_severe = sum(severe_self_reports)

        self.QI = self.QI + new_QI
        self.QI_mild += new_QI_mild
        self.QI_severe += new_QI_severe

        self.new_QI_from_self_reports = new_QI
        return new_QI



    def step(self):
        """ simulate a single day in the progression of the disease """

        new_QI = 0
        new_contact_traces = 0
        # do testing logic first
        if self.current_day - self.last_test_day >= self.days_between_tests:
            self.last_test_day = self.current_day
            new_QI += self.run_test()
            new_contact_traces += int(self.contact_trace_testing_frac * new_QI)


        # resolve symptomatic self-reporting
        new_self_reports = self.isolate_self_reports()
        new_QI += new_self_reports
        new_contact_traces += new_self_reports

        # do contact tracing
        if self.perform_contact_tracing:
            self.step_contact_trace(new_contact_traces)

        # simulate number of contacts between free infectious & free susceptible:
        free_infectious = 0

        if self.pre_ID_state == 'infectious':
            free_infectious += sum(self.pre_ID)

        free_infectious += sum(self.ID) + sum(self.SyID_mild) + sum(self.SyID_severe)

        free_susceptible = self.S
        free_tot = free_infectious + free_susceptible + self.R + sum(self.E) 

        if self.pre_ID_state == 'detectable':
            free_tot += sum(self.pre_ID)

        if free_tot == 0:
            poisson_param = 0
        else:
            poisson_param = free_infectious * self.daily_contacts_lambda * free_susceptible / free_tot
        n_contacts = np.random.poisson(poisson_param)
        #n_contacts = int(free_infectious * free_susceptible / free_tot * np.random.geometric(1/self.daily_contacts_lambda))

        # sample number of new E cases from 'inside' contacts
        new_E_from_inside = min(np.random.binomial(n_contacts, self.exposed_infection_p), self.S)
        
        # sample number of new E cases from 'outside' infection
        new_E_from_outside = np.random.binomial(self.S - new_E_from_inside, self.daily_outside_infection_p)
        self.cumulative_outside_infections += new_E_from_outside

        new_E = new_E_from_inside + new_E_from_outside
        self.S -= new_E

        # update E queue and record new pre-ID cases
        new_E_times = self.sample_E_times(new_E)
        new_pre_ID = self.E[0] + new_E_times[0]
        self._shift_E_queue()
        self.E = self.E + new_E_times[1:]

        # sample times of new pre-ID cases / update pre-ID queue/ record new ID cases
        new_pre_ID_times = self.sample_pre_ID_times(new_pre_ID)
        new_ID = self.pre_ID[0] + new_pre_ID_times[0]
        self._shift_pre_ID_queue()
        self.pre_ID = self.pre_ID + new_pre_ID_times[1:]


        # sample times of new ID cases / update ID queue/ record new SyID cases
        new_ID_times = self.sample_ID_times(new_ID)
        new_SyID = self.ID[0] + new_ID_times[0]
        self._shift_ID_queue()
        self.ID = self.ID + new_ID_times[1:]

        # decompose new_SyID into mild and severe
        new_SyID_mild = np.random.binomial(new_SyID, self.mild_symptoms_p)
        new_SyID_severe = new_SyID - new_SyID_mild

        # samples times of new SyID mild cases/ update mild queue/ record new R cases
        new_SyID_mild_times = self.sample_SyID_mild_times(new_SyID_mild)
        new_R_from_mild = self.SyID_mild[0] + new_SyID_mild_times[0]
        self._shift_SyID_mild_queue()
        self.SyID_mild = self.SyID_mild + new_SyID_mild_times[1:]

        # same as above, but for the severe symptom queue
        new_SyID_severe_times = self.sample_SyID_severe_times(new_SyID_severe)
        new_R_from_severe = self.SyID_severe[0] + new_SyID_severe_times[0]
        self._shift_SyID_severe_queue()
        self.SyID_severe = self.SyID_severe + new_SyID_severe_times[1:]


        # sample number of people who leave quarantine-I/ resolve new R cases
        leave_QI = self.sample_QI_exit_count(self.QI)
        if leave_QI == 0:
            leave_QI_mild = 0
            leave_QI_severe = 0
        else:
            leave_QI_mild = min(np.random.binomial(leave_QI, self.QI_mild / self.QI), self.QI_mild)
            leave_QI_severe = min(leave_QI - leave_QI_mild, self.QI_severe)
            leave_QI = leave_QI_mild + leave_QI_severe
        self.QI -= leave_QI
        self.QI_mild -= leave_QI_mild
        self.QI_severe -= leave_QI_severe
        self.R += leave_QI + new_R_from_mild + new_R_from_severe
        self.R_mild += leave_QI_mild + new_R_from_mild
        self.R_severe += leave_QI_severe + new_R_from_severe
        leave_QS = self.sample_QS_exit_count(self.QS)
        self.QS -= leave_QS
        self.S += leave_QS

        self._append_sim_df()

        self.current_day += 1


    ## add new_E people to the infections queue from the S queue.
    ## this function is written to support the companion multi-group simulation
    def add_new_infections(self, new_E):

        new_E = min(self.S, new_E)

        self.S = self.S - new_E

        # in theory it is possible for someone to go from new_E to R in a single step,
        # so we have to pass through all the states...
        new_E_times = self.sample_E_times(new_E)
        new_pre_ID = new_E_times[0]
        self.E = self.E + new_E_times[1:]

        # sample times of new pre-ID cases / update pre-ID queue/ record new ID cases
        new_pre_ID_times = self.sample_pre_ID_times(new_pre_ID)
        new_ID =  new_pre_ID_times[0]
        self.pre_ID = self.pre_ID + new_pre_ID_times[1:]


        # sample times of new ID cases / update ID queue/ record new SyID cases
        new_ID_times = self.sample_ID_times(new_ID)
        new_SyID = new_ID_times[0]
        self.ID = self.ID + new_ID_times[1:]

        # decompose new_SyID into mild and severe
        new_SyID_mild = np.random.binomial(new_SyID, self.mild_symptoms_p)
        new_SyID_severe = new_SyID - new_SyID_mild

        # samples times of new SyID mild cases/ update mild queue/ record new R cases
        new_SyID_mild_times = self.sample_SyID_mild_times(new_SyID_mild)
        new_R_from_mild = new_SyID_mild_times[0]
        self.SyID_mild = self.SyID_mild + new_SyID_mild_times[1:]

        # same as above, but for the severe symptom queue
        new_SyID_severe_times = self.sample_SyID_severe_times(new_SyID_severe)
        new_R_from_severe = new_SyID_severe_times[0]
        self.SyID_severe = self.SyID_severe + new_SyID_severe_times[1:]

        self.R += new_R_from_mild + new_R_from_severe
        self.R_mild += new_R_from_mild
        self.R_severe += new_R_from_severe


    def _append_sim_df(self):
        self.generate_cumulative_stats()
        data = self.get_current_state_vector()
        labels = self.get_state_vector_labels()
        new_row_df = pd.DataFrame([data], columns=labels)
        self.sim_df = self.sim_df.append(new_row_df, ignore_index=True)
        # print(sum(data), sum(data[-1*(len(self.severity_prevalence)+2):]))
        # print(labels[-1*(len(self.severity_prevalence)+2):])
        assert(self.QI_mild + self.QI_severe == self.QI)
        assert(self.R_mild + self.R_severe == self.R)
        assert(min(self.QI_mild, self.QI_severe, self.R_mild, self.R_severe) >= 0)
        if abs(sum(data) - sum(data[-3:]) - self.pop_size) > 0.0001:
            raise(Exception("population has shrunk"))
        if np.sum(data < 0) > 0:
            raise(Exception("negative category size"))


    def _shift_E_queue(self):
        idx = 0
        while idx <= self.max_time_E - 2:
            self.E[idx] = self.E[idx+1]
            idx += 1
        self.E[self.max_time_E - 1] = 0

    def _shift_pre_ID_queue(self):
        idx = 0
        while idx <= self.max_time_pre_ID - 2:
            self.pre_ID[idx] = self.pre_ID[idx+1]
            idx += 1
        self.pre_ID[self.max_time_pre_ID - 1] = 0

    def _shift_ID_queue(self):
        idx = 0
        while idx <= self.max_time_ID - 2:
            self.ID[idx] = self.ID[idx+1]
            idx += 1
        self.ID[self.max_time_ID - 1] = 0

    def _shift_SyID_mild_queue(self):
        idx = 0
        while idx <= self.max_time_SyID_mild - 2:
            self.SyID_mild[idx] = self.SyID_mild[idx+1]
            idx += 1
        self.SyID_mild[self.max_time_SyID_mild - 1] = 0

    def _shift_SyID_severe_queue(self):
        idx = 0
        while idx <= self.max_time_SyID_severe - 2:
            self.SyID_severe[idx] = self.SyID_severe[idx+1]
            idx += 1
        self.SyID_severe[self.max_time_SyID_severe - 1] = 0


    def get_current_state_vector(self):
        return np.concatenate([
            [self.S], [self.QS], [self.QI], [self.R],
            self.E, self.pre_ID, self.ID, self.SyID_mild, self.SyID_severe,
            [self.cumulative_mild], [self.cumulative_severe], [self.cumulative_outside_infections]
            ])

    def get_state_vector_labels(self):
        return ['S', 'QS', 'QI', 'R'] + \
                ['E_{}'.format(x) for x in range(self.max_time_E)] + \
                ['pre_ID_{}'.format(x) for x in range(self.max_time_pre_ID)] + \
                ['ID_{}'.format(x) for x in range(self.max_time_ID)] + \
                ['SyID_mild_{}'.format(x) for x in range(self.max_time_SyID_mild)] + \
                ['SyID_severe_{}'.format(x) for x in range(self.max_time_SyID_severe)] + \
                ['cumulative_mild', 'cumulative_severe', 'cumulative_outside_infections']


    def generate_cumulative_stats(self):
        self.cumulative_mild = self.QI_mild + sum(self.SyID_mild) + self.R_mild
        self.cumulative_severe = self.QI_severe + sum(self.SyID_severe) + self.R_severe

        # self.severity = list()
        # for i in range(self.mild_severity_levels):
        #     self.severity.append((self.severity_prevalence[i] / self.mild_symptoms_p ) * self.cumulative_mild)
        # for i in range(len(self.severity_prevalence) - self.mild_severity_levels):
        #     self.severity.append((self.severity_prevalence[i + self.mild_severity_levels]) / (1 - self.mild_symptoms_p) * self.cumulative_severe)

    def update_severity_levels(self):
        for i in range(len(self.severity_prevalence)):
            if i < self.mild_severity_levels:
                self.sim_df['severity_'+str(i)] = self.sim_df['cumulative_mild'] * (self.severity_prevalence[i] / self.mild_symptoms_p)
            else:
                self.sim_df['severity_'+str(i)] = self.sim_df['cumulative_severe'] * (self.severity_prevalence[i] / (1 - self.mild_symptoms_p))

