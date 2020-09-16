import numpy as np

def debug(txt):
    print(txt)

class SurveillanceTesting:
    
    def __init__(self, 
            agents,
            test_FNR,
            mean_test_delay, # params governing the distribution of test-delay-times
            mean_contact_trace_delay,
            test_schedule_proportions, # a dict (a,b,...c) => z
                #the first components specify days of the week to be tested, the final component specifies the proportion of agents
                # who get this test schedule
            non_compliance_params, # a length-2 tuple (a,b) parameterizing a Beta distribution which governs each individual's non-compliance rate
            contact_trace_time_window,
            contact_trace_recall_pct, # how many contacts can we recall (this param also controls double-counting in the simple formula for how many contacts a person has) -- this is really a rate, not a pct
            adaptive_testing_time_window,
            mean_adaptive_testing_delay,
            adaptive_testing_recall_pct,
            contact_inner_products
                ):

        self.adaptive_testing_time_window = adaptive_testing_time_window
        self.mean_adaptive_testing_delay = mean_adaptive_testing_delay
        self.adaptive_testing_recall_pct = adaptive_testing_recall_pct
        self.adaptive_testing_queue = {} 

        self.contact_trace_recall_pct = contact_trace_recall_pct
        self.contact_trace_time_window = contact_trace_time_window
        self.test_FNR = test_FNR
        self.mean_test_delay = mean_test_delay
        self.mean_contact_trace_delay = mean_contact_trace_delay
        self.test_FPR = 0 # currently not using this but writing it here to call attention to that fact
        self.agents = agents

        self.test_reporting = {}
        self.contact_tracing = {}

        # process the test-schedule & compliance arguments
        test_schedules = []
        probs = []
        for schedule, proportion in test_schedule_proportions.items():
            test_schedules.append(schedule)
            probs.append(proportion)

        # sample test schedules for the agents
        normalizer = sum(probs)
        probs = np.array(probs) / normalizer
        test_schedule_idxs = list(range(len(test_schedules)))
        n_agents = len(self.agents)

        test_sched_sample = np.random.choice(test_schedule_idxs, n_agents, replace=True, p=probs)

        # sample non-compliance rates
        non_compliance_sample = np.random.beta(non_compliance_params[0], 
                                                non_compliance_params[1],
                                                n_agents)

        # configure testing params for all agents 
        for agent_id, test_sched_idx, non_compliance in zip(self.agents.keys(),
                                                        test_sched_sample,
                                                        non_compliance_sample):
            test_sched = test_schedules[test_sched_idx]
            self.agents[agent_id].configure_test_params(test_sched, non_compliance)

        # sort closest-contacts for contact tracing
        self.sorted_contacts = {}
        for agent_id, agent in self.agents.items():
            inner_product_vec = np.array(contact_inner_products[agent_id,:].T).flatten()
            ip_idx_list = [(idx, ip) for idx, ip in enumerate(inner_product_vec) if idx != agent_id]
            self.sorted_contacts[agent_id] = [idx for (idx, ip) in sorted(ip_idx_list, key=lambda x: x[1], reverse=True)]



    def step_test(self, t, agent_id_subset=None):
        # run testing logic for day t
        # if you assume FPR = 0 then can use agent_id_subset to only run testing sim code for infected individuals
        if agent_id_subset == None:
            agent_id_subset = self.agents.keys()

        # run the test incorporating time-dependent detectability and non-compliance
        test_results = []
        for agent_id in agent_id_subset:
            schedule = self.agents[agent_id].test_sched
            if any([(t - x) % 7 == 0 for x in schedule]):
                result = self.get_test_result(t, agent_id)
                if result != None:
                    test_results.append((agent_id, result))

        debug("performed {} tests on day {}".format(len(test_results), t))

        return self.step_test_result_reporting(t, test_results)


    def step_test_result_reporting(self, t, test_results):
        # record results for future
        for agent_id, result in test_results:
            delay = np.random.geometric(1/(self.mean_test_delay+1)) - 1
            new_t = t + delay
            if new_t not in self.test_reporting:
                self.test_reporting[new_t] = []
            self.test_reporting[new_t].append((agent_id, result))

        # report results from today
        new_positives = []
        results_for_today = self.test_reporting.pop(t, [])
        for agent_id, result in results_for_today:
            self.agents[agent_id].record_test_result(result)
            if result:
                new_positives.append(result)
                self.agents[agent_id].isolate()

        debug("observed {} test results on day {}; {} were positive".format(len(set(results_for_today)), t, len(set(new_positives))))


        return new_positives


    def step_adaptive_testing(self, t, new_positives):
        # mean of a geometric
        avg_daily_contacts = 1 / self.agents[agent_id].contact_magnitude

        contacts_to_recall = int(avg_daily_contacts * self.adaptive_testing_time_window)

        self.step_followup_testing(t, new_positives,
                                    self.adaptive_testing_queue,
                                    self.mean_adaptive_testing_delay,
                                    contacts_to_recall,
                                    self.adaptive_testing_recall_pct)


    def step_contact_trace(self, t, new_positives):
        # start by getting closest contacts
        # mean of a geometric
        avg_daily_contacts = 1 / self.agents[agent_id].contact_magnitude

        contacts_to_recall = int(avg_daily_contacts * self.contact_trace_time_window)

        self.step_followup_testing(t, new_positives,
                                    self.contact_tracing,
                                    self.mean_contact_trace_delay,
                                    contacts_to_recall,
                                    self.contact_trace_recall_pct)

                
    def step_followup_testing(self, t, new_positives, 
                                followup_test_queue, 
                                mean_delay,
                                contacts_to_recall,
                                recall_pct):
        # add people to followup-test queue (either contact-trace or adaptive testing)
        for agent_id in new_positives:
            if self.agents[agent_id].is_in_isolation:
                continue
            
            contacts = [contact_idx for contact_idx in self.sorted_contacts[agent_id][0:contacts_to_recall] \
                            if np.random.uniform() < recall_pct]
     
            for contact_idx in contacts:
                delay = np.random.geometric(1/(mean_delay + 1)) - 1
                new_t = t + delay
                if new_t not in followup_test_queue:
                    followup_test_queue[new_t] = []
                followup_test_queue[new_t].append(contact_idx)


        # process today's followup tests
        followup_tests_for_today = followup_test_queue.pop(t, [])
        for agent_id in set(followup_tests_for_today):
            self.run_followup_test(t, agent_id)


    def run_followup_test(self, t, agent_id):
        if self.agents[contact_idx].is_in_isolation:
            return
        result = self.get_test_result(t, agent_id)
        if result == None:
            return
        self.agents[agent_id].record_test_result(result)
        if result:
            self.agents[agent_id].isolate()

                

    def get_test_result(self, t, agent_id):
        detectability = self.agents[agent_id].get_detectability(t)
        non_compliance = self.agents[agent_id].get_non_compliance()
        if np.random.uniform() > non_compliance:
            prob_of_detection = detectability * (1-self.test_FNR)
            result = np.random.uniform() < prob_of_detection
        else:
            result = None
        return result





