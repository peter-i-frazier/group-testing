import numpy as np

def debug(txt):
    pass
    #print(txt)

def sample_delay(distn):
    values = list(range(len(distn)))
    return np.random.choice(values, 1, p=distn)[0]

class SurveillanceTesting:
    
    def __init__(self, 
            agents,
            test_FNR,
            surveillance_test_delay_distn, # params governing the distribution of test-delay-times
            contact_trace_delay_distn,
            test_schedule_proportions, # a dict (a,b,...c) => z
                #the first components specify days of the week to be tested, the final component specifies the proportion of agents
                # who get this test schedule
            non_compliance_params, # a length-2 tuple (a,b) parameterizing a Beta distribution which governs each individual's non-compliance rate
            adaptive_testing_time_window,
            adaptive_testing_delay_distn,
            adaptive_testing_recall_rate,
            contact_inner_products
                ):

        self.adaptive_testing_time_window = adaptive_testing_time_window
        self.adaptive_testing_delay_distn = adaptive_testing_delay_distn
        self.adaptive_testing_recall_rate = adaptive_testing_recall_rate
        self.adaptive_testing_queue = {} 

        self.test_FNR = test_FNR
        self.surveillance_test_delay_distn = surveillance_test_delay_distn
        self.contact_trace_delay_distn = contact_trace_delay_distn
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

        # record results for future
        for agent_id, result in test_results:
            delay = sample_delay(self.surveillance_test_delay_distn)
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
        agents_for_followup = set([])


        for agent_id in new_positives:
            num_contacts_in_net = int(self.agents[agent_id].get_avg_contacts() * self.adaptive_testing_time_window)
            contacts = [contact_idx for contact_idx in self.sorted_contacts[agent_id][0:num_contacts_in_net] \
                            if np.random.uniform() < self.adaptive_testing_recall_rate]
            agents_for_followup = agents_for_followup.union(set(contacts))


        num_tests_processed = self.step_followup_testing(t, 
                                    self.adaptive_testing_queue,
                                    self.adaptive_testing_delay_distn,
                                    agents_for_followup,
                                    force_test_compliance=True)

        debug("adaptive testing on day {}: submitted {} new agents to adaptive testing queue, performed {} followup tests".format(t, 
                len(agents_for_followup), num_tests_processed))


    def step_contact_trace(self, t, new_positives):
        
        new_agents_for_followup = set([])
        for agent_id in new_positives:
            for contact_set in self.agents[agent_id].previous_contacts:
                new_agents_for_followup = new_agents_for_followup.union(contact_set)
        

        num_tests_processed = self.step_followup_testing(t, 
                                    self.contact_tracing,
                                    self.contact_trace_delay_distn,
                                    new_agents_for_followup,
                                    isolate_for_followup=True)
        
        debug("contact tracing on day {}: submitted {} new agents to the contact trace queue, performed {} followup tests".format(t, 
                len(new_agents_for_followup), num_tests_processed))

                
    def step_followup_testing(self, t,
                                followup_test_queue, 
                                delay_distn,
                                new_agents_for_followup,
                                isolate_for_followup=False,
                                force_test_compliance=False
                                ):
        # add people to followup-test queue (either contact-trace or adaptive testing)
        for agent_id in new_agents_for_followup:
            if self.agents[agent_id].is_in_isolation and not self.agents[agent_id].is_isolated_for_followup:
                continue
            
            if isolate_for_followup:
                self.agents[agent_id].isolate_for_followup()
     
            delay = sample_delay(delay_distn)
            new_t = t + delay
            if new_t not in followup_test_queue:
                followup_test_queue[new_t] = []
            followup_test_queue[new_t].append(agent_id)


        # process today's followup tests
        followup_tests_for_today = followup_test_queue.pop(t, [])
        for agent_id in set(followup_tests_for_today):
            result = self.get_test_result(t, agent_id, force_compliance=force_test_compliance)
            if result == None:
                return
            self.agents[agent_id].record_test_result(result)
            if result:
                self.agents[agent_id].isolate()
            else:
                self.agents[agent_id].remove_from_followup_isolation()

        return len(set(followup_tests_for_today))

                

    def get_test_result(self, t, agent_id, force_compliance=False):
        detectability = self.agents[agent_id].get_detectability(t)
        if force_compliance:
            non_compliance = -1
        else:
            non_compliance = self.agents[agent_id].get_non_compliance()
        if np.random.uniform() > non_compliance:
            prob_of_detection = detectability * (1-self.test_FNR)
            result = np.random.uniform() < prob_of_detection
        else:
            result = None
        return result





