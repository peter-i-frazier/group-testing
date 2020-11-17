import numpy as np

def sample_delay(distn):
    values = list(range(len(distn)))
    return np.random.choice(values, 1, p=distn)[0]

class ContactTracing:

    def __init__(self, agents, quarantine_manager, params):
        self.agents=agents
        self.quarantine_manager = quarantine_manager
        self.recall_window = params['ct_recall_window']
        self.delay_dist = params['ct_delay_distribution']
        self.recall_rate = params['ct_recall_rate']

        self.recorded_contacts = {}

        self.ct_start_days = {}
        self.ct_recalled_contacts = {}


    def record_contacts(self, agent_id, contact_ids, day):
        """
        This function should be called once for each free & infected agent 
        on each day they are free and infectious

        The function subsamples the contacts (to model the imperfect recall process)
        and saves the subsampled "recalled" contacts in a dictionary with keys
        corresponding to the day and the infectious agent's ID

        agent_id: the id of the infectious agent
        contact_ids: the ids of all agents with whom they come into contact
        day: the day the interactions occur
        """
        if agent_id not in self.recorded_contacts:
            self.recorded_contacts[agent_id] = {}

        subsampled_contacts = [agent_id for agent_id in contact_ids \
                        if np.random.uniform() <= self.recall_rate]

        self.recorded_contacts[agent_id][day] = subsampled_contacts


    def add_agent_to_trace_queue(self, agent_id, day):
        """
        This function should be called on the day that agent_id tests positive for the virus
        and the plan to trace their contacts is put into place

        The function does two things: 
        1) it gathers their (subsampled) contacts from the recall_window time period preceding
           the current day
        2) it samples the random delay between when a case is identified and when contact tracing
           actually happens and it places the agent_id in a queue of traces to resolve, bucketed into
           the delay-corrected day
        """
        delay = sample_delay(self.delay_dist)
        trace_day = day + delay

        if trace_day not in self.ct_start_days:
            self.ct_start_days[trace_day] = set([agent_id])
        else:
            self.ct_start_days[trace_day].add(agent_id)

        recalled_contacts = set()
        recorded_contacts = self.recorded_contacts.get(agent_id, {})
        for t in range(self.recall_window+1):
            recall_day = day - t
            contacts_on_day = recorded_contacts.get(recall_day, [])
            recalled_contacts = recalled_contacts.union(set(contacts_on_day))

        self.ct_recalled_contacts[agent_id] = recalled_contacts


    def step_trace_queue(self, day):
        """
        This function should be called once a day, 
        ***after agents are added to the contact trace queue for that day***

        The function processes contact traces for the current day

        For now the CT logic is simple - all close contacts are quarantined for 21 days and we never test them
        - but the code is modular and we can handle more sophisticated quarantine dynamics in the future
        """
        agents_to_trace = self.ct_start_days.get(day, [])

        for agent_id in agents_to_trace:
            recalled_contacts = self.ct_recalled_contacts[agent_id]
            
            for recalled_contact_id in recalled_contacts:
                if not self.quarantine_manager.is_agent_isolated(recalled_contact_id):
                    self.quarantine_manager.isolate_agent(recalled_contact_id, day)

                





