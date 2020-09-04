from single_dorm_discrete_test_schedule import SingleDormSimulation
from multi_group_simulation import MultiGroupSimulation
from load_params import load_params
import numpy as np

class DormStaggeredTestSimulation(MultiGroupSimulation):
    def __init__(self,
            dorm_population,
            dorm_test_rate,
            high_alert_dorm_test_rate,
            contacts_per_day,
            initial_cases=1,
            safe_days_to_exit_high_alert=3,
            use_default_outside_infection_p=False,
            base_test_schedule=(1, 3),
            high_alert_test_schedule=(1,3,5),
            base_config="/home/jmc678/covid_data/group-testing/src/simulations_v2/"
            "params/june8params/noreopen/nominal_students.yaml"):


        base_schedules = [1,3,5]
        self.sims = [SingleDormSimulation(
                        int(dorm_population/3),
                        dorm_test_rate,
                        high_alert_dorm_test_rate,
                        contacts_per_day / 3,
                        contacts_per_day / 3,
                        base_test_schedule = (base_test_day),
                        high_alert_test_schedule = (1,3,5),
                        base_config=base_config,
                        use_default_outside_infection_p=use_default_outside_infection_p)
                        for base_test_day in base_schedules]
        self.N=3
        self.interaction_matrix = np.matrix([[contacts_per_day / 3] * 3]*3)

    def reset_initial_state(self):
        for sim in self.sims:
            sim.reset()
                
    def update_high_alert_mode(self):
        if any([sim.high_alert_mode for sim in self.sims]):
            for sim in self.sims:
                sim.update_high_alert_mode(force_on=True)
            self.high_alert_mode= True
        else:
            self.high_alert_mode=False

    def run_new_trajectory(self, max_days=112):
        
        days_controlled = []
        free_infected_counts = []
        high_alert_statuses = []
        while len(free_infected_counts) <= max_days:
            self.step()
            self.update_high_alert_mode()

            free_infected = sum([sim.get_free_infected() for sim in self.sims])
            free_infected_counts.append(free_infected)
            high_alert_statuses.append(self.high_alert_mode)

            if free_infected == 0 and not self.high_alert_mode:
                days_controlled.append(True)
            else:
                days_controlled.append(False)

        return days_controlled, sum([sim.get_total_infected() for sim in self.sims]), free_infected_counts, high_alert_statuses
    
    def run_multiple_new_trajectories(self, ntrajectories=100):
        all_trajectories_days_controlled = []
        all_trajectories_total_infected = []
        all_trajectories_free_infected_counts = []
        all_trajectories_high_alert_statuses = []
        for _ in range(ntrajectories):
            self.reset_initial_state()
            d, t, f, h = self.run_new_trajectory()
            all_trajectories_days_controlled.append(d)
            all_trajectories_total_infected.append(t)
            all_trajectories_free_infected_counts.append(f)
            all_trajectories_high_alert_statuses.append(h)
        return all_trajectories_days_controlled, \
                all_trajectories_total_infected, \
                all_trajectories_free_infected_counts, \
                all_trajectories_high_alert_statuses

