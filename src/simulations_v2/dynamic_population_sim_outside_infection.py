from multi_group_simulation import MultiGroupSimulation
from stochastic_simulation import StochasticSimulation
import numpy as np

class DynamicPopulationSim:

    def __init__(self, 
            movein_free_params, # initial simulation parameters for free group during movein period
            movein_selfiso_params, # initial simulation params for self-isolate group during movein period
            post_movein_params, # sim params for post-movein period -- the initialization params here will be irrelevant, b/c they
                                # will be copied over from the state of the move-in simulation at the end of the move-in period
            movein_contact_matrix, # 2x2 interaction matrix for movein period -- order of groups is [free, selfiso]
            movein_time_horizon, # number of time periods in the movein period
            free_group_dynamics, # dictionary mapping time periods t -> 
                                         # (dictionary specifying how many people join each compartment on day t in the free group)
            selfiso_group_dynamics, # dictionary mapping time periods t -> 
                                         # (dictionary specifying how many people join each compartment on day t in the self-isolation group)
                                         ):

        self.movein_sim = MultiGroupSimulation([movein_free_params, movein_selfiso_params], movein_contact_matrix)
        self.post_movein_sim = None

        self.movein_time_horizon = movein_time_horizon

        self.free_group_dynamics = free_group_dynamics
        self.selfiso_group_dynamics = selfiso_group_dynamics

        self.post_movein_params = post_movein_params

        self.current_t = 0
        self.current_week = 0

    def get_sim_obj(self):
        if self.current_t >= self.movein_time_horizon:
            return self.post_movein_sim
        else:
            return self.movein_sim.sims[0]


    def add_new_infections(self, new_E):
        sim_obj = self.get_sim_obj()
        sim_obj.add_new_infections(new_E)


    def step(self):
        outside_df = self.post_movein_params['outside_infection_p_array']
        self.current_week = max(round((self.current_t - 18)/7),0)
        
        if self.current_week in outside_df['week_since_sem_start'].values:
            self.post_movein_params['daily_outside_infection_p'] = outside_df[outside_df['week_since_sem_start'] == self.current_week]['weekly_outside_cases'].values[0]
                
        else:
            self.post_movein_params['daily_outside_infection_p'] = 0
            
        print(self.current_week, self.post_movein_params['daily_outside_infection_p'])
        
        if self.current_t >= self.movein_time_horizon:
            self.post_movein_sim.step()
        else:
            #self.update_movein_populations()
            # update move-in populations
            self.update_populations(self.movein_sim.sims[0], self.free_group_dynamics.get(self.current_t, {}))
            self.update_populations(self.movein_sim.sims[1], self.selfiso_group_dynamics.get(self.current_t, {}))
            self.movein_sim.step()

        self.current_t += 1
        if self.current_t == self.movein_time_horizon:
            self.initialize_post_movein_sim()
        
        

    def initialize_post_movein_sim(self):
        assert(self.current_t == self.movein_time_horizon)
        self.post_movein_sim = StochasticSimulation(self.post_movein_params)

        pop_size = 0 

        # give vars shorter names to save typing effort...
        sims = self.movein_sim.sims
        new_sim = self.post_movein_sim

        new_sim.S = sims[0].S + sims[1].S
        pop_size += new_sim.S
        
        new_sim.E = sims[0].E + sims[1].E
        pop_size += sum(new_sim.E)
        
        new_sim.pre_ID = sims[0].pre_ID + sims[1].pre_ID
        pop_size += sum(new_sim.pre_ID)
        
        new_sim.ID = sims[0].ID + sims[1].ID
        pop_size += sum(new_sim.ID)

        new_sim.SyID_mild = sims[0].SyID_mild + sims[1].SyID_mild
        pop_size += sum(new_sim.SyID_mild)

        new_sim.SyID_severe = sims[0].SyID_severe + sims[1].SyID_severe
        pop_size += sum(new_sim.SyID_severe)

        new_sim.R = sims[0].R + sims[1].R
        pop_size += new_sim.R

        new_sim.QI = sims[0].QI + sims[1].QI
        pop_size += new_sim.QI

        new_sim.QS = sims[0].QS + sims[1].QS
        pop_size += new_sim.QS

        new_sim.R_mild = sims[0].R_mild + sims[1].R_mild
        new_sim.R_severe = sims[0].R_severe + sims[1].R_severe

        new_sim.QI_mild = sims[0].QI_mild + sims[1].QI_mild
        new_sim.QI_severe = sims[0].QI_severe + sims[1].QI_severe

        new_sim.pop_size = pop_size


    # the logic for this function follows closely what is done in reset_initial_state()
    # for the main population-level simulation code
    def update_populations(self, sim, pop_updates):
        
        pop_increase = 0

        new_S = pop_updates.get('S', 0)
        pop_increase += new_S
        sim.S += new_S

        new_E = pop_updates.get('E', 0)
        pop_increase += new_E
        E_sample = sim.sample_E_times(new_E)
        sim.E = sim.E + E_sample[1:]
        E_carryover = E_sample[0]

        new_pre_ID = pop_updates.get('pre_ID', 0)
        pop_increase += new_pre_ID
        pre_ID_sample = sim.sample_pre_ID_times(new_pre_ID + E_carryover)
        sim.pre_ID = sim.pre_ID + pre_ID_sample[1:]
        pre_ID_carryover = pre_ID_sample[0]

        new_ID = pop_updates.get('ID', 0)
        pop_increase += new_ID
        ID_sample = sim.sample_ID_times(new_ID + pre_ID_carryover)
        sim.ID = sim.ID + ID_sample[1:]
        ID_carryover = ID_sample[0]

        ID_carryover_mild = np.random.binomial(ID_carryover, sim.mild_symptoms_p)
        ID_carryover_severe = ID_carryover - ID_carryover_mild

        new_SyID_mild = pop_updates.get('SyID_mild', 0)
        pop_increase += new_SyID_mild
        SyID_mild_sample = sim.sample_SyID_mild_times(new_SyID_mild + ID_carryover_mild)
        sim.SyID_mild = sim.SyID_mild + SyID_mild_sample[1:]
        SyID_mild_carryover = SyID_mild_sample[0]

        new_SyID_severe = pop_updates.get('SyID_severe', 0)
        pop_increase += new_SyID_severe
        SyID_severe_sample = sim.sample_SyID_severe_times(new_SyID_severe + ID_carryover_severe)
        sim.SyID_severe = sim.SyID_severe + SyID_severe_sample[1:]
        SyID_severe_carryover = SyID_severe_sample[0]

        new_R = pop_updates.get('new_R', 0)
        pop_increase += new_R
        sim.R = sim.R + new_R + SyID_mild_carryover + SyID_severe_carryover
        sim.R_mild = sim.R_mild + new_R + SyID_mild_carryover
        sim.R_severe = sim.R_severe + SyID_severe_carryover

        new_QI = pop_updates.get('QI', 0)
        pop_increase += new_QI
        sim.QI = sim.QI + new_QI
        sim.QI_mild = sim.QI_mild + new_QI

        new_QS = pop_updates.get('QS', 0)
        pop_increase += new_QS
        sim.QS = sim.QS + new_QS

        sim.pop_size = sim.pop_size + pop_increase

            
