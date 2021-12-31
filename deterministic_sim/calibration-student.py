import yaml
import json
import numpy as np
import pandas as pd
import micro
from sim import sim
from groups import meta_group, population
import matplotlib.pyplot as plt
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

NOMINAL = yaml.safe_load(open("nominal.yaml", "r"))

T = 6                           # 4 generations (20 days)
INITIAL_INFECTIOUS = [1, 0, 0]  # initial infections for each group
PAST_INFECTIONS = [0, 0, 0]     # recovered for each group

# This is the parameter we aim to calibrate. It is the number of infections per
# day per (metagroup-specific contact unit). This script calibrates this
# parameter to the December Omicron surge.
# [UG, GR, PR]
# TODO (hwr26): Investigate this. Seems a little strange right now.
# TODO (pf98): I agree.  I tuned this to minimize MSE and came up with [0.27, 0.13, 0.05] ---
#  surprising that GR has more contact than PR. But the MSE changes only by a small amount
#  as we change the parameters for GR and PR.  Leaving them both at 0.1 for now.
infections_per_day_per_contact_unit = np.array([0.28, 0.1, 0.1])

# UG and PR were in 1x / wk surveillance
# GR were not in surveillance and so had an infinite number of days between scheduled tests
days_between_scheduled_tests = np.array([7, np.inf, 7])


def main():

    params = NOMINAL
    # Include the parameters defined in the JSON file
    json_params = json.load(open(params["json_path"]))
    params.update(json_params)

    GENERATION_TIME = params["generation_time"]

    # ===================================================================
    # [Initialize] Assume no recovered and set initial Omicron infections
    # ===================================================================

    # ONLY INCLUDE THE STUDENT POPULATIONS
    population_count = params["population_count"][:3]
    population_names = params["population_names"][:3]
    initial_infections = INITIAL_INFECTIOUS
    past_infections = PAST_INFECTIONS
    meta_groups = []
    for i in range(len(population_count)):
        name = population_names[i]
        pop = population_count[i] * np.array(params['pop_fracs'][i])
        contact_units = np.arange(1, len(pop) + 1)
        meta_groups.append(meta_group(name, pop, contact_units))

    popul = population(meta_groups, np.array(params['meta_matrix'])[:3,:3])
    # Assume that infections were among the most social in the population
    S0, I0, R0 = popul.get_init_SIR_vec(initial_infections, past_infections,
                                        weight="most_social")

    # ========================================================================
    # [Run] Increase testing delay and reduce interaction over duration of sim
    # ========================================================================

    # [11/27 to 12/9] 36hr testing delay
    days_infectious = [micro.days_infectious(d,1.5) for d in days_between_scheduled_tests ]
    infections_per_contact = infections_per_day_per_contact_unit * days_infectious
    infection_rate = popul.infection_matrix(infections_per_contact)
    s = sim(T, S0, I0, R0, infection_rate=infection_rate, generation_time=GENERATION_TIME)
    s.step(3)

    # [12/10 to 12/16] 3 day testing delay.
    # In the first generation, we model contacts as being 50% of what they were during study
    # week, because Cornell moved to Yellow COVID status and final exams were beginning.
    # In the second generation, we model contacts as being 33% of what they were during study
    # week, because Cornell went to Red COVID status.
    days_infectious = [micro.days_infectious(d,3) for d in days_between_scheduled_tests ]
    infections_per_contact = infections_per_day_per_contact_unit * days_infectious
    infection_rate = popul.infection_matrix(0.5 * infections_per_contact)
    s.step(1, infection_rate=infection_rate)

    infection_rate = popul.infection_matrix(0.33 * infections_per_contact)
    s.step(1, infection_rate=infection_rate)

    # ==================================================
    # [Plot] Actual Infections vs. Simulation Infections
    # ==================================================

    colors = ["navy", "royalblue", "powderblue"]
    X = np.arange(s.max_T)*s.generation_time

    groups = popul.metagroup_indices(params["population_names"][:3])
    total_MSE = 0
    for i in range(len(groups)):
        # plot simulated
        group_name = params["population_names"][i]
        infectious = s.get_total_infected_for_different_groups(groups[i], cumulative=True)
        plt.plot(X, infectious, 'k--', label=f"simulated {group_name}", color=colors[i])

        # plot actual counts
        group_name = params["population_names"][i]
        dec_daily_positives = list(pd.read_csv(f"data/dec_infections_{group_name}.csv")['positives'])
        dec_positives = np.cumsum(dec_daily_positives)
        plt.plot(np.arange(20), dec_positives[:20], label=f"actual {group_name}", color=colors[i])

        # print mean squared errors
        MSE = 0
        for i in range(len(X)):
            # The fact that infectious[i] and dec_positives[X[i]] are referring to the same day
            # is because dec_positives[X[i]] happened on day X[i]
            MSE = MSE + np.power(infectious[i] - dec_positives[X[i]],2)
        MSE = MSE / len(X)
        total_MSE = total_MSE + MSE
        print('{} MSE={}'.format(group_name,MSE))
    total_MSE = total_MSE / 3 # Since we summed the mean MSE from 3 groups with equally many datapoints
    print(infections_per_day_per_contact_unit )
    print('Total MSE={}'.format(total_MSE))

    plt.title(f"Actual vs. Simulated Infection Trajectories [Students]\n" + \
              f"infections_per_day_per_contact_unit: {str(infections_per_day_per_contact_unit)}")
    plt.rcParams.update({'font.size': 8})
    plt.xlabel('Days Since Nov. 27')
    plt.ylabel('Cumulative Infected')
    plt.legend()
    plt.savefig(f'calibration-student.png', facecolor='w')
    plt.close()


if __name__ == "__main__":
    main()
