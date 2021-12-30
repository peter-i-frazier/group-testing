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
infections_per_day_per_contact_unit = np.array([0.28, 0.1, 0.05])


def main():

    params = NOMINAL
    # Include the parameters defined in the JSON file
    json_params = json.load(open(params["json_path"]))
    params.update(json_params)

    GENRATION_TIME = params["generation_time"]

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

    # [11/27 to 12/9] 1x / week testing with 36hr delay
    infections_per_contact = infections_per_day_per_contact_unit * micro.days_infectious(7,1.5)
    infection_rate = popul.infection_matrix(infections_per_contact)
    s = sim(T, S0, I0, R0, infection_rate=infection_rate, generation_time=GENRATION_TIME)
    s.step(3)

    # [12/10 to 12/16] 1x / week testing with 3 day delay
    infections_per_contact = infections_per_day_per_contact_unit * micro.days_infectious(7,3)
    infection_rate = popul.infection_matrix(0.5 * infections_per_contact)
    s.step(1, infection_rate=infection_rate)

    infections_per_contact = infections_per_day_per_contact_unit * micro.days_infectious(7,3)
    infection_rate = popul.infection_matrix(0.33 * infections_per_contact)
    s.step(1, infection_rate=infection_rate)

    # ==================================================
    # [Plot] Actual Infections vs. Simulation Infections
    # ==================================================

    colors = ["navy", "royalblue", "powderblue"]
    X = np.arange(s.max_T)*s.generation_time

    # plot simulated
    groups = popul.metagroup_indices(params["population_names"][:3])
    for i in range(len(groups)):
        group_name = params["population_names"][i]
        infectious = s.get_total_infected_for_different_groups(groups[i], cumulative=True)
        plt.plot(X, infectious, 'k--', label=f"simulated {group_name}", color=colors[i])

    # plot actual counts
    for i in range(3):
        group_name = params["population_names"][i]
        dec_daily_positives = list(pd.read_csv(f"data/dec_infections_{group_name}.csv")['positives'])
        dec_positives = np.cumsum(dec_daily_positives)
        plt.plot(np.arange(20), dec_positives[:20], label=f"actual {group_name}", color=colors[i])

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
