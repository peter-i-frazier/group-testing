import yaml
import json
import numpy as np
import pandas as pd
import micro
from sim import sim
from groups import meta_group
import matplotlib.pyplot as plt

NOMINAL = yaml.safe_load(open("nominal.yaml", "r"))

T = 5                           # 4 generations (16 days)
INITIAL_INFECTIOUS = [6, 1, 1]  # initial infections for each group

# This is the parameter we aim to calibrate. It is the number of infections per
# day per (metagroup-specific contact unit). This script calibrates this
# parameter to the December Omicron surge.
infections_per_day_per_contact_unit = [0.22, 0.15, 0.13]


def main():

    params = NOMINAL
    # Include the parameters defined in the JSON file
    json_params = json.load(open(params["json_path"]))
    params.update(json_params)

    GENRATION_TIME = params["generation_time"]
    SYMPTOMATIC_RATE = params["symptomatic_rate"]

    # only look at student groups
    for i in range(3):

        # ===================================================================
        # [Initialize] Assume no recovered and set initial Omicron infections
        # ===================================================================

        pop = params["population_count"][i] * np.array(params["pop_fracs"][i])
        marginal_contacts = np.arange(1, len(params["pop_fracs"][i]) + 1)
        K = len(marginal_contacts)
        R0 = np.zeros(K)
        I0 = np.zeros(K)
        # Assume the most social of the population of gets infected first.
        I0[K - 1] = INITIAL_INFECTIOUS[i]
        S0 = np.maximum(pop - R0 - I0, 0)

        # ========================================================================
        # [Run] Increase testing delay and reduce interaction over duration of sim
        # ========================================================================

        # [12/1 to 12/9] 1x / week testing with 36hr delay
        infections_per_contact = infections_per_day_per_contact_unit[i] * micro.days_infectious(7,1.5)
        infection_rate = meta_group("UG", pop, marginal_contacts).infection_matrix(infections_per_contact)
        s = sim(T, S0, I0, R0, infection_rate=infection_rate, generation_time=GENRATION_TIME)
        s.step(2)

        # [12/10 to 12/16] 1x / week testing with 3 day delay
        infections_per_contact = infections_per_day_per_contact_unit[i] * micro.days_infectious(7,3)
        infection_rate = meta_group("UG", pop, 0.5 * marginal_contacts).infection_matrix(infections_per_contact)
        s.step(1, infection_rate=infection_rate)

        infections_per_contact = infections_per_day_per_contact_unit[i] * micro.days_infectious(7,3)
        infection_rate = meta_group("UG", pop, 0.33 * marginal_contacts).infection_matrix(infections_per_contact)
        s.step(1, infection_rate=infection_rate)

        # ==================================================
        # [Plot] Actual Infections vs. Simulation Infections
        # ==================================================

        group_name = params["population_names"][i]
        I = s.get_metric('I', aggregate=True, cumulative=True)
        plt.plot(np.arange(T)*GENRATION_TIME, I, label="simulated")

        # plot actual counts
        dec_daily_positives = list(pd.read_csv(f"data/dec_infections_{group_name}.csv")['positives'])
        dec_positives = np.cumsum(dec_daily_positives)
        plt.plot(np.arange(15), dec_positives[:15], label="actual")

        plt.title(f"Actual vs. Simulated Infection Trajectories [{group_name}]\n"
                f"infections_per_day_per_contact_unit={infections_per_day_per_contact_unit[i]}, Symptomatic Rate = {SYMPTOMATIC_RATE}")
        plt.rcParams.update({'font.size': 8})
        plt.xlabel('Days Since Dec. 1')
        plt.ylabel('Cumulative Infected')
        plt.legend()
        plt.savefig(f'calibration-student-{group_name}.png', facecolor='w')
        plt.close()


if __name__ == "__main__":
    main()
