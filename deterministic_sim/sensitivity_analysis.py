import yaml
import json
import numpy as np
from groups import population
import metrics
from transform import transform
from sim_helper import sim_test_strategy
from sp22_strategies import surge_testing_strategy
import plotting

COLORS = ['#084594', '#2171b5', '#4292c6', '#6baed6',
          '#9ecae1', '#c6dbef', '#deebf7', '#f7fbff']

# =======================
# [Initialize Parameters]
# =======================

nominal_scenario = yaml.safe_load(open("nominal.yaml", "r"))
nominal_scenario["meta_matrix"] = \
    np.array([list(row.values()) for row in nominal_scenario["meta_matrix"].values()])

# ==========================
# [Run] Compare trajectories
# ==========================

trajectories = []
scalers = np.linspace(0.25,4,8)
for i in range(len(scalers)):
    scenario = transform(nominal_scenario,
                         {"winter_break_infections_global": scalers[i]})
    traj = sim_test_strategy(scenario=scenario,
                             strategy=surge_testing_strategy(scenario),
                             color=COLORS[i],
                             name=str(scalers[i]))
    trajectories.append(traj)

# ======
# [Plot]
# ======

popul = population.from_scenario(nominal_scenario)
plotting.plot_metric_over_time(outfile="sensitivity_analysis.png",
                                trajectories=trajectories,
                                metric_name="Cumulative Hospitalizatins",
                                metric=metrics.get_cumulative_hospitalizations)
