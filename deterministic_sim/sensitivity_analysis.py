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
scalers = np.linspace(0.0,1.0,20)
for i in range(len(scalers)):
    scenario = transform(nominal_scenario,
                         {"booster_rate/FS_multiply_linear_scale": scalers[i]})
    traj = sim_test_strategy(scenario=scenario,
                             strategy=surge_testing_strategy(scenario),
                             color="white",
                             name=str(scenario["booster_rate"]["FS"]))
    trajectories.append(traj)

# ======
# [Plot]
# ======

popul = population.from_scenario(nominal_scenario)
plotting.plot_parameter_sensitivity(outfile="sensitivity_analysis.png",
                                trajectories=trajectories,
                                param_name="Employee Booster Rate",
                                metric_name="Total Employee Hospitalizations",
                                metric=metrics.get_total_employee_hospitalizations)
