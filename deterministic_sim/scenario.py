import numpy as np
import pandas as pd
from typing import Dict, Union


# https://stackoverflow.com/a/54829922
def set_for_keys(my_dict, key_arr, value):
    """Set value at path in my_dict defined by key array."""
    current = my_dict
    for i in range(len(key_arr)):
        key = key_arr[i]
        if key not in current:
            if i == len(key_arr)-1:
                current[key] = value
            else:
                current[key] = {}
        else:
            if type(current[key]) is not dict:
                print("Given dictionary is not compatible with key structure requested")
                raise ValueError("Dictionary key already occupied")
        current = current[key]
    return my_dict


class ScenarioFamily:

    def __init__(self, nominal: Dict, prior: Dict):
        """Initialize a scenario family.

        Args:
            nominal (Dict): Parameters of the nominal scenario (x0).
            prior (Dict): Independent truncated normal marginals (w) used to \
                scale nominal parameters. For each scaler, the dictionary
                should indicate which parameters it scales (B).

        """
        self.flattened = pd.json_normalize(nominal, sep='/').iloc[0].to_dict()
        scaled_params = [i for _, v in prior.items() for i in v["scales"]]

        # initialize dictionary from sim parameters to indices
        self.raw_param_to_index = \
            {scaled_params[i] : i for i in range(len(scaled_params))}
        self.x0 = np.array([self.flattened[i] for i in scaled_params])

        # initialize dictionary from scale parameters to indices
        scale_params = list(prior.keys())
        self.scale_param_to_index = \
            {scale_params[i] : i for i in range(len(scale_params))}

        # initialize B matrix
        self.B = np.zeros((len(scaled_params), len(prior)))
        params = list(prior.items())
        for j in range(len(params)):
            for scales in params[j][1]["scales"]:
                i = self.raw_param_to_index[scales]
                self.B[i,j] = 1

    def _get_scenario_raw(self, x):
        """Return parameters dictionary consumable by the simulator."""
        raw_params = self.flattened.copy()
        for k,v in self.raw_param_to_index.items():
            raw_params[k] = x[v]

        result = {}
        for k,v in raw_params.items():
            set_for_keys(result, k.split('/'), v)

        return result

    def get_scenario_from_w(self, w: Union[np.ndarray, Dict]):
        """Return scenario with w sampled from the prior [x0 + (x0 * Bw)].

        Args:
            w (Union[np.ndarray, Dict]): Vector w sampled from the prior either \
                as a NumPy array or dictionary. If a dictionary, not all scale \
                parameters need be passed. Default to 0 if not passed."""
        if type(w) == dict:
            w_tmp = np.zeros(len(self.scale_param_to_index))
            for k,v in w.items():
                w_tmp[self.scale_param_to_index[k]] = v
            w = w_tmp
        x = self.x0 + (self.x0 * np.matmul(self.B, w))
        return self._get_scenario_raw(x)

    def get_scenario(self, direction: Union[np.ndarray, Dict], distance: float):
        """Return a scenario some [distance] in [direction] from the nominal.

        Args:
            direction (Union[np.ndarray, Dict]): Direction from the nominal. \
                The vector is normalized automatically.
            distance (float): Distance from the nominal scenario.
        """
        if type(direction) == dict:
            direction_tmp = np.zeros(len(self.scale_param_to_index))
            for k,v in direction.items():
                direction_tmp[self.scale_param_to_index[k]] = v
            direction = direction_tmp

        n = direction / np.linalg.norm(direction)
        w = distance * n
        return self.get_scenario_from_w(w)

# TODO (hwr26): Implement this
class PessimismLevelScenarioFamily:

    def __init__():
        pass
