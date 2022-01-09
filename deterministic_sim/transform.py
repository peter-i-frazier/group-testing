import yaml
import pandas as pd
from typing import Dict


TRANSFORMATIONS = yaml.safe_load(open("transformations.yaml", "r"))


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


def transform(scenario: Dict, transformations: Dict) -> Dict:
    """Transform the given scenario according to the given transformations.

    Args:
        scenario (Dict): Scenario to transform.
        transformations (Dict): Transformations to apply.

    Returns:
        Dict: Scenario with transformaions applied.
    """
    flattened_scenario = pd.json_normalize(scenario, sep='/').iloc[0].to_dict()

    # TODO (hwr26): move out to only be computed once
    # define recognized transformations
    recognized = TRANSFORMATIONS
    for k in flattened_scenario.keys():
        # add operation for individual parameter
        recognized[f"{k}_add"] = {}
        recognized[f"{k}_add"]["op"] = "add"
        recognized[f"{k}_add"]["affected"] = {k: 1}
        # scale operation for individual parameter
        recognized[f"{k}_scale"] = {}
        recognized[f"{k}_scale"]["op"] = "scale"
        recognized[f"{k}_scale"]["affected"] = {k: 1}

    # get list of transformations for each transformed parameter
    transformed_parameters = {}
    for k, value in transformations.items():
        transformation = recognized[k]
        for affected_param, weight in transformation["affected"].items():
            tmp = [transformation["op"], weight * value]
            if affected_param in transformed_parameters:
                transformed_parameters[affected_param].append(tmp)
            else:
                transformed_parameters[affected_param] = [tmp]

    # apply transformations
    for k,v in transformed_parameters.items():
        first_op = v[0][0]
        for op, value in v:
            if op == first_op:
                if op == "add":
                    flattened_scenario[k] += value
                elif op == "scale":
                    flattened_scenario[k] *= value
                else:
                    raise ValueError(f"{op} is an unsupported operation.")
            else:
                raise ValueError(f"Unsupported: multiple ops provided for {k}")

    # un-flatten
    transformed_scenario = {}
    for k,v in flattened_scenario.items():
        set_for_keys(transformed_scenario, k.split('/'), v)

    return transformed_scenario
