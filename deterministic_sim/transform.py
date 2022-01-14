import yaml
import numpy as np
import pandas as pd
from enum import Enum
from typing import Dict

class Op(Enum):
    ADD = "add"
    MULTIPLY = "multiply"


class Scale(Enum):
    LINEAR = "linear"
    LOG = "log"


def flatten_dict(my_dict: Dict) -> Dict:
    """Flatten a dict with the / separator."""
    return pd.json_normalize(my_dict, sep='/').iloc[0].to_dict()


# define recognized transformations
TRANSFORMATIONS = yaml.safe_load(open("transformations.yaml", "r"))
parameters = flatten_dict(yaml.safe_load(open("nominal.yaml", "r"))).keys()
for parameter in parameters:
    for op in Op:
        for scale in Scale:
            name = f"{parameter}_{op.value}_{scale.value}_scale"
            TRANSFORMATIONS[name] = {}
            TRANSFORMATIONS[name]["op"] = op.value
            TRANSFORMATIONS[name]["scale"] = scale.value
            TRANSFORMATIONS[name]["affected"] = {parameter: 1}


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
    flattened_scenario = flatten_dict(scenario)

    # get list of transformations for each transformed parameter
    transformed_parameters = {}
    for k, value in transformations.items():
        transformation = TRANSFORMATIONS[k]
        op = transformation["op"]
        scale = transformation["scale"]
        for affected_param, weight in transformation["affected"].items():
            tmp = [op, scale, weight * value]
            if affected_param in transformed_parameters:
                transformed_parameters[affected_param].append(tmp)
            else:
                transformed_parameters[affected_param] = [tmp]

    # apply transformations
    for k,v in transformed_parameters.items():
        # TODO (hwr26): Rule for limiting to operations on one scale?
        first_op = v[0][0]
        for op, scale, value in v:
            if op == first_op:
                if scale == Scale.LOG:
                    flattened_scenario[k] = np.log(flattened_scenario[k])
                if op == Op.ADD.value:
                    flattened_scenario[k] += value
                elif op == Op.MULTIPLY.value:
                    flattened_scenario[k] *= value
                else:
                    raise ValueError(f"{op} is an unsupported operation.")
                if scale == Scale.LOG:
                    flattened_scenario[k] = np.exp(flattened_scenario[k])
            else:
                raise ValueError(f"Unsupported: multiple ops provided for {k}")

    # un-flatten
    transformed_scenario = {}
    for k,v in flattened_scenario.items():
        set_for_keys(transformed_scenario, k.split('/'), v)

    return transformed_scenario
