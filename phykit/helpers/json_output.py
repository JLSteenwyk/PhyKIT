import json

import numpy as np


def to_builtin_json_types(value):
    if isinstance(value, dict):
        return {key: to_builtin_json_types(sub_value) for key, sub_value in value.items()}
    if isinstance(value, list):
        return [to_builtin_json_types(sub_value) for sub_value in value]
    if isinstance(value, tuple):
        return [to_builtin_json_types(sub_value) for sub_value in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return [to_builtin_json_types(sub_value) for sub_value in value.tolist()]
    return value


def print_json(payload, sort_keys=True):
    try:
        print(json.dumps(to_builtin_json_types(payload), sort_keys=sort_keys))
    except BrokenPipeError:
        pass
