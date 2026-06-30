def to_builtin_json_types(value):
    if value is None or type(value) in (str, int, float, bool):
        return value
    if isinstance(value, dict):
        converted = None
        for key, sub_value in value.items():
            converted_value = to_builtin_json_types(sub_value)
            if converted is None:
                if converted_value is sub_value:
                    continue
                converted = {}
                for previous_key, previous_value in value.items():
                    if previous_key == key:
                        break
                    converted[previous_key] = previous_value
            converted[key] = converted_value
        return value if converted is None else converted
    if isinstance(value, list):
        converted = None
        for index, sub_value in enumerate(value):
            converted_value = to_builtin_json_types(sub_value)
            if converted is None:
                if converted_value is sub_value:
                    continue
                converted = value[:index]
            converted.append(converted_value)
        return value if converted is None else converted
    if isinstance(value, tuple):
        return [to_builtin_json_types(sub_value) for sub_value in value]

    value_type_module = type(value).__module__
    if value_type_module == "numpy" or value_type_module.startswith("numpy."):
        if hasattr(value, "tolist"):
            return to_builtin_json_types(value.tolist())
        if hasattr(value, "item"):
            return value.item()
    if hasattr(value, "tolist") and value.__class__.__name__ == "ndarray":
        return [to_builtin_json_types(sub_value) for sub_value in value.tolist()]
    return value


def _json_default(value):
    value_type_module = type(value).__module__
    if value_type_module == "numpy" or value_type_module.startswith("numpy."):
        if hasattr(value, "tolist"):
            return value.tolist()
        if hasattr(value, "item"):
            return value.item()
    if hasattr(value, "tolist") and value.__class__.__name__ == "ndarray":
        return value.tolist()
    raise TypeError(
        f"Object of type {value.__class__.__name__} is not JSON serializable"
    )


def print_json(payload, sort_keys=True):
    import json

    try:
        output = json.dumps(
            payload,
            sort_keys=sort_keys,
            default=_json_default,
        )
        print(output)
    except BrokenPipeError:
        pass
