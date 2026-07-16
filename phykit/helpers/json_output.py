from math import isfinite


_JSON_SCALAR_TYPES = (str, int, float, bool)
_LARGE_DICT_COPY_THRESHOLD = 50_000
_JSON_DUMPS = None


def to_builtin_json_types(value):
    value_type = type(value)
    if value is None:
        return value
    if value_type is float:
        return value if isfinite(value) else None
    if value_type in _JSON_SCALAR_TYPES:
        return value
    if isinstance(value, dict):
        converted = None
        use_copy_on_convert = (
            type(value) is dict and len(value) >= _LARGE_DICT_COPY_THRESHOLD
        )
        for key, sub_value in value.items():
            sub_value_type = type(sub_value)
            if sub_value_type is float:
                converted_value = (
                    sub_value if isfinite(sub_value) else None
                )
            elif sub_value is None or sub_value_type in _JSON_SCALAR_TYPES:
                converted_value = sub_value
            else:
                converted_value = to_builtin_json_types(sub_value)
            if converted is None:
                if converted_value is sub_value:
                    continue
                if use_copy_on_convert:
                    converted = value.copy()
                else:
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
            sub_value_type = type(sub_value)
            if sub_value_type is float:
                converted_value = (
                    sub_value if isfinite(sub_value) else None
                )
            elif sub_value is None or sub_value_type in _JSON_SCALAR_TYPES:
                converted_value = sub_value
            else:
                converted_value = to_builtin_json_types(sub_value)
            if converted is None:
                if converted_value is sub_value:
                    continue
                converted = value[:index]
            converted.append(converted_value)
        return value if converted is None else converted
    if isinstance(value, tuple):
        converted = None
        for index, sub_value in enumerate(value):
            sub_value_type = type(sub_value)
            if sub_value_type is float:
                converted_value = (
                    sub_value if isfinite(sub_value) else None
                )
            elif sub_value is None or sub_value_type in _JSON_SCALAR_TYPES:
                converted_value = sub_value
            else:
                converted_value = to_builtin_json_types(sub_value)
            if converted is None:
                if converted_value is sub_value:
                    continue
                converted = list(value[:index])
            converted.append(converted_value)
        return list(value) if converted is None else converted

    value_type_module = type(value).__module__
    if value_type_module == "numpy" or value_type_module.startswith("numpy."):
        if hasattr(value, "tolist"):
            return to_builtin_json_types(value.tolist())
        if hasattr(value, "item"):
            return to_builtin_json_types(value.item())
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
    global _JSON_DUMPS
    dumps = _JSON_DUMPS
    if dumps is None:
        import json

        dumps = json.dumps
        _JSON_DUMPS = dumps

    try:
        try:
            output = dumps(
                payload,
                sort_keys=sort_keys,
                default=_json_default,
                allow_nan=False,
            )
        except ValueError as error:
            if "Out of range float values" not in str(error):
                raise
            output = dumps(
                to_builtin_json_types(payload),
                sort_keys=sort_keys,
                default=_json_default,
                allow_nan=False,
            )
        print(output)
    except BrokenPipeError:
        pass
