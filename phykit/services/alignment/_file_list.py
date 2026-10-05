"""Shared path spelling rules for alignment and tree file lists."""

import os


def _normalize_list_path(path: str) -> str:
    separator = os.sep
    if not (
        separator + separator in path
        or path == "."
        or path.startswith("." + separator)
        or separator + "." + separator in path
        or path.endswith(separator + ".")
    ):
        return path

    is_absolute = os.path.isabs(path)
    if is_absolute:
        if (
            separator == "/"
            and path.startswith("//")
            and not path.startswith("///")
        ):
            prefix = "//"
            rest = path[2:]
        else:
            prefix = separator
            rest = path.lstrip(separator)
    else:
        prefix = ""
        rest = path

    parts = [
        part for part in rest.split(separator)
        if part and part != "."
    ]
    normalized = separator.join(parts)
    if prefix:
        return prefix + normalized if normalized else prefix
    return normalized or "."
