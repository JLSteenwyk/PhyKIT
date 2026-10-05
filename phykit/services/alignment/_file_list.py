"""Shared path spelling rules for alignment and tree file lists."""

import os

from ...errors import PhykitUserError


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


def read_file_list(path, *, path_factory, normalize_path):
    """Read paths relative to the list, retaining service-specific lazy hooks."""
    source = path_factory(path)
    if not source.exists():
        raise PhykitUserError(
            [f"{path} corresponds to no such file or directory."],
            code=2,
        )

    paths = []
    append = paths.append
    parent_str = str(source.parent)
    parent_prefix = "" if parent_str == "." else parent_str + os.sep
    separator = os.sep
    double_separator = separator + separator
    dot_prefix = "." + separator
    dot_segment = separator + "." + separator
    dot_suffix = separator + "."
    check_isabs = os.path.isabs if separator == "\\" else None
    with source.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line[0] == "#":
                continue
            if line == ".":
                append(parent_str)
                continue
            has_separator = separator in line
            if not has_separator:
                append(parent_prefix + line)
                continue
            needs_normalization = (
                line.startswith(dot_prefix)
                or double_separator in line
                or dot_segment in line
                or line.endswith(dot_suffix)
            )
            if line[0] == separator or (
                check_isabs is not None and check_isabs(line)
            ):
                append(
                    normalize_path(line)
                    if needs_normalization
                    else line
                )
            elif needs_normalization:
                normalized = normalize_path(line)
                if normalized == ".":
                    append(parent_str)
                else:
                    append(parent_prefix + normalized)
            else:
                append(parent_prefix + line)
    return paths
