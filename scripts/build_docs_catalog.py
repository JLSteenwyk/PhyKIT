#!/usr/bin/env python3
"""Build runtime command tables and the machine-readable command catalog."""

from __future__ import annotations

import argparse
import contextlib
import io
import json
from pathlib import Path
from typing import Any

import phykit.phykit as cli
from phykit.cli_registry import COMMAND_IDENTITIES
from phykit.version import __version__


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
SOURCE_SPEC = DOCS / "_data" / "command_spec.json"
GENERATED = DOCS / "_generated" / "commands"
CATALOG = DOCS / "_static" / "commands.json"


def json_value(value: Any) -> Any:
    if value is argparse.SUPPRESS:
        return None
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, (list, tuple)):
        return [json_value(item) for item in value]
    return str(value)


def type_name(action: argparse.Action) -> str:
    if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
        return "boolean"
    if action.type is None:
        return "string"
    return getattr(action.type, "__name__", str(action.type))


def action_record(action: argparse.Action) -> dict[str, Any]:
    positional = not action.option_strings
    required = action.required or (positional and action.nargs not in ("?", "*"))
    return {
        "name": action.dest,
        "flags": list(action.option_strings),
        "positional": positional,
        "required": required,
        "type": type_name(action),
        "nargs": json_value(action.nargs),
        "default": json_value(action.default),
        "choices": json_value(list(action.choices)) if action.choices is not None else None,
        "metavar": json_value(action.metavar),
        "action": action.__class__.__name__,
    }


def capture_parser(identity) -> argparse.ArgumentParser | None:
    if identity.handler == "version":
        return None
    captured = []
    original = cli._new_parser

    def recording_parser(*, description):
        parser = original(description=description)
        captured.append(parser)
        return parser

    cli._new_parser = recording_parser
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            try:
                getattr(cli.Phykit, identity.handler)(["--help"])
            except SystemExit:
                pass
    finally:
        cli._new_parser = original
    if not captured:
        raise RuntimeError(f"unable to capture parser for {identity.canonical}")
    return captured[-1]


def argument_token(argument: dict[str, Any]) -> str:
    if argument["positional"]:
        token = f"<{argument['name']}>"
    else:
        long_flags = [flag for flag in argument["flags"] if flag.startswith("--")]
        token = long_flags[0] if long_flags else argument["flags"][0]
        if argument["action"] not in ("_StoreTrueAction", "_StoreFalseAction"):
            token += f" <{argument['name']}>"
    return token if argument["required"] else f"[{token}]"


def rst_value(value: Any) -> str:
    if value is None:
        return "none"
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, list):
        return ", ".join(str(item) for item in value)
    return str(value).replace("`", "\\`")


def runtime_include(command: dict[str, Any]) -> str:
    arguments = command["arguments"]
    synopsis = " ".join(
        ["phykit", command["canonical"]]
        + [argument_token(argument) for argument in arguments]
    )
    rows = []
    for argument in arguments:
        names = argument["name"] if argument["positional"] else ", ".join(argument["flags"])
        default = "required" if argument["required"] else rst_value(argument["default"])
        choices = rst_value(argument["choices"]) if argument["choices"] else "any"
        rows.append(
            "\n".join(
                [
                    f"   * - ``{names}``",
                    f"     - {str(argument['required']).lower()}",
                    f"     - {argument['type']}",
                    f"     - {default}",
                    f"     - {choices}",
                ]
            )
        )
    if not rows:
        rows.append("   * - none\n     - false\n     - none\n     - none\n     - none")
    json_note = (
        "``--json`` provides the command's structured JSON representation."
        if any("--json" in argument["flags"] for argument in arguments)
        else "This command does not expose a ``--json`` option."
    )
    output_options = [
        argument for argument in arguments
        if any("output" in flag for flag in argument["flags"])
    ]
    output_note = (
        "Output-file options: "
        + ", ".join(
            f"``{next(flag for flag in argument['flags'] if flag.startswith('--'))}``"
            for argument in output_options
        )
        + "."
        if output_options
        else "Unless the guidance below states otherwise, results are emitted as command output."
    )
    return f"""Synopsis
^^^^^^^^

.. code-block:: text

   {synopsis}

Arguments
^^^^^^^^^

This table is generated from the live command parser. It is the authoritative
source for accepted spellings, required arguments, types, defaults, and choices.

.. list-table::
   :header-rows: 1
   :widths: 28 10 12 18 20

   * - Argument
     - Required
     - Type
     - Default
     - Choices
{chr(10).join(rows)}

Output and errors
^^^^^^^^^^^^^^^^^

{json_note} {output_note} Invalid command syntax exits with status 2. Input
validation and scientific limitations are described in the guidance below.
"""


def build() -> tuple[dict[str, str], str]:
    source = json.loads(SOURCE_SPEC.read_text())
    records = {item["canonical"]: item for item in source["commands"]}
    identities = {item.canonical: item for item in COMMAND_IDENTITIES}
    if set(records) != set(identities):
        raise RuntimeError("command_spec.json and cli_registry.py command sets differ")

    includes = {}
    catalog_commands = []
    for canonical in sorted(records):
        identity = identities[canonical]
        parser = capture_parser(identity)
        arguments = [] if parser is None else [
            action_record(action) for action in parser._actions if action.dest != "help"
        ]
        command = dict(records[canonical])
        command["arguments"] = arguments
        command["documentation_url"] = f"reference/commands/{canonical}.html"
        includes[canonical] = runtime_include(command)
        catalog_commands.append(command)
    catalog = {
        "schema_version": 1,
        "phykit_version": __version__,
        "generated_from": ["phykit/cli_registry.py", "live argparse parsers", "docs/_data/command_spec.json"],
        "commands": catalog_commands,
    }
    return includes, json.dumps(catalog, indent=2) + "\n"


def write_or_check(check: bool) -> int:
    includes, catalog = build()
    mismatches = []
    for canonical, content in includes.items():
        path = GENERATED / f"{canonical}.inc"
        if check:
            if not path.exists() or path.read_text() != content:
                mismatches.append(str(path.relative_to(ROOT)))
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content)
    if check:
        if not CATALOG.exists() or CATALOG.read_text() != catalog:
            mismatches.append(str(CATALOG.relative_to(ROOT)))
        if mismatches:
            print("Generated documentation is stale:")
            print("\n".join(f"- {path}" for path in mismatches))
            return 1
    else:
        CATALOG.parent.mkdir(parents=True, exist_ok=True)
        CATALOG.write_text(catalog)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    return write_or_check(args.check)


if __name__ == "__main__":
    raise SystemExit(main())
