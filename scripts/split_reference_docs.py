#!/usr/bin/env python3
"""One-time migration from the monolithic usage page to command pages."""

from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path

from phykit.cli_registry import COMMAND_IDENTITIES


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
USAGE = DOCS / "usage" / "index.rst"
REFERENCE = DOCS / "reference"
COMMANDS = REFERENCE / "commands"
CATEGORIES = REFERENCE / "categories"
SPEC_PATH = DOCS / "_data" / "command_spec.json"


def slugify(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")


def parse_categories(text: str) -> tuple[dict[str, list[str]], dict[str, str]]:
    start = text.index("Functions by analytical category")
    end = text.index("Alignment-based functions\n-------------------------")
    lines = text[start:end].splitlines()
    categories: dict[str, list[str]] = defaultdict(list)
    summaries: dict[str, str] = {}
    category = "Other utilities"
    for index, line in enumerate(lines):
        if index + 1 < len(lines) and set(lines[index + 1]) == {"#"}:
            category = line.strip()
        match = re.match(
            r"- :ref:`[^<]+<cmd-([\w]+)>`:\s*(.+)", line.strip()
        )
        if match:
            command, summary = match.groups()
            categories[command].append(category)
            summaries.setdefault(command, summary.strip())
    return dict(categories), summaries


def parse_command_blocks(text: str) -> dict[str, str]:
    matches = list(re.finditer(r"^\.\. _cmd-([\w]+):\n", text, re.MULTILINE))
    blocks = {}
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        block = text[match.start():end]
        block = re.sub(
            r"\n\|\n+(?:[A-Za-z][^\n]+\n-+\n+)?\Z", "\n", block.rstrip() + "\n"
        )
        blocks[match.group(1)] = block.rstrip() + "\n"
    return blocks


def command_title(block: str, fallback: str) -> str:
    lines = block.splitlines()
    for index, line in enumerate(lines[:-1]):
        if line.strip() and set(lines[index + 1].strip()) == {"#"}:
            return line.strip()
    return fallback.replace("_", " ").title()


def command_page(identity, block: str, summary: str, categories: list[str]) -> str:
    title = command_title(block, identity.canonical)
    lines = block.splitlines()
    title_index = next(
        index
        for index, line in enumerate(lines[:-1])
        if line.strip() == title and set(lines[index + 1].strip()) == {"#"}
    )
    body = "\n".join(lines[title_index + 2:]).strip()
    body = re.sub(
        r"^Function names:.*?\nCommand line interface:.*?(?:\n|$)",
        "",
        body,
        count=1,
        flags=re.MULTILINE,
    ).strip()
    body = body.replace(".. image:: ../_static/", ".. image:: /_static/")
    aliases = ", ".join(identity.aliases) or "none"
    entry_points = ", ".join(identity.entry_points)
    category_text = ", ".join(categories) or "Other utilities"
    return f""".. _cmd-{identity.canonical}:
.. _command-{identity.canonical}:

{title}
{'=' * len(title)}

{summary}

Command identity
----------------

:Canonical command: ``{identity.canonical}``
:Handler: ``{identity.handler}``
:Aliases: {aliases}
:Standalone executables: {entry_points}
:Categories: {category_text}

Runtime interface
-----------------

.. include:: /_generated/commands/{identity.canonical}.inc

Guidance, interpretation, and examples
--------------------------------------

{body}
"""


def category_page(title: str, commands: list[dict]) -> str:
    rows = "\n".join(
        f"- :doc:`{item['title']} <../commands/{item['canonical']}>` - {item['summary']}"
        for item in sorted(commands, key=lambda item: item["title"].lower())
    )
    toctree = "\n".join(
        f"   ../commands/{item['canonical']}"
        for item in sorted(commands, key=lambda item: item["canonical"])
    )
    return f"""{title}
{'=' * len(title)}

{rows}

.. toctree::
   :hidden:

{toctree}
"""


def main() -> None:
    text = USAGE.read_text()
    if "Alignment-based functions\n-------------------------" not in text:
        raise SystemExit("usage page has already been migrated")

    category_map, summaries = parse_categories(text)
    blocks = parse_command_blocks(text)
    identities = {identity.canonical: identity for identity in COMMAND_IDENTITIES}
    unknown = sorted(set(blocks) - set(identities))
    if unknown:
        raise SystemExit(f"documented commands missing from registry: {unknown}")

    COMMANDS.mkdir(parents=True, exist_ok=True)
    CATEGORIES.mkdir(parents=True, exist_ok=True)
    SPEC_PATH.parent.mkdir(parents=True, exist_ok=True)

    records = []
    by_category: dict[str, list[dict]] = defaultdict(list)
    for canonical, identity in sorted(identities.items()):
        block = blocks.get(canonical)
        title = command_title(block, canonical) if block else canonical.replace("_", " ").title()
        summary = summaries.get(canonical, f"Run the {title} analysis or utility.")
        categories = category_map.get(canonical, ["Other utilities"])
        if block is None and canonical == "version":
            block = """.. _cmd-version:\n\nVersion\n#######\n\nPrint the installed PhyKIT version.\n\n.. code-block:: shell\n\n   phykit version\n"""
        elif block is None:
            raise SystemExit(f"registry command has no documentation: {canonical}")
        (COMMANDS / f"{canonical}.rst").write_text(
            command_page(identity, block, summary, categories)
        )
        record = {
            "canonical": canonical,
            "handler": identity.handler,
            "aliases": list(identity.aliases),
            "entry_points": list(identity.entry_points),
            "title": title,
            "summary": summary,
            "categories": categories,
            "page": f"reference/commands/{canonical}.rst",
            "documentation": {
                "purpose": summary,
                "use_cases": [summary],
                "input_formats": [],
                "outputs": [],
                "errors": [],
                "limitations": [],
                "scientific_interpretation": "See the canonical command page.",
                "validation": [],
                "citations": [],
                "examples": [],
                "related_commands": [],
                "related_tutorials": [],
            },
        }
        records.append(record)
        for category in categories:
            by_category[category].append(record)

    SPEC_PATH.write_text(
        json.dumps({"schema_version": 1, "commands": records}, indent=2) + "\n"
    )

    category_links = []
    for category, commands in sorted(by_category.items()):
        slug = slugify(category)
        (CATEGORIES / f"{slug}.rst").write_text(category_page(category, commands))
        category_links.append((category, slug, len(commands)))

    links = "\n".join(
        f"- :doc:`{title} <categories/{slug}>` ({count} commands)"
        for title, slug, count in category_links
    )
    category_toctree = "\n".join(f"   categories/{slug}" for _, slug, _ in category_links)
    REFERENCE.joinpath("index.rst").write_text(
        f"""Command reference
=================

PhyKIT commands are organized below by analytical task. Each command has one
canonical page containing its identity, live runtime interface, scientific
guidance, output interpretation, examples, limitations, validation, and citations.

Choose by task
--------------

{links}

All commands
------------

Use :download:`the machine-readable command catalog </_static/commands.json>`
for programmatic discovery of commands, aliases, arguments, and documentation URLs.

.. toctree::
   :hidden:

{category_toctree}
"""
    )

    general_end = text.index("Functions by analytical category")
    general = text[:general_end].rstrip()
    compatibility = []
    for record in records:
        canonical = record["canonical"]
        compatibility.append(
            f".. raw:: html\n\n   <span id=\"cmd-{canonical}\"></span>\n\n"
            f"- :doc:`{record['title']} </reference/commands/{canonical}>`"
        )
    USAGE.write_text(
        general
        + "\n\nCommand reference\n-----------------\n\n"
        + "The complete reference is organized by task on the "
          ":doc:`command reference page </reference/index>`.\n\n"
        + "Legacy command anchors\n######################\n\n"
        + "These links preserve the command fragments used by earlier versions of the documentation.\n\n"
        + "\n\n".join(compatibility)
        + "\n"
    )


if __name__ == "__main__":
    main()
