#!/usr/bin/env python3
"""One-time migration from the monolithic tutorial page to stable pages."""

from __future__ import annotations

import re
from pathlib import Path

from phykit.cli_registry import COMMAND_IDENTITIES


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
INDEX = DOCS / "tutorials" / "index.rst"
PAGES = DOCS / "tutorials" / "pages"


PUBLIC_TO_CANONICAL = {
    name: identity.canonical
    for identity in COMMAND_IDENTITIES
    for name in (identity.canonical, *identity.aliases)
}


def slugify(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")


def tutorial_sections(text: str) -> list[tuple[int, str, str]]:
    lines = text.splitlines()
    starts = []
    for index in range(len(lines) - 1):
        match = re.match(r"^(\d+)\.\s+(.+)$", lines[index].strip())
        if match and lines[index + 1].strip() and set(lines[index + 1].strip()) == {"#"}:
            starts.append((index, int(match.group(1)), match.group(2).strip()))
    sections = []
    for position, (start, number, title) in enumerate(starts):
        end = starts[position + 1][0] if position + 1 < len(starts) else len(lines)
        sections.append((number, title, "\n".join(lines[start + 2:end]).strip()))
    return sections


def related_commands(body: str) -> list[str]:
    found = []
    for name in re.findall(r"\bphykit\s+([a-z][a-z0-9_]*)", body):
        canonical = PUBLIC_TO_CANONICAL.get(name)
        if canonical and canonical not in found:
            found.append(canonical)
    return found


def portable_body(body: str) -> str:
    body = body.replace("tests/sample_files/", "")
    body = body.replace(".. image:: ../_static/", ".. image:: /_static/")
    body = re.sub(
        r"Both files are included in the PhyKIT test suite at ``[^`]+``\s+and\s+``[^`]+``\.",
        "Download both files above into the tutorial working directory.",
        body,
    )
    return body


def page(number: int, title: str, body: str) -> str:
    slug = slugify(title)
    body = portable_body(body)
    commands = related_commands(body)
    command_links = "\n".join(
        f"- :doc:`{command.replace('_', ' ').title()} </reference/commands/{command}>`"
        for command in commands
    ) or "- No PhyKIT command is required for this conceptual tutorial."
    full_title = f"Tutorial {number}: {title}"
    return f""".. _tutorial-{number:02d}:
.. _tutorial-{slug}:

{full_title}
{'=' * len(full_title)}

Objectives
----------

- Complete the {title.lower()} workflow.
- Interpret the reported values and generated artifacts in their scientific context.
- Identify the canonical command references for each analysis step.

Prerequisites and working directory
-----------------------------------

Install the current PhyKIT release and create a dedicated working directory.
Download the data linked in this tutorial into that directory before running
the commands. All paths below are relative to this directory.

.. code-block:: shell

   mkdir phykit-tutorial-{number:02d}
   cd phykit-tutorial-{number:02d}

Related command references
--------------------------

{command_links}

Workflow
--------

{body}

Expected artifacts
------------------

Each step identifies its expected terminal output or generated files. Confirm
that those artifacts exist before continuing to the next step; filenames are
relative to the tutorial working directory unless an absolute path is shown.

Troubleshooting
---------------

- Run ``phykit <command> --help`` to compare an invocation with the live interface.
- Confirm that downloaded files are in the current working directory and retain
  the filenames shown in the tutorial.
- For parsing errors, compare taxon names exactly across alignments, trees, and
  trait tables, including capitalization and underscores.
- See :doc:`Troubleshooting </troubleshooting/index>` for installation, format,
  and error-reporting guidance.
"""


def main() -> None:
    text = INDEX.read_text()
    sections = tutorial_sections(text)
    if len(sections) != 21:
        raise SystemExit(f"expected 21 tutorials, found {len(sections)}")
    PAGES.mkdir(parents=True, exist_ok=True)
    entries = []
    for number, title, body in sections:
        slug = f"{number:02d}-{slugify(title)}"
        (PAGES / f"{slug}.rst").write_text(page(number, title, body))
        entries.append((number, title, slug, slugify(title)))
    links = "\n".join(
        f"{number}. :doc:`{title} <pages/{slug}>`" for number, title, slug, _ in entries
    )
    toctree = "\n".join(f"   pages/{slug}" for _, _, slug, _ in entries)
    legacy = "\n\n".join(
        f".. raw:: html\n\n   <span id=\"{old_slug}\"></span>\n\n"
        f":doc:`Tutorial {number}: {title} <pages/{slug}>`"
        for number, title, slug, old_slug in entries
    )
    INDEX.write_text(
        f""".. _tutorials:

Tutorials
=========

These task-oriented tutorials use downloadable, versioned sample data and link
every workflow back to the canonical command reference.

Tutorial index
--------------

{links}

.. toctree::
   :hidden:

{toctree}

Legacy tutorial anchors
-----------------------

{legacy}
"""
    )


if __name__ == "__main__":
    main()
