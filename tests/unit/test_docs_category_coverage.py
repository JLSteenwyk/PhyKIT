import re
from pathlib import Path

import pytest


DOCS_PATH = Path(__file__).resolve().parents[2] / "docs" / "usage" / "index.rst"


def _parse_docs():
    """Return (labels, category_refs) from the usage docs."""
    content = DOCS_PATH.read_text()

    # Every function heading has a ``.. _cmd-<name>:`` label above it
    labels = set(re.findall(r"\.\. _cmd-(\w+):", content))

    # The analytical-category section sits between these two headings
    cat_start = content.index("Functions by analytical category")
    cat_end = content.index(
        "Alignment-based functions\n-------------------------"
    )
    category_section = content[cat_start:cat_end]
    refs = set(re.findall(r":ref:.*?<cmd-(\w+)>", category_section))

    return labels, refs


class TestDocsCategoryCoverage:
    """Ensure every documented function appears in the analytical-category section."""

    def test_every_function_has_a_category(self):
        labels, refs = _parse_docs()
        missing = sorted(labels - refs)
        assert not missing, (
            f"Functions missing from 'Functions by analytical category': {missing}"
        )

    def test_no_stale_category_refs(self):
        labels, refs = _parse_docs()
        stale = sorted(refs - labels)
        assert not stale, (
            f"Category refs pointing to non-existent function labels: {stale}"
        )

    def test_at_least_one_function_documented(self):
        labels, refs = _parse_docs()
        assert len(labels) > 0, "No function labels found in docs"
        assert len(refs) > 0, "No category refs found in docs"
