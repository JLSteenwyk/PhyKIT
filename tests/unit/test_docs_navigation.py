import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"


def test_primary_navigation_pages_exist_and_are_linked():
    index = (DOCS / "index.rst").read_text()
    required = [
        "getting_started/index",
        "reference/index",
        "tutorials/index",
        "formats/index",
        "glossary/index",
        "troubleshooting/index",
        "frequently_asked_questions/index",
    ]
    for document in required:
        assert (DOCS / f"{document}.rst").is_file()
        assert document in index


def test_homepage_contains_a_runnable_first_analysis():
    index = (DOCS / "index.rst").read_text()
    assert "PhyKIT\n======" in index
    assert "python -m pip install phykit" in index
    assert "phykit alignment_length" in index
    assert "Expected output" in index


def test_shared_format_contract_covers_core_inputs():
    formats = (DOCS / "formats" / "index.rst").read_text()
    for heading in (
        "Multiple sequence alignments",
        "Phylogenetic trees",
        "Taxon names",
        "Single-trait tables",
        "Multi-trait tables",
        "Missing values and missing taxa",
        "Paths and list files",
    ):
        assert heading in formats


def test_every_image_and_image_substitution_has_alt_text():
    image_pattern = re.compile(r"^\.\. (?:\|[^|]+\| )?image::", re.MULTILINE)
    for page in DOCS.rglob("*.rst"):
        lines = page.read_text().splitlines()
        for index, line in enumerate(lines):
            if image_pattern.match(line):
                options = lines[index + 1:index + 6]
                assert any(option.strip().startswith(":alt:") for option in options), page
