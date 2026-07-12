import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
TUTORIALS = DOCS / "tutorials" / "pages"
INDEX = DOCS / "tutorials" / "index.rst"
MANIFEST = DOCS / "_data" / "tutorial_smoke.json"


def test_tutorial_pages_and_smoke_manifest_cover_all_tutorials():
    pages = sorted(TUTORIALS.glob("*.rst"))
    manifest = json.loads(MANIFEST.read_text())["tutorials"]
    assert len(pages) == 21
    assert [item["number"] for item in manifest] == list(range(1, 22))


def test_every_tutorial_has_the_standard_contract():
    required = [
        "Objectives\n----------",
        "Prerequisites and working directory",
        "Related command references",
        "Workflow\n--------",
        "Expected artifacts",
        "Troubleshooting\n---------------",
    ]
    for number, page in enumerate(sorted(TUTORIALS.glob("*.rst")), start=1):
        content = page.read_text()
        assert f".. _tutorial-{number:02d}:" in content
        for heading in required:
            assert heading in content
        assert "tests/sample_files" not in content
        assert "raw.githubusercontent.com/JLSteenwyk/PhyKIT/master" not in content


def test_all_tutorial_downloads_resolve_to_versioned_data():
    pattern = re.compile(r":download:`[^`]+ </data/([^>]+)>`")
    for page in TUTORIALS.glob("*.rst"):
        downloads = pattern.findall(page.read_text())
        assert downloads, f"{page.name} has no downloadable data"
        for filename in downloads:
            assert (DOCS / "data" / filename).is_file(), (page.name, filename)


def test_legacy_tutorial_fragments_are_present():
    content = INDEX.read_text()
    fragments = re.findall(r'<span id="([a-z0-9-]+)"></span>', content)
    assert len(fragments) == 21
    assert len(set(fragments)) == 21


def test_commands_and_expected_output_use_separate_blocks():
    suspicious_output = re.compile(
        r"(?m)^\s+(?:mean:|median:|p-value:|chi-squared:|total genes:|"
        r"\d+(?:\.\d+)?(?:\s+\d+(?:\.\d+)?)+)$"
    )
    for page in TUTORIALS.glob("*.rst"):
        content = page.read_text()
        for block in re.findall(
            r"\.\. code-block:: shell\n\n((?:   .*\n|\n)+)", content
        ):
            assert not suspicious_output.search(block), page.name


def test_every_documentation_image_has_alt_text():
    for page in DOCS.rglob("*.rst"):
        lines = page.read_text().splitlines()
        for index, line in enumerate(lines):
            if line.startswith(".. image::"):
                options = lines[index + 1:index + 5]
                assert any(option.strip().startswith(":alt:") for option in options), page
