import json
import re
import subprocess
import sys
from pathlib import Path

from phykit.cli_registry import COMMAND_IDENTITIES, PUBLIC_COMMAND_TO_HANDLER


ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
COMMANDS_DIR = DOCS / "reference" / "commands"
SPEC_PATH = DOCS / "_data" / "command_spec.json"
CATALOG_PATH = DOCS / "_static" / "commands.json"
USAGE_PATH = DOCS / "usage" / "index.rst"


def load_commands(path):
    return json.loads(path.read_text())["commands"]


def test_registry_spec_catalog_and_pages_have_identical_command_coverage():
    registry = {identity.canonical for identity in COMMAND_IDENTITIES}
    spec = {command["canonical"] for command in load_commands(SPEC_PATH)}
    catalog = {command["canonical"] for command in load_commands(CATALOG_PATH)}
    pages = {path.stem for path in COMMANDS_DIR.glob("*.rst")}
    assert registry == spec == catalog == pages


def test_every_command_has_one_canonical_page_and_generated_interface():
    for identity in COMMAND_IDENTITIES:
        page = COMMANDS_DIR / f"{identity.canonical}.rst"
        content = page.read_text()
        label = f".. _cmd-{identity.canonical}:"
        include = f".. include:: /_generated/commands/{identity.canonical}.inc"
        assert content.count(label) == 1
        assert content.count(include) == 1
        assert "Command identity\n----------------" in content
        assert "Runtime interface\n-----------------" in content
        assert "Guidance, interpretation, and examples" in content


def test_spec_identity_and_entry_points_match_runtime_registry():
    spec = {command["canonical"]: command for command in load_commands(SPEC_PATH)}
    for identity in COMMAND_IDENTITIES:
        command = spec[identity.canonical]
        assert command["handler"] == identity.handler
        assert command["aliases"] == list(identity.aliases)
        assert command["entry_points"] == list(identity.entry_points)
        for name in (identity.canonical, *identity.aliases):
            assert PUBLIC_COMMAND_TO_HANDLER[name] == identity.handler


def test_spec_contains_complete_documentation_contract():
    required = {
        "purpose",
        "use_cases",
        "input_formats",
        "outputs",
        "errors",
        "limitations",
        "scientific_interpretation",
        "validation",
        "citations",
        "examples",
        "related_commands",
        "related_tutorials",
    }
    for command in load_commands(SPEC_PATH):
        assert command["summary"].strip()
        assert command["categories"]
        assert set(command["documentation"]) == required


def test_catalog_contains_runtime_argument_contracts():
    for command in load_commands(CATALOG_PATH):
        assert "arguments" in command
        for argument in command["arguments"]:
            assert set(argument) == {
                "name",
                "flags",
                "positional",
                "required",
                "type",
                "nargs",
                "default",
                "choices",
                "metavar",
                "action",
            }


def test_generated_catalog_is_current():
    result = subprocess.run(
        [sys.executable, "scripts/build_docs_catalog.py", "--check"],
        cwd=ROOT,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_legacy_usage_fragments_are_preserved():
    content = USAGE_PATH.read_text()
    fragments = set(re.findall(r'<span id="cmd-([a-z0-9_]+)"></span>', content))
    assert fragments == {identity.canonical for identity in COMMAND_IDENTITIES}
