#!/usr/bin/env python3
"""Build the traceable command and service test-coverage inventory."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import phykit.phykit as cli
from phykit.cli_registry import COMMAND_IDENTITIES, PUBLIC_COMMAND_TO_HANDLER
from phykit.service_factories import SERVICE_FACTORIES


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "tests" / "coverage_inventory.json"
UNIT_ROOT = ROOT / "tests" / "unit"
INTEGRATION_ROOT = ROOT / "tests" / "integration"
TUTORIAL_MANIFEST = ROOT / "docs" / "_data" / "tutorial_smoke.json"

UNIT_OVERRIDES = {
    "phykit/services/alignment/_fasta.py": [
        "tests/unit/services/alignment/test_fasta_helpers.py"
    ],
    "phykit/services/alignment/base.py": [
        "tests/unit/services/alignment/test_alignment_length.py",
        "tests/unit/services/alignment/test_create_concatenation_matrix.py",
    ],
    "phykit/services/base.py": ["tests/unit/services/test_base.py"],
    "phykit/services/tree/base.py": [
        "tests/unit/services/tree/test_base_tree_service.py"
    ],
}

INTEGRATION_OVERRIDES = {
    "phykit/services/alignment/alignment_length_no_gaps.py": [
        "tests/integration/alignment/test_alignment_length_integration_no_gaps.py"
    ],
    "phykit/services/tree/prune_tree.py": [
        "tests/integration/tree/test_prune_integration.py"
    ],
}


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def test_files(root: Path) -> list[Path]:
    return sorted(root.rglob("test_*.py"))


def normalized_test_name(path: Path, *, integration: bool) -> str:
    name = path.stem.removeprefix("test_")
    if integration:
        return name.replace("_integration", "")
    return name.removesuffix("_service")


def named_tests(
    module_path: str,
    tests: list[Path],
    overrides: dict,
    *,
    integration: bool,
) -> list[str]:
    if module_path in overrides:
        return overrides[module_path]
    stem = Path(module_path).stem
    accepted = {stem, f"{stem}_invariant"}
    matches = [
        relative(path)
        for path in tests
        if normalized_test_name(path, integration=integration) in accepted
    ]
    return sorted(matches)


def command_service(identity) -> tuple[str | None, str | None]:
    if identity.handler == "version":
        return None, None
    method = getattr(cli.Phykit, identity.handler)
    factories = sorted(
        set(method.__code__.co_names).intersection(SERVICE_FACTORIES)
    )
    if len(factories) != 1:
        raise RuntimeError(
            f"expected one service factory for {identity.handler}, found {factories}"
        )
    factory_name = factories[0]
    factory = SERVICE_FACTORIES[factory_name]
    return factory_name, factory.module_path.replace(".", "/") + ".py"


def tutorial_coverage() -> dict[str, list[int]]:
    manifest = json.loads(TUTORIAL_MANIFEST.read_text())
    canonical_by_handler = {
        identity.handler: identity.canonical for identity in COMMAND_IDENTITIES
    }
    tutorials: dict[str, list[int]] = {}
    for tutorial in manifest["tutorials"]:
        public_name = tutorial["command"][0]
        handler = PUBLIC_COMMAND_TO_HANDLER[public_name]
        canonical = canonical_by_handler[handler]
        tutorials.setdefault(canonical, []).append(tutorial["number"])
    return tutorials


def build_inventory() -> dict:
    unit_files = test_files(UNIT_ROOT)
    integration_files = test_files(INTEGRATION_ROOT)
    tutorials = tutorial_coverage()
    commands = []
    commands_by_module: dict[str, list[str]] = {}

    for identity in COMMAND_IDENTITIES:
        factory_name, module_path = command_service(identity)
        if module_path:
            unit = named_tests(
                module_path, unit_files, UNIT_OVERRIDES, integration=False
            )
            integration = named_tests(
                module_path,
                integration_files,
                INTEGRATION_OVERRIDES,
                integration=True,
            )
            commands_by_module.setdefault(module_path, []).append(identity.canonical)
        else:
            unit = [
                "tests/unit/test_cli_registry.py",
                "tests/unit/test_main_entrypoint.py",
            ]
            integration = []
        tutorial_numbers = tutorials.get(identity.canonical, [])
        wrapper_callable = callable(getattr(cli, identity.handler, None))
        if not wrapper_callable:
            status = "defect"
            disposition = "Registered console entry-point target is not callable."
        elif not unit:
            status = "review_required"
            disposition = "No direct unit-test module was located."
        elif integration or tutorial_numbers:
            status = "unit_and_workflow"
            disposition = "Direct unit coverage plus an integration or wheel workflow."
        else:
            status = "unit_only_documented"
            disposition = (
                "Direct unit tests exercise the service, including run/output paths; "
                "add integration coverage only for a concrete cross-layer contract."
            )
        commands.append(
            {
                "canonical": identity.canonical,
                "handler": identity.handler,
                "aliases": list(identity.aliases),
                "entry_points": list(identity.entry_points),
                "factory": factory_name,
                "service_module": module_path,
                "unit_tests": unit,
                "integration_tests": integration,
                "tutorial_smoke": tutorial_numbers,
                "parser_contract": "tests/unit/test_docs_category_coverage.py",
                "module_wrapper_callable": wrapper_callable,
                "status": status,
                "disposition": disposition,
            }
        )

    services = []
    service_files = sorted(
        path
        for path in (ROOT / "phykit" / "services").rglob("*.py")
        if path.name != "__init__.py"
    )
    for path in service_files:
        module_path = relative(path)
        unit = named_tests(
            module_path, unit_files, UNIT_OVERRIDES, integration=False
        )
        integration = named_tests(
            module_path,
            integration_files,
            INTEGRATION_OVERRIDES,
            integration=True,
        )
        commands_for_module = sorted(commands_by_module.get(module_path, []))
        if unit:
            disposition = "Direct unit coverage is present."
        else:
            disposition = (
                "Shared base behavior is exercised through downstream service tests; "
                "direct tests remain a documented follow-up where branch evidence warrants."
            )
        services.append(
            {
                "module": module_path,
                "commands": commands_for_module,
                "unit_tests": unit,
                "integration_tests": integration,
                "disposition": disposition,
            }
        )

    return {
        "schema_version": 1,
        "authoritative_sources": [
            "phykit/cli_registry.py",
            "phykit/phykit.py",
            "phykit/service_factories.py",
            "tests/unit",
            "tests/integration",
            "docs/_data/tutorial_smoke.json",
        ],
        "summary": {
            "canonical_commands": len(commands),
            "public_command_names": len(PUBLIC_COMMAND_TO_HANDLER),
            "service_modules": len(services),
            "commands_with_unit_tests": sum(bool(item["unit_tests"]) for item in commands),
            "commands_with_integration_tests": sum(
                bool(item["integration_tests"]) for item in commands
            ),
            "commands_with_tutorial_smoke": sum(
                bool(item["tutorial_smoke"]) for item in commands
            ),
            "uncallable_module_wrappers": sum(
                not item["module_wrapper_callable"] for item in commands
            ),
        },
        "commands": commands,
        "services": services,
    }


def render(data: dict) -> str:
    return json.dumps(data, indent=2, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    expected = render(build_inventory())
    if args.check:
        actual = OUTPUT.read_text() if OUTPUT.exists() else ""
        if actual != expected:
            print(f"test coverage inventory is stale: {relative(OUTPUT)}")
            return 1
        print(f"test coverage inventory is current: {relative(OUTPUT)}")
        return 0
    OUTPUT.write_text(expected)
    print(f"wrote {relative(OUTPUT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
