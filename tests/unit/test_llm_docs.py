import json
import subprocess
import sys
from pathlib import Path

from phykit.cli_registry import COMMAND_IDENTITIES


ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"


def test_llm_artifacts_cover_every_canonical_command():
    index = (DOCS / "llms.txt").read_text()
    full = (DOCS / "llms-full.txt").read_text()
    for identity in COMMAND_IDENTITIES:
        url = f"reference/commands/{identity.canonical}.html"
        assert url in index
        assert f"## {identity.canonical}\n" in full


def test_llm_artifacts_are_exposed_at_documentation_root():
    config = (DOCS / "conf.py").read_text()
    assert 'html_extra_path = ["llms.txt", "llms-full.txt"]' in config


def test_retrieval_benchmark_has_at_least_fifty_gold_questions():
    benchmark = json.loads((DOCS / "_data" / "retrieval_benchmark.json").read_text())
    questions = benchmark["questions"]
    assert len(questions) >= 50
    assert len({item["id"] for item in questions}) == len(questions)
    for item in questions:
        assert item["question"].strip()
        assert item["answer"].strip()
        assert item["answer_terms"]
        assert item["sources"]


def test_retrieval_benchmark_meets_accuracy_thresholds():
    result = subprocess.run(
        [sys.executable, "scripts/benchmark_docs_retrieval.py"],
        cwd=ROOT,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "gold-answer factual support: 100.0%" in result.stdout
    assert "canonical source in top five: 100.0%" in result.stdout
