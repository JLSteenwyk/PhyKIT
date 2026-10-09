"""
CLI wiring for the global --threads option
"""

import os

import pytest

import phykit.phykit as phykit_module
from phykit.helpers import threads


@pytest.fixture(autouse=True)
def clean_thread_env(monkeypatch):
    # setenv first so monkeypatch restores (or removes) every variable the
    # CLI writes directly to os.environ.
    for name in (threads.PHYKIT_THREADS_ENV,) + threads.THREAD_ENV_VARS:
        monkeypatch.setenv(name, "")
        monkeypatch.delenv(name)


@pytest.fixture
def captured_alignment_length(monkeypatch):
    calls = {}

    def fake_method(argv):
        calls["argv"] = argv

    monkeypatch.setattr(
        phykit_module.Phykit, "alignment_length", staticmethod(fake_method)
    )
    return calls


@pytest.mark.parametrize(
    "argv",
    [
        ["phykit", "alignment_length", "aln.fa", "--threads", "2"],
        ["phykit", "--threads", "2", "alignment_length", "aln.fa"],
        ["phykit", "alignment_length", "--threads=2", "aln.fa"],
    ],
)
def test_phykit_strips_threads_and_applies_limit(
    monkeypatch, captured_alignment_length, argv
):
    monkeypatch.setattr(phykit_module.sys, "argv", argv)

    phykit_module.Phykit()

    assert captured_alignment_length["argv"] == ["aln.fa"]
    assert os.environ[threads.PHYKIT_THREADS_ENV] == "2"
    assert os.environ["OMP_NUM_THREADS"] == "2"


def test_phykit_alias_receives_argv_without_threads(
    monkeypatch, captured_alignment_length
):
    monkeypatch.setattr(
        phykit_module.sys, "argv", ["phykit", "aln_len", "aln.fa", "--threads", "3"]
    )

    phykit_module.Phykit()

    assert captured_alignment_length["argv"] == ["aln.fa"]
    assert os.environ[threads.PHYKIT_THREADS_ENV] == "3"


def test_pk_entry_point_strips_threads_and_applies_limit(
    monkeypatch, captured_alignment_length
):
    monkeypatch.setattr(
        phykit_module.sys, "argv", ["pk_alignment_length", "aln.fa", "--threads=4"]
    )

    phykit_module.alignment_length()

    assert captured_alignment_length["argv"] == ["aln.fa"]
    assert os.environ[threads.PHYKIT_THREADS_ENV] == "4"


def test_pk_entry_point_with_explicit_argv_strips_threads(monkeypatch):
    calls = {}
    monkeypatch.setattr(
        phykit_module.Phykit,
        "topology_landscape",
        staticmethod(lambda argv: calls.setdefault("argv", argv)),
    )

    phykit_module.topology_landscape(["-t", "trees.txt", "--threads", "2"])

    assert calls["argv"] == ["-t", "trees.txt"]


@pytest.mark.parametrize(
    "argv",
    [
        ["phykit", "alignment_length", "aln.fa", "--threads", "0"],
        ["phykit", "alignment_length", "aln.fa", "--threads"],
    ],
)
def test_invalid_threads_value_exits_with_message(
    monkeypatch, capsys, captured_alignment_length, argv
):
    monkeypatch.setattr(phykit_module.sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        phykit_module.Phykit()

    assert excinfo.value.code == 2
    assert "--threads requires a positive integer" in capsys.readouterr().err
    assert "argv" not in captured_alignment_length


def test_top_level_help_documents_threads(monkeypatch, capsys):
    monkeypatch.setattr(phykit_module.sys, "argv", ["phykit"])

    with pytest.raises(SystemExit):
        phykit_module.Phykit()

    out = capsys.readouterr().out
    assert "--threads" in out
    assert "PHYKIT_THREADS" in out
