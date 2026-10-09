"""
Unit tests for the global --threads option helpers
"""

import subprocess
import sys

import pytest

from phykit.helpers import threads
from phykit.errors import PhykitUserError


@pytest.fixture(autouse=True)
def clean_thread_env(monkeypatch):
    # setenv first so monkeypatch restores (or removes) every variable the
    # code under test writes directly to os.environ.
    for name in (threads.PHYKIT_THREADS_ENV,) + threads.THREAD_ENV_VARS:
        monkeypatch.setenv(name, "")
        monkeypatch.delenv(name)


def test_module_import_does_not_import_heavy_dependencies():
    code = """
import sys
import phykit.helpers.threads

assert "multiprocessing" not in sys.modules
assert "numpy" not in sys.modules
assert "concurrent.futures" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


class TestExtractThreadsOption:
    def test_returns_argv_unchanged_without_option(self):
        assert threads.extract_threads_option(["-a", "aln.fa"]) == (None, ["-a", "aln.fa"])

    def test_strips_space_separated_option(self):
        assert threads.extract_threads_option(["aln.fa", "--threads", "4"]) == (4, ["aln.fa"])

    def test_strips_equals_option(self):
        assert threads.extract_threads_option(["--threads=2", "aln.fa"]) == (2, ["aln.fa"])

    def test_last_occurrence_wins(self):
        value, rest = threads.extract_threads_option(
            ["--threads", "2", "x", "--threads=3"]
        )
        assert (value, rest) == (3, ["x"])

    def test_does_not_touch_similar_option_names(self):
        argv = ["--threshold", "0.5", "--threads-extra", "1"]
        assert threads.extract_threads_option(argv) == (None, argv)

    @pytest.mark.parametrize(
        "argv",
        [["--threads"], ["--threads", "0"], ["--threads", "-1"], ["--threads=abc"]],
    )
    def test_rejects_missing_or_invalid_values(self, argv):
        with pytest.raises(PhykitUserError) as excinfo:
            threads.extract_threads_option(argv)
        assert excinfo.value.code == 2
        assert "--threads" in excinfo.value.messages[0]


class TestApplyThreadLimit:
    def test_sets_phykit_and_numeric_library_env_vars(self, monkeypatch):
        import os

        monkeypatch.setenv("OMP_NUM_THREADS", "16")
        threads.apply_thread_limit(3)

        assert os.environ[threads.PHYKIT_THREADS_ENV] == "3"
        for name in threads.THREAD_ENV_VARS:
            assert os.environ[name] == "3"


class TestThreadLimit:
    def test_none_when_unset(self):
        assert threads.thread_limit() is None

    def test_reads_env_var(self, monkeypatch):
        monkeypatch.setenv(threads.PHYKIT_THREADS_ENV, "5")
        assert threads.thread_limit() == 5

    @pytest.mark.parametrize("value", ["", "0", "-2", "many"])
    def test_ignores_invalid_env_values(self, monkeypatch, value):
        monkeypatch.setenv(threads.PHYKIT_THREADS_ENV, value)
        assert threads.thread_limit() is None


class TestLimitWorkers:
    def test_returns_count_unchanged_without_limit(self):
        assert threads.limit_workers(8) == 8

    def test_caps_count_at_limit(self, monkeypatch):
        monkeypatch.setenv(threads.PHYKIT_THREADS_ENV, "2")
        assert threads.limit_workers(8) == 2

    def test_does_not_raise_count_above_requested(self, monkeypatch):
        monkeypatch.setenv(threads.PHYKIT_THREADS_ENV, "16")
        assert threads.limit_workers(4) == 4

    def test_never_returns_less_than_one(self):
        assert threads.limit_workers(0) == 1


class TestConfigureThreads:
    def test_flag_is_stripped_and_applied(self):
        import os

        rest = threads.configure_threads(["saturation", "--threads", "2", "-a", "x"])

        assert rest == ["saturation", "-a", "x"]
        assert os.environ[threads.PHYKIT_THREADS_ENV] == "2"
        assert os.environ["OMP_NUM_THREADS"] == "2"

    def test_env_var_alone_limits_numeric_libraries(self, monkeypatch):
        import os

        monkeypatch.setenv(threads.PHYKIT_THREADS_ENV, "3")
        assert threads.configure_threads(["x"]) == ["x"]
        assert os.environ["OMP_NUM_THREADS"] == "3"

    def test_env_var_does_not_override_explicit_library_settings(self, monkeypatch):
        import os

        monkeypatch.setenv(threads.PHYKIT_THREADS_ENV, "3")
        monkeypatch.setenv("OMP_NUM_THREADS", "1")
        threads.configure_threads(["x"])
        assert os.environ["OMP_NUM_THREADS"] == "1"

    def test_no_flag_and_no_env_leaves_environment_alone(self):
        import os

        assert threads.configure_threads(["x"]) == ["x"]
        assert threads.PHYKIT_THREADS_ENV not in os.environ
        for name in threads.THREAD_ENV_VARS:
            assert name not in os.environ
