from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import check_wheel_entry_points


def test_resolve_wheel_accepts_file_and_single_wheel_directory(tmp_path):
    wheel = tmp_path / "phykit-1-py3-none-any.whl"
    wheel.touch()

    assert check_wheel_entry_points.resolve_wheel(wheel) == wheel
    assert check_wheel_entry_points.resolve_wheel(tmp_path) == wheel


@pytest.mark.parametrize("wheel_count", [0, 2])
def test_resolve_wheel_rejects_ambiguous_directory(tmp_path, wheel_count):
    for index in range(wheel_count):
        (tmp_path / f"phykit-{index}-py3-none-any.whl").touch()

    with pytest.raises(ValueError, match=f"found {wheel_count}"):
        check_wheel_entry_points.resolve_wheel(tmp_path)


def test_installed_target_probe_runs_outside_checkout(mocker):
    completed = SimpleNamespace(returncode=0, stdout="{}", stderr="")
    run = mocker.patch.object(
        check_wheel_entry_points.subprocess, "run", return_value=completed
    )

    assert check_wheel_entry_points.installed_target_load_errors(Path("python")) == {}
    assert Path(run.call_args.kwargs["cwd"]) != Path.cwd()


def test_smoke_installed_targets_combines_failures(mocker):
    mocker.patch.object(
        check_wheel_entry_points,
        "installed_target_load_errors",
        return_value={"pk_one": "import failed"},
    )
    mocker.patch.object(
        check_wheel_entry_points,
        "canonical_help_errors",
        return_value={"pk_two": "exit 1"},
    )

    assert check_wheel_entry_points.smoke_installed_targets(Path("python")) == {
        "installed target pk_one": "import failed",
        "installed help pk_two": "exit 1",
    }
