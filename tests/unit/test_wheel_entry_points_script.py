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


def test_install_smoke_does_not_require_checkout_dependencies(
    mocker, monkeypatch, tmp_path
):
    wheel = tmp_path / "phykit-1-py3-none-any.whl"
    mocker.patch.object(
        check_wheel_entry_points, "resolve_wheel", return_value=wheel
    )
    archive = mocker.patch.object(check_wheel_entry_points, "ZipFile")
    archive.return_value.__enter__.return_value.namelist.return_value = [
        "phykit.dist-info/entry_points.txt"
    ]
    archive.return_value.__enter__.return_value.read.return_value = (
        b"[console_scripts]\n"
    )
    expected = {"phykit": "phykit.phykit:main"}
    expected.update(
        {
            f"pk_{command}": f"phykit.phykit:{handler}"
            for command, handler in (
                check_wheel_entry_points.PUBLIC_COMMAND_TO_HANDLER.items()
            )
        }
    )
    config = mocker.patch.object(check_wheel_entry_points, "configparser")
    config.ConfigParser.return_value.__getitem__.return_value = expected
    local_check = mocker.patch.object(
        check_wheel_entry_points,
        "target_load_errors",
        side_effect=AssertionError("checkout targets must not load"),
    )
    mocker.patch.object(
        check_wheel_entry_points,
        "create_installed_wheel_environment",
        return_value=Path("python"),
    )
    mocker.patch.object(
        check_wheel_entry_points, "smoke_installed_targets", return_value={}
    )
    monkeypatch.setattr(
        check_wheel_entry_points.sys,
        "argv",
        ["check_wheel_entry_points.py", str(tmp_path), "--install-smoke"],
    )

    assert check_wheel_entry_points.main() == 0
    local_check.assert_not_called()
