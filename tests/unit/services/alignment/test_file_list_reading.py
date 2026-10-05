"""Characterize shared file-list semantics independently of file contents."""

import importlib
import os

import pytest

from phykit.errors import PhykitUserError


@pytest.fixture(params=[("taxon_groups", "TaxonGroups"),
                       ("occupancy_filter", "OccupancyFilter")])
def read_list(request):
    name, class_name = request.param
    module = importlib.import_module(f"phykit.services.alignment.{name}")
    cls = getattr(module, class_name)
    return cls.__new__(cls)._read_file_list


def test_relative_list_preserves_duplicates_and_parent_segments(read_list, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "inputs").mkdir()
    path = tmp_path / "inputs" / "files.txt"
    path.write_text(" # comment\n\n . \n./\nfile.fa\nfile.fa\n"
                    "./nested/../file.fa\nwith spaces.fa\nname#part.fa".replace("/", os.sep))
    prefix = "inputs" + os.sep
    assert read_list(os.path.join("inputs", "files.txt")) == [
        "inputs", "inputs", prefix + "file.fa", prefix + "file.fa",
        prefix + os.sep.join(["nested", "..", "file.fa"]),
        prefix + "with spaces.fa", prefix + "name#part.fa",
    ]


def test_current_directory_list_does_not_add_dot_prefix(read_list, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "files.txt").write_text(".\nfile.fa\n")
    assert read_list("files.txt") == [".", "file.fa"]


def test_empty_list_is_returned_to_caller(read_list, tmp_path):
    path = tmp_path / "files.txt"
    path.write_text("# comment\n\n")
    assert read_list(path) == []


def test_missing_list_preserves_original_path_in_error(read_list, tmp_path):
    path = str(tmp_path) + os.sep + "." + os.sep + "missing.txt"
    with pytest.raises(PhykitUserError) as error:
        read_list(path)
    assert error.value.code == 2
    assert error.value.messages == [f"{path} corresponds to no such file or directory."]


@pytest.mark.skipif(os.sep != "/", reason="POSIX absolute-path spelling")
def test_absolute_paths_keep_double_slash_root(read_list, tmp_path):
    path = tmp_path / "files.txt"
    path.write_text("//server//./file.fa\n///server//file.fa\n/dir/../file.fa\n")
    assert read_list(path) == ["//server/file.fa", "/server/file.fa", "/dir/../file.fa"]
