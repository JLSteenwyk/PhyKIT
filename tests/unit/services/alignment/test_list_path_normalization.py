"""Preserve list-path spelling without resolving parents or symlinks."""

import importlib
import os

import pytest


@pytest.fixture(params=["_file_list", "taxon_groups", "occupancy_filter"])
def normalize(request):
    module = importlib.import_module(f"phykit.services.alignment.{request.param}")
    return module._normalize_list_path


@pytest.mark.parametrize("path,expected", [
    ("", ""),
    ("file.fa", "file.fa"),
    (".", "."),
    ("./", "."),
    ("././file.fa", "file.fa"),
    ("dir//./file.fa", "dir/file.fa"),
    ("dir/.", "dir"),
    ("dir/../file.fa", "dir/../file.fa"),
    ("./dir/../file.fa", "dir/../file.fa"),
    ("dir//", "dir"),
    ("dir/", "dir/"),
    ("dir with spaces/./file.fa", "dir with spaces/file.fa"),
])
def test_relative_path_spelling(normalize, path, expected):
    assert normalize(path.replace("/", os.sep)) == expected.replace("/", os.sep)


@pytest.mark.skipif(os.sep != "/", reason="POSIX root spelling")
@pytest.mark.parametrize("path,expected", [
    ("/", "/"),
    ("//", "//"),
    ("///", "/"),
    ("//server//./file.fa", "//server/file.fa"),
    ("///server//./file.fa", "/server/file.fa"),
    ("/dir/.././file.fa", "/dir/../file.fa"),
])
def test_posix_root_spelling(normalize, path, expected):
    assert normalize(path) == expected


def test_simple_path_returns_original_string(normalize):
    path = os.sep.join(["directory", "file.fa"])
    assert normalize(path) is path


@pytest.mark.parametrize("name", ["taxon_groups", "occupancy_filter"])
def test_services_alias_shared_function_without_wrapper(name):
    from phykit.services.alignment._file_list import _normalize_list_path

    module = importlib.import_module(f"phykit.services.alignment.{name}")
    assert module._normalize_list_path is _normalize_list_path
