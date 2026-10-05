"""Characterize the existing headerless format without tightening validation."""

import math

import pytest

from phykit.errors import PhykitUserError
from phykit.services.tree.cont_map import ContMap
from phykit.services.tree.phenogram import Phenogram


@pytest.fixture(params=["shared", ContMap, Phenogram], ids=["shared", "cont_map", "phenogram"])
def parse_traits(request):
    if request.param == "shared":
        from phykit.helpers.trait_parsing import parse_single_trait_file

        return parse_single_trait_file
    service = request.param.__new__(request.param)
    return service._parse_single_trait_data


@pytest.mark.parametrize("tips", [["C", "A", "B"], ["A", "B", "C"], ["A", "A", "B", "C"]])
def test_all_shared_preserves_file_order(parse_traits, tmp_path, capsys, tips):
    path = tmp_path / "traits.tsv"
    path.write_text("  # comment\n\n C\t3\nA\t 1.5 \nB\t-2e1\n")
    result = parse_traits(str(path), tips)
    assert list(result.items()) == [("C", 3.0), ("A", 1.5), ("B", -20.0)]
    assert capsys.readouterr() == ("", "")


@pytest.mark.parametrize("service_class", [ContMap, Phenogram])
def test_service_method_delegates_to_shared_parser(service_class, monkeypatch):
    from phykit.helpers import trait_parsing

    tips = ["A", "B", "C"]
    expected = {"A": 1.0, "B": 2.0, "C": 3.0}
    calls = []

    def parser(path, tree_tips):
        calls.append((path, tree_tips))
        return expected

    monkeypatch.setattr(trait_parsing, "parse_single_trait_file", parser)
    service = service_class.__new__(service_class)
    assert service._parse_single_trait_data("traits.tsv", tips) is expected
    assert calls == [("traits.tsv", tips)]
    assert calls[0][1] is tips


def test_duplicate_rows_overwrite_without_moving_taxon(parse_traits, tmp_path, capsys):
    path = tmp_path / "traits.tsv"
    path.write_text("C\t3\nA\t1\nB\t2\nA\t9\n")
    assert list(parse_traits(str(path), ["A", "B", "C"]).items()) == [
        ("C", 3.0), ("A", 9.0), ("B", 2.0),
    ]
    assert capsys.readouterr() == ("", "")


def test_nonfinite_values_remain_accepted(parse_traits, tmp_path, capsys):
    path = tmp_path / "traits.tsv"
    path.write_text("A\tnan\nB\tinf\nC\t-inf\n")
    result = parse_traits(str(path), ["A", "B", "C"])
    assert math.isnan(result["A"])
    assert result["B"] == math.inf
    assert result["C"] == -math.inf
    assert capsys.readouterr() == ("", "")


def test_taxon_whitespace_and_case_are_not_normalized(parse_traits, tmp_path):
    path = tmp_path / "traits.tsv"
    path.write_text(" A \t1\na\t2\nC\t3\n")
    assert list(parse_traits(str(path), ["A ", "a", "C"])) == ["A ", "a", "C"]


def test_partial_overlap_preserves_set_iteration_and_sorted_warnings(parse_traits, tmp_path, capsys):
    path = tmp_path / "traits.tsv"
    path.write_text("C\t3\nA\t1\nY\t8\nB\t2\nX\t9\n")
    tips = ["D", "B", "E", "A", "C"]
    result = parse_traits(str(path), tips)
    shared = set(tips) & set(["C", "A", "Y", "B", "X"])
    assert list(result) == list(shared)
    assert result == {"A": 1.0, "B": 2.0, "C": 3.0}
    assert capsys.readouterr() == ("", (
        "Warning: 2 taxa in tree but not in trait file: D, E\n"
        "Warning: 2 taxa in trait file but not in tree: X, Y\n"
    ))


@pytest.mark.parametrize("text,messages", [
    ("# first\n\nA\t1\textra\n", [
        "Line 3 in trait file has 3 columns; expected 2.",
        "Each line should be: taxon_name<tab>trait_value",
    ]),
    ("A\t\n", [
        "Line 1 in trait file has 1 columns; expected 2.",
        "Each line should be: taxon_name<tab>trait_value",
    ]),
    ("A\t1\tx\ty\n", [
        "Line 1 in trait file has 4 columns; expected 2.",
        "Each line should be: taxon_name<tab>trait_value",
    ]),
    ("# comment\nA\tnot-a-number\n", [
        "Non-numeric trait value 'not-a-number' for taxon 'A' on line 2.",
    ]),
    ("taxon\ttrait\nA\t1\n", [
        "Non-numeric trait value 'trait' for taxon 'taxon' on line 1.",
    ]),
])
def test_exact_parse_errors(parse_traits, tmp_path, capsys, text, messages):
    path = tmp_path / "traits.tsv"
    path.write_text(text)
    with pytest.raises(PhykitUserError) as error:
        parse_traits(str(path), ["A", "B", "C"])
    assert error.value.code == 2
    assert error.value.messages == messages
    assert capsys.readouterr() == ("", "")


@pytest.mark.parametrize("text,tips,shared,warnings", [
    ("# no rows\n\n", [], 0, ""),
    ("", ["C", "A", "B"], 0,
     "Warning: 3 taxa in tree but not in trait file: A, B, C\n"),
    ("A\t1\nB\t2\nX\t3\n", ["A", "B", "C"], 2,
     "Warning: 1 taxa in tree but not in trait file: C\n"
     "Warning: 1 taxa in trait file but not in tree: X\n"),
    ("A\t1\nB\t2\n", ["A", "A", "B"], 2, ""),
])
def test_insufficient_overlap_errors_follow_warnings(parse_traits, tmp_path, capsys, text, tips, shared, warnings):
    path = tmp_path / "traits.tsv"
    path.write_text(text)
    with pytest.raises(PhykitUserError) as error:
        parse_traits(str(path), tips)
    assert error.value.code == 2
    assert error.value.messages == [
        f"Only {shared} shared taxa between tree and trait file.",
        "At least 3 shared taxa are required.",
    ]
    assert capsys.readouterr() == ("", warnings)


def test_missing_file_error_is_unchanged(parse_traits, tmp_path, capsys):
    path = tmp_path / "missing.tsv"
    with pytest.raises(PhykitUserError) as error:
        parse_traits(str(path), ["A", "B", "C"])
    assert error.value.code == 2
    assert error.value.messages == [
        f"{path} corresponds to no such file or directory.",
        "Please check filename and pathing",
    ]
    assert capsys.readouterr() == ("", "")
