from io import StringIO

import pytest
from Bio import Phylo

from phykit.errors import PhykitUserError
from phykit.helpers.topology_landscape import (
    classify_tree,
    groups_from_branch,
    read_tree,
    validate_groups,
)


def tree(text):
    return Phylo.read(StringIO(text), "newick")


GROUPS = {"A": ["a"], "B": ["b"], "C": ["c"], "D": ["d"]}


@pytest.mark.parametrize("newick,expected", [
    ("((a,b),(c,d));", "topology_1"),
    ("(a,(b,(c,d)));", "topology_1"),
    ("((a,c),(b,d));", "topology_2"),
    ("((a,d),(b,c));", "topology_3"),
    ("(a,b,c,d);", "unresolved"),
    ("(a,b,c);", "insufficient_sampling"),
    ("((a,b,x),(c,d,y));", "topology_1"),
])
def test_classification(newick, expected):
    assert classify_tree(tree(newick), GROUPS)["classification"] == expected


@pytest.mark.parametrize("newick,expected", [
    ("(((a,a2),b),(c,d));", "topology_1"),
    ("((a,b),(a2,c,d));", "incompatible_groups"),
    ("(a,a2,b,(c,d));", "unresolved"),
    ("((a,a2),b,c,d);", "unresolved"),
])
def test_all_representatives(newick, expected):
    groups = dict(GROUPS, A=["a", "a2"])
    assert classify_tree(tree(newick), groups)["classification"] == expected


def test_root_and_order_invariance():
    original = tree("(((a,a2),b),(c,(d,d2)));")
    groups = dict(GROUPS, A=["a", "a2"], D=["d", "d2"])
    for tip in [t.name for t in original.get_terminals()]:
        original.root_with_outgroup(tip)
        original.ladderize(reverse=True)
        assert classify_tree(original, groups)["classification"] == "topology_1"


def test_support_policies_and_root_duplicate():
    t = tree("((a,b)95,(c,d)60);")
    assert classify_tree(t, GROUPS)["resolution_support"] == 60
    assert classify_tree(t, GROUPS, 70)["classification"] == "unresolved"
    t = tree("((a,b),(c,d));")
    assert classify_tree(t, GROUPS, 70)["classification"] == "unresolved"
    assert classify_tree(t, GROUPS, 70, missing_support="keep")["classification"] == "topology_1"
    with pytest.raises(PhykitUserError):
        classify_tree(t, GROUPS, 70, missing_support="error")
    assert classify_tree(tree("((a,b)0.9,(c,d)0.9);"), GROUPS, .8, 1)["classification"] == "topology_1"


@pytest.mark.parametrize("kwargs", [
    {"min_support": -1}, {"min_support": float("nan")},
    {"min_support": 101}, {"support_scale": 10}, {"missing_support": "ignore"},
])
def test_bad_support_arguments(kwargs):
    with pytest.raises(PhykitUserError):
        classify_tree(tree("((a,b),(c,d));"), GROUPS, **kwargs)


@pytest.mark.parametrize("groups", [
    {}, dict(GROUPS, A=[]), dict(GROUPS, A=["b"]),
    dict(GROUPS, A=["a", "a"]), dict(GROUPS, A="a"),
    dict(GROUPS, A=[None]),
])
def test_invalid_groups(groups):
    with pytest.raises(PhykitUserError):
        validate_groups(groups)


def test_invalid_trees(tmp_path):
    for content in ("(a,a,b,c);", "(,a,b,c);", "(a,b);\n(c,d);"):
        path = tmp_path / "tree"
        path.write_text(content)
        with pytest.raises(PhykitUserError):
            read_tree(path)
    with pytest.raises(PhykitUserError):
        read_tree(tmp_path / "absent")
    with pytest.raises(PhykitUserError):
        classify_tree(tree("(a,a,b,c);"), GROUPS)
    with pytest.raises(PhykitUserError):
        classify_tree(tree("((a,b)101,(c,d));"), GROUPS)


def test_reference_branch_root_invariance():
    t = tree("(((a,a2),b),(c,(d,d2)));")
    expected = {"A": frozenset(["a", "a2"]), "B": frozenset(["b"]),
                    "C": frozenset(["c"]), "D": frozenset(["d", "d2"])}
    for tip in [x.name for x in t.get_terminals()]:
        t.root_with_outgroup(tip)
        assert groups_from_branch(t, ["a", "a2", "b"]) == expected


@pytest.mark.parametrize("newick,taxa", [
    ("(a,b,c,d);", ["a"]), ("((a,b),(c,d));", ["a", "c"]),
    ("((a,b,c),(d,e));", ["a", "b", "c"]),
])
def test_bad_reference_branch(newick, taxa):
    with pytest.raises(PhykitUserError):
        groups_from_branch(tree(newick), taxa)
