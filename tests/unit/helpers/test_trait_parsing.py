import builtins
import subprocess
import sys

import pytest

from phykit.errors import PhykitUserError
from phykit.helpers.trait_parsing import (
    parse_multi_trait_file,
    response_predictor_arrays,
    subset_traits_to_ordered_shared_taxa,
    trait_column_from_rows,
    trait_matrix_from_rows,
)


def test_module_import_does_not_import_typing():
    code = """
import sys
import phykit.helpers.trait_parsing
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_parse_multi_trait_file_skips_comments_and_filters_shared_taxa(
    tmp_path, capsys
):
    trait_file = tmp_path / "traits.tsv"
    trait_file.write_text(
        "\n".join(
            [
                "   # ignored before header",
                "taxon\tbody_mass\tlength",
                "",
                "A\t1.0\t10.0",
                "\t# ignored between rows",
                "B\t2.0\t20.0",
                "C\t3.0\t30.0",
                "off_tree\t4.0\t40.0",
            ]
        )
        + "\n"
    )

    trait_names, traits = parse_multi_trait_file(
        str(trait_file), ["A", "B", "C", "missing"]
    )

    assert trait_names == ["body_mass", "length"]
    assert traits == {"A": [1.0, 10.0], "B": [2.0, 20.0], "C": [3.0, 30.0]}
    err = capsys.readouterr().err
    assert "1 taxa in tree but not in trait file: missing" in err
    assert "1 taxa in trait file but not in tree: off_tree" in err


def test_parse_multi_trait_file_all_shared_emits_no_warnings(tmp_path, capsys):
    trait_file = tmp_path / "traits.tsv"
    trait_file.write_text(
        "\n".join(
            [
                "taxon\tbody_mass\tlength",
                "A\t1.0\t10.0",
                "B\t2.0\t20.0",
                "C\t3.0\t30.0",
            ]
        )
        + "\n"
    )

    trait_names, traits = parse_multi_trait_file(str(trait_file), ["A", "B", "C"])

    assert trait_names == ["body_mass", "length"]
    assert traits == {"A": [1.0, 10.0], "B": [2.0, 20.0], "C": [3.0, 30.0]}
    assert capsys.readouterr().err == ""


def test_parse_multi_trait_file_four_trait_fast_path(tmp_path, capsys):
    trait_file = tmp_path / "traits.tsv"
    trait_file.write_text(
        "\n".join(
            [
                "taxon\ta\tb\tc\td",
                "A\t1.0\t2.0\t3.0\t4.0",
                "B\t5.0\t6.0\t7.0\t8.0",
                "C\t9.0\t10.0\t11.0\t12.0",
            ]
        )
        + "\n"
    )

    trait_names, traits = parse_multi_trait_file(str(trait_file), ["A", "B", "C"])

    assert trait_names == ["a", "b", "c", "d"]
    assert traits == {
        "A": [1.0, 2.0, 3.0, 4.0],
        "B": [5.0, 6.0, 7.0, 8.0],
        "C": [9.0, 10.0, 11.0, 12.0],
    }
    assert capsys.readouterr().err == ""


def test_parse_multi_trait_file_ordered_all_shared_skips_set_validation(
    tmp_path, monkeypatch, capsys
):
    trait_file = tmp_path / "traits.tsv"
    trait_file.write_text(
        "\n".join(
            [
                "taxon\tbody_mass\tlength",
                "A\t1.0\t10.0",
                "B\t2.0\t20.0",
                "C\t3.0\t30.0",
            ]
        )
        + "\n"
    )

    def fail_set(*_args, **_kwargs):
        raise AssertionError("ordered all-shared traits should skip set validation")

    monkeypatch.setattr(builtins, "set", fail_set)

    trait_names, traits = parse_multi_trait_file(str(trait_file), ["A", "B", "C"])

    assert trait_names == ["body_mass", "length"]
    assert traits == {"A": [1.0, 10.0], "B": [2.0, 20.0], "C": [3.0, 30.0]}
    assert capsys.readouterr().err == ""


def test_parse_multi_trait_file_preserves_logical_data_row_errors(tmp_path):
    trait_file = tmp_path / "traits.tsv"
    trait_file.write_text(
        "\n".join(
            [
                "taxon\tbody_mass\tlength",
                "A\t1.0\t10.0",
                "# ignored between rows",
                "B\t2.0",
            ]
        )
        + "\n"
    )

    with pytest.raises(PhykitUserError) as excinfo:
        parse_multi_trait_file(str(trait_file), ["A", "B", "C"])

    assert "Line 3 has 2 columns; expected 3." in excinfo.value.messages


def test_parse_multi_trait_file_reports_non_numeric_trait(tmp_path):
    trait_file = tmp_path / "traits.tsv"
    trait_file.write_text(
        "\n".join(
            [
                "taxon\tbody_mass\tlength",
                "A\t1.0\t10.0",
                "B\t2.0\tbad",
                "C\t3.0\t30.0",
            ]
        )
        + "\n"
    )

    with pytest.raises(PhykitUserError) as excinfo:
        parse_multi_trait_file(str(trait_file), ["A", "B", "C"])

    assert (
        "Non-numeric trait value 'bad' for taxon 'B' "
        "(trait 'length') on line 3."
    ) in excinfo.value.messages


def test_parse_multi_trait_file_reports_non_numeric_third_trait(tmp_path):
    trait_file = tmp_path / "traits.tsv"
    trait_file.write_text(
        "\n".join(
            [
                "taxon\tbody_mass\tlength\theight",
                "A\t1.0\t10.0\t100.0",
                "B\t2.0\t20.0\tbad",
                "C\t3.0\t30.0\t300.0",
            ]
        )
        + "\n"
    )

    with pytest.raises(PhykitUserError) as excinfo:
        parse_multi_trait_file(str(trait_file), ["A", "B", "C"])

    assert (
        "Non-numeric trait value 'bad' for taxon 'B' "
        "(trait 'height') on line 3."
    ) in excinfo.value.messages


def test_trait_matrix_from_rows_preserves_requested_order_and_numeric_dtype():
    import numpy as np

    traits = {
        "taxon_b": [2, 20.5],
        "taxon_a": [1, 10],
        "taxon_c": [3, 30],
    }

    matrix = trait_matrix_from_rows(traits, ["taxon_a", "taxon_b", "taxon_c"])

    assert matrix.dtype == np.dtype(float)
    assert matrix.tolist() == [[1.0, 10.0], [2.0, 20.5], [3.0, 30.0]]


def test_trait_column_from_rows_preserves_requested_order_and_numeric_dtype():
    import numpy as np

    traits = {
        "taxon_b": [2, 20.5, 200],
        "taxon_a": [1, 10.0, 100],
        "taxon_c": [3, 30.0, 300],
    }

    column = trait_column_from_rows(
        traits,
        ["taxon_a", "taxon_b", "taxon_c"],
        1,
    )

    assert column.dtype == np.dtype(float)
    assert column.tolist() == [10.0, 20.5, 30.0]


def test_response_predictor_arrays_builds_intercept_design_matrix():
    import numpy as np

    traits = {
        "taxon_b": [2, 20.5, 200],
        "taxon_a": [1, 10.0, 100],
        "taxon_c": [3, 30.0, 300],
    }

    y, X = response_predictor_arrays(
        traits,
        ["taxon_a", "taxon_b", "taxon_c"],
        1,
        [2, 0],
    )

    assert y.dtype == np.dtype(float)
    assert X.dtype == np.dtype(float)
    assert y.tolist() == [10.0, 20.5, 30.0]
    assert X.tolist() == [[1.0, 100.0, 1.0], [1.0, 200.0, 2.0], [1.0, 300.0, 3.0]]


def test_response_predictor_arrays_narrow_design_uses_single_row_pass(mocker):
    import numpy as np

    traits = {
        "taxon_b": [2, 20.5, 200],
        "taxon_a": [1, 10.0, 100],
        "taxon_c": [3, 30.0, 300],
    }
    mocked_asarray = mocker.patch.object(
        np,
        "asarray",
        side_effect=AssertionError("narrow selected columns should fill rows directly"),
    )
    mocked_fromiter = mocker.patch.object(
        np,
        "fromiter",
        side_effect=AssertionError("narrow selected columns should fill rows directly"),
    )

    y, X = response_predictor_arrays(
        traits,
        ["taxon_a", "taxon_b", "taxon_c"],
        1,
        [2],
    )

    mocked_asarray.assert_not_called()
    mocked_fromiter.assert_not_called()
    assert y.tolist() == [10.0, 20.5, 30.0]
    assert X.tolist() == [[1.0, 100.0], [1.0, 200.0], [1.0, 300.0]]


def test_response_predictor_arrays_wide_design_keeps_selected_matrix(mocker):
    import numpy as np

    traits = {
        "taxon_b": [2, 20.5, 200, 2000, 20000, 200000, 2_000_000],
        "taxon_a": [1, 10.0, 100, 1000, 10000, 100000, 1_000_000],
        "taxon_c": [3, 30.0, 300, 3000, 30000, 300000, 3_000_000],
    }
    asarray_spy = mocker.spy(np, "asarray")

    y, X = response_predictor_arrays(
        traits,
        ["taxon_a", "taxon_b", "taxon_c"],
        1,
        [0, 2, 3, 4, 5],
    )

    asarray_spy.assert_called_once()
    assert y.tolist() == [10.0, 20.5, 30.0]
    assert X.tolist() == [
        [1.0, 1.0, 100.0, 1000.0, 10000.0, 100000.0],
        [1.0, 2.0, 200.0, 2000.0, 20000.0, 200000.0],
        [1.0, 3.0, 300.0, 3000.0, 30000.0, 300000.0],
    ]


def test_response_predictor_arrays_handles_no_predictors():
    y, X = response_predictor_arrays(
        {"taxon_b": [2.0], "taxon_a": [1.0]},
        ["taxon_a", "taxon_b"],
        0,
        [],
    )

    assert y.tolist() == [1.0, 2.0]
    assert X.tolist() == [[1.0], [1.0]]


def test_response_predictor_arrays_intercept_only_avoids_fromiter(monkeypatch):
    import numpy as np

    traits = {
        "taxon_b": [2.0, 20.0],
        "taxon_a": [1.0, 10.0],
        "taxon_c": [3.0, 30.0],
    }

    def fail_fromiter(*_args, **_kwargs):
        raise AssertionError(
            "intercept-only GLM arrays should use list-backed array construction"
        )

    monkeypatch.setattr(np, "fromiter", fail_fromiter)

    y, X = response_predictor_arrays(
        traits,
        ["taxon_a", "taxon_b", "taxon_c"],
        1,
        [],
    )

    assert y.tolist() == [10.0, 20.0, 30.0]
    assert X.tolist() == [[1.0], [1.0], [1.0]]


def test_subset_traits_to_ordered_shared_taxa_returns_original_when_all_shared():
    traits = {"A": [1.0], "B": [2.0], "C": [3.0]}
    ordered_names = ["A", "B", "C"]

    subset, names = subset_traits_to_ordered_shared_taxa(
        traits, ordered_names, ["A", "B", "C"]
    )

    assert subset is traits
    assert names is ordered_names


def test_subset_traits_to_ordered_shared_taxa_filters_to_shared_order():
    traits = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}

    subset, names = subset_traits_to_ordered_shared_taxa(
        traits, ["A", "B", "C", "D"], ["B", "D"]
    )

    assert subset == {"B": 2.0, "D": 4.0}
    assert names == ["B", "D"]
