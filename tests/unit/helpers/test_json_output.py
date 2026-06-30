import json
import subprocess
import sys
from unittest.mock import patch

import numpy as np

import phykit.helpers.json_output as json_output_module
from phykit.helpers.json_output import print_json, to_builtin_json_types


class TestJsonOutput:
    def test_module_import_does_not_import_numpy(self):
        code = (
            "import sys; "
            "import phykit.helpers.json_output; "
            "assert 'json' not in sys.modules; "
            "assert 'numpy' not in sys.modules"
        )

        subprocess.run([sys.executable, "-c", code], check=True)

    def test_to_builtin_json_types_converts_numpy_scalars(self):
        payload = {
            "i": np.int64(7),
            "f": np.float64(1.25),
            "arr": np.array([1, 2], dtype=np.int64),
        }

        converted = to_builtin_json_types(payload)
        assert converted == {"i": 7, "f": 1.25, "arr": [1, 2]}

    def test_to_builtin_json_types_returns_builtin_scalars_directly(self):
        assert to_builtin_json_types(None) is None
        assert to_builtin_json_types(True) is True
        assert to_builtin_json_types("value") == "value"
        assert to_builtin_json_types(7) == 7
        assert to_builtin_json_types(1.25) == 1.25

    def test_to_builtin_json_types_reuses_builtin_containers(self):
        payload = {"rows": [{"a": 1, "b": [2, 3]}], "ok": True}

        converted = to_builtin_json_types(payload)

        assert converted is payload

    def test_to_builtin_json_types_reuses_unchanged_nested_branches(self):
        unchanged_row = {"a": 1, "b": [2, 3]}
        payload = {
            "rows": [unchanged_row],
            "summary": {"mean": np.float64(1.25)},
        }

        converted = to_builtin_json_types(payload)

        assert converted == {
            "rows": [{"a": 1, "b": [2, 3]}],
            "summary": {"mean": 1.25},
        }
        assert converted is not payload
        assert converted["rows"] is payload["rows"]
        assert converted["rows"][0] is unchanged_row
        assert converted["summary"] is not payload["summary"]

    def test_to_builtin_json_types_preserves_unchanged_list_prefix(self):
        unchanged_row = {"a": 1, "b": [2, 3]}
        payload = [unchanged_row, np.int64(4)]

        converted = to_builtin_json_types(payload)

        assert converted == [{"a": 1, "b": [2, 3]}, 4]
        assert converted is not payload
        assert converted[0] is unchanged_row

    @patch("builtins.print")
    def test_print_json(self, mocked_print):
        payload = {"a": np.int64(1), "b": np.float64(2.5)}
        print_json(payload)
        parsed = json.loads(mocked_print.call_args.args[0])
        assert parsed == {"a": 1, "b": 2.5}

    @patch("builtins.print")
    def test_print_json_skips_conversion_for_builtin_payload(
        self, mocked_print, monkeypatch
    ):
        def fail_conversion(_payload):
            raise AssertionError("builtin payloads should serialize directly")

        monkeypatch.setattr(
            json_output_module,
            "to_builtin_json_types",
            fail_conversion,
        )

        print_json({"rows": [{"a": 1, "b": [2, 3]}]})

        assert json.loads(mocked_print.call_args.args[0]) == {
            "rows": [{"a": 1, "b": [2, 3]}],
        }

    @patch("builtins.print")
    def test_print_json_converts_numpy_without_recursive_payload_conversion(
        self, mocked_print, monkeypatch
    ):
        def fail_conversion(_payload):
            raise AssertionError("print_json should use the JSON default hook")

        monkeypatch.setattr(
            json_output_module,
            "to_builtin_json_types",
            fail_conversion,
        )

        print_json(
            {
                "rows": [{"a": 1, "b": [2, 3]}],
                "summary": {"n": np.int64(2)},
            }
        )

        assert json.loads(mocked_print.call_args.args[0]) == {
            "rows": [{"a": 1, "b": [2, 3]}],
            "summary": {"n": 2},
        }

    def test_print_json_raises_type_error_for_unknown_objects(self):
        class Unknown:
            pass

        with patch("builtins.print"):
            try:
                print_json({"value": Unknown()})
            except TypeError as exc:
                assert "Object of type Unknown is not JSON serializable" in str(exc)
            else:
                raise AssertionError("unsupported objects should raise TypeError")

    def test_to_builtin_json_types_converts_tuple_and_nested_values(self):
        payload = {
            "t": (np.int64(1), np.float64(2.0)),
            "nested": [{"x": np.array([np.int64(3)])}],
        }
        assert to_builtin_json_types(payload) == {
            "t": [1, 2.0],
            "nested": [{"x": [3]}],
        }

    @patch("builtins.print", side_effect=BrokenPipeError)
    def test_print_json_handles_broken_pipe(self, _mocked_print):
        print_json({"a": 1})
