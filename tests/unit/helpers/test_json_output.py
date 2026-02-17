import json
from unittest.mock import patch

import numpy as np

from phykit.helpers.json_output import print_json, to_builtin_json_types


class TestJsonOutput:
    def test_to_builtin_json_types_converts_numpy_scalars(self):
        payload = {
            "i": np.int64(7),
            "f": np.float64(1.25),
            "arr": np.array([1, 2], dtype=np.int64),
        }

        converted = to_builtin_json_types(payload)
        assert converted == {"i": 7, "f": 1.25, "arr": [1, 2]}

    @patch("builtins.print")
    def test_print_json(self, mocked_print):
        payload = {"a": np.int64(1), "b": np.float64(2.5)}
        print_json(payload)
        parsed = json.loads(mocked_print.call_args.args[0])
        assert parsed == {"a": 1, "b": 2.5}
