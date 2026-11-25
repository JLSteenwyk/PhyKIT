import subprocess
import sys


def test_service_numpy_loaders_keep_separate_caches():
    subprocess.run(
        [
            sys.executable,
            "-c",
            """
import sys

from phykit.services.alignment import alignment_entropy
from phykit.services.tree import phylogenetic_signal

alignment_numpy = alignment_entropy.np
tree_numpy = phylogenetic_signal.np
assert "numpy" not in sys.modules
assert alignment_numpy is not tree_numpy

zeros = alignment_numpy.zeros
assert alignment_numpy.__dict__["zeros"] is zeros
assert tree_numpy._module is None
assert "zeros" not in vars(tree_numpy)
assert tree_numpy.zeros is zeros
""",
        ],
        check=True,
    )
