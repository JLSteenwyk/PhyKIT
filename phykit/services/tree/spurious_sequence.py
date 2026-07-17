from __future__ import annotations

from math import erfc, isfinite, log, sqrt

from .base import Tree
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyNumpy:
    _module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy as _np

            module = _np
            self._module = module

        value = getattr(module, name)
        setattr(self, name, value)
        return value


np = _LazyNumpy()

_MEDIAN_NUMPY_THRESHOLD = 1024
_DEFAULT_ALPHA = 0.05
_MIN_KDE_OBSERVATIONS = 4
_SQRT_TWO = sqrt(2.0)


class SpuriousSequence(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            factor=parsed["factor"],
        )
        self.json_output = parsed["json_output"]
        self.method = parsed["method"]
        self.alpha = parsed["alpha"]
        self.max_remove = parsed["max_remove"]
        self.tree_list = parsed["tree_list"]
        self.per_species = parsed["per_species"]

    def run(self) -> None:
        if self.method == "diameter-impact":
            self._run_diameter_impact()
            return

        self._run_median_factor()

    def _run_median_factor(self) -> None:
        tree = self.read_tree_file_unmodified()
        name_and_branch_len, threshold, median = \
            self.identify_spurious_sequence(
                tree, self.factor
            )

        rounded_threshold = round(threshold, 4)
        rounded_median = round(median, 4)
        has_spurious_sequence = self._has_spurious_sequence(
            name_and_branch_len,
            threshold,
        )

        if self.json_output:
            if not has_spurious_sequence:
                print_json(dict(rows=[], spurious_sequences=[]))
                return

            spurious_rows = [
                {
                    "taxon": name,
                    "branch_length": round(length, 4),
                    "threshold": rounded_threshold,
                    "median": rounded_median,
                }
                for name, length in name_and_branch_len.items()
                if length >= threshold
            ]
            print_json(dict(rows=spurious_rows, spurious_sequences=spurious_rows))
            return

        if not has_spurious_sequence:
            print("None")
            return

        try:
            print(
                "\n".join(
                    f"{name}\t{round(length, 4)}\t{rounded_threshold}\t{rounded_median}"
                    for name, length in name_and_branch_len.items()
                    if length >= threshold
                )
            )
        except BrokenPipeError:
            pass

    def process_args(self, args) -> dict:
        method = getattr(args, "method", "median-factor") or "median-factor"
        factor_arg = getattr(args, "factor", None)
        alpha = getattr(args, "alpha", _DEFAULT_ALPHA)
        max_remove = getattr(args, "max_remove", None)
        tree_list = getattr(args, "tree_list", False)
        per_species = getattr(args, "per_species", False)

        if method not in {"median-factor", "diameter-impact"}:
            raise PhykitUserError(
                ["--method must be median-factor or diameter-impact."],
                code=2,
            )
        if method == "median-factor":
            if (
                tree_list
                or per_species
                or max_remove is not None
                or alpha != _DEFAULT_ALPHA
            ):
                raise PhykitUserError(
                    [
                        "--alpha, --tree-list, --per-species, and --max-remove "
                        "require --method diameter-impact."
                    ],
                    code=2,
                )
        else:
            if factor_arg is not None:
                raise PhykitUserError(
                    ["--factor cannot be used with --method diameter-impact."],
                    code=2,
                )
            if not isfinite(alpha) or not 0.0 < alpha < 1.0:
                raise PhykitUserError(
                    ["--alpha must be a finite value between 0 and 1."],
                    code=2,
                )
            if max_remove is not None and max_remove < 1:
                raise PhykitUserError(
                    ["--max-remove must be at least 1."],
                    code=2,
                )
            if per_species and not tree_list:
                raise PhykitUserError(
                    ["--per-species requires --tree-list."],
                    code=2,
                )

        return dict(
            tree_file_path=args.tree,
            factor=factor_arg or 20,
            json_output=getattr(args, "json", False),
            method=method,
            alpha=alpha,
            max_remove=max_remove,
            tree_list=tree_list,
            per_species=per_species,
        )

    def _run_diameter_impact(self) -> None:
        tree_paths = (
            self._read_tree_list(self.tree_file_path)
            if self.tree_list
            else [self.tree_file_path]
        )
        analyses = [self._analyze_tree_path(path) for path in tree_paths]

        if self.tree_list:
            scope = "per-species" if self.per_species else "all-genes"
            self._assign_collection_p_values(analyses, per_species=self.per_species)
        else:
            scope = "per-gene"
            self._assign_per_gene_p_values(analyses[0])

        rows = self._diameter_impact_rows(analyses)
        if self.json_output:
            print_json(
                {
                    "method": self.method,
                    "scope": scope,
                    "alpha": self.alpha,
                    "tree_count": len(analyses),
                    "rows": rows,
                    "spurious_sequences": rows,
                }
            )
            return

        if not rows:
            print("None")
            return

        include_tree = self.tree_list
        try:
            print(
                "\n".join(
                    self._format_diameter_impact_row(row, include_tree)
                    for row in rows
                )
            )
        except BrokenPipeError:
            pass

    @staticmethod
    def _read_tree_list(list_path: str) -> list[str]:
        from pathlib import Path

        source = Path(list_path)
        try:
            handle = source.open()
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{list_path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            ) from None

        paths = []
        with handle:
            for line in handle:
                value = line.strip()
                if not value or value.startswith("#"):
                    continue
                path = Path(value)
                if not path.is_absolute():
                    path = source.parent / path
                paths.append(str(path))

        if not paths:
            raise PhykitUserError(
                [f"Tree list {list_path} contains no tree paths."],
                code=2,
            )
        return paths

    def _analyze_tree_path(self, tree_path: str) -> dict:
        original_path = self.tree_file_path
        self.tree_file_path = tree_path
        try:
            tree = self.read_tree_file_unmodified()
        finally:
            self.tree_file_path = original_path

        self.validate_tree(
            tree,
            min_tips=4,
            require_branch_lengths=True,
            context="diameter-impact detection",
        )
        adjacency, terminal_names, terminal_lengths = self._tree_graph(tree)
        tip_count = len(terminal_names)
        max_remove = (
            self.max_remove
            if self.max_remove is not None
            else self._default_max_remove(tip_count)
        )
        if max_remove > tip_count - 2:
            raise PhykitUserError(
                [
                    f"--max-remove must be at most {tip_count - 2} for tree "
                    f"{tree_path}."
                ],
                code=2,
            )

        trajectory, initial_diameter = self._diameter_shrink_trajectory(
            adjacency,
            terminal_names,
            max_remove,
        )
        scores = {name: 0.0 for name in terminal_names.values()}
        details = {}
        for step in trajectory:
            taxon = step["taxon"]
            scores[taxon] = step["signature"]
            step["terminal_branch_length"] = terminal_lengths[taxon]
            details[taxon] = step

        return {
            "tree": tree_path,
            "tip_count": tip_count,
            "max_remove": max_remove,
            "initial_diameter": initial_diameter,
            "scores": scores,
            "details": details,
            "p_values": {},
        }

    @staticmethod
    def _default_max_remove(tip_count: int) -> int:
        return max(
            1,
            min(tip_count // 4, int(5 * sqrt(tip_count)), tip_count - 2),
        )

    @staticmethod
    def _tree_graph(
        tree,
    ) -> tuple[
        list[list[tuple[int, float]]],
        dict[int, str],
        dict[str, float],
    ]:
        adjacency = [[]]
        terminal_names = {}
        terminal_lengths = {}
        stack = [(tree.root, 0)]
        while stack:
            clade, node_index = stack.pop()
            children = clade.clades
            if not children:
                terminal_names[node_index] = clade.name
                terminal_lengths[clade.name] = float(clade.branch_length)
                continue

            for child in reversed(children):
                child_index = len(adjacency)
                adjacency.append([])
                branch_length = float(child.branch_length)
                adjacency[node_index].append((child_index, branch_length))
                adjacency[child_index].append((node_index, branch_length))
                stack.append((child, child_index))
        return adjacency, terminal_names, terminal_lengths

    @classmethod
    def _diameter_shrink_trajectory(
        cls,
        adjacency: list[list[tuple[int, float]]],
        terminal_names: dict[int, str],
        max_remove: int,
    ) -> tuple[list[dict], float]:
        active = set(terminal_names)
        current_diameter, endpoint_a, endpoint_b = cls._active_diameter(
            adjacency,
            active,
        )
        initial_diameter = current_diameter
        trajectory = []

        for rank in range(1, max_remove + 1):
            if current_diameter <= 0.0:
                break

            candidates = []
            for endpoint in (endpoint_a, endpoint_b):
                active.remove(endpoint)
                diameter, next_a, next_b = cls._active_diameter(adjacency, active)
                active.add(endpoint)
                candidates.append(
                    (
                        diameter,
                        terminal_names[endpoint],
                        endpoint,
                        next_a,
                        next_b,
                    )
                )

            next_diameter, taxon, removed, next_a, next_b = min(candidates)
            if next_diameter <= 0.0:
                signature = -log(1e-15)
            else:
                signature = max(0.0, log(current_diameter / next_diameter))

            trajectory.append(
                {
                    "taxon": taxon,
                    "signature": signature,
                    "removal_rank": rank,
                    "diameter_before": current_diameter,
                    "diameter_after": next_diameter,
                }
            )
            active.remove(removed)
            current_diameter = next_diameter
            endpoint_a, endpoint_b = next_a, next_b

        return trajectory, initial_diameter

    @classmethod
    def _active_diameter(
        cls,
        adjacency: list[list[tuple[int, float]]],
        active: set[int],
    ) -> tuple[float, int, int]:
        start = min(active)
        endpoint_a, _distance = cls._farthest_active(adjacency, active, start)
        endpoint_b, diameter = cls._farthest_active(
            adjacency,
            active,
            endpoint_a,
        )
        return diameter, endpoint_a, endpoint_b

    @staticmethod
    def _farthest_active(
        adjacency: list[list[tuple[int, float]]],
        active: set[int],
        start: int,
    ) -> tuple[int, float]:
        farthest = start
        farthest_distance = 0.0
        stack = [(start, -1, 0.0)]
        while stack:
            node, parent, distance = stack.pop()
            if node in active and (
                distance > farthest_distance
                or (distance == farthest_distance and node < farthest)
            ):
                farthest = node
                farthest_distance = distance
            for neighbor, branch_length in adjacency[node]:
                if neighbor != parent:
                    stack.append((neighbor, node, distance + branch_length))
        return farthest, farthest_distance

    def _assign_per_gene_p_values(self, analysis: dict) -> None:
        values = list(analysis["scores"].values())
        mean = sum(values) / len(values)
        variance = sum((value - mean) ** 2 for value in values) / len(values)
        if variance <= 0.0:
            analysis["p_values"] = {name: 1.0 for name in analysis["scores"]}
            return

        standard_deviation = sqrt(variance)
        analysis["p_values"] = {
            name: self._normal_survival(value, mean, standard_deviation)
            for name, value in analysis["scores"].items()
        }

    def _assign_collection_p_values(
        self,
        analyses: list[dict],
        per_species: bool,
    ) -> None:
        observations = {}
        for analysis_index, analysis in enumerate(analyses):
            for taxon, score in analysis["scores"].items():
                key = taxon if per_species else "all"
                observations.setdefault(key, []).append(
                    (analysis_index, taxon, score)
                )

        for group in observations.values():
            values = [observation[2] for observation in group]
            probabilities = self._kde_survival_probabilities(values)
            for observation, probability in zip(group, probabilities):
                analysis_index, taxon, _score = observation
                analyses[analysis_index]["p_values"][taxon] = probability

    @classmethod
    def _kde_survival_probabilities(cls, values: list[float]) -> list[float]:
        count = len(values)
        if count < _MIN_KDE_OBSERVATIONS:
            return [1.0] * count

        mean = sum(values) / count
        variance = sum((value - mean) ** 2 for value in values) / (count - 1)
        if variance <= 0.0:
            return [1.0] * count

        bandwidth = 1.06 * sqrt(variance) * count ** -0.2
        if bandwidth <= 0.0 or not isfinite(bandwidth):
            return [1.0] * count

        value_counts = {}
        for value in values:
            value_counts[value] = value_counts.get(value, 0) + 1

        probability_by_value = {}
        denominator = count - 1
        for value in value_counts:
            cumulative = 0.0
            for other, frequency in value_counts.items():
                cumulative += frequency * cls._normal_cdf(
                    (value - other) / bandwidth
                )
            cumulative -= 0.5
            probability = 1.0 - cumulative / denominator
            probability_by_value[value] = min(1.0, max(0.0, probability))
        return [probability_by_value[value] for value in values]

    @staticmethod
    def _normal_cdf(value: float) -> float:
        return 0.5 * erfc(-value / _SQRT_TWO)

    @staticmethod
    def _normal_survival(
        value: float,
        mean: float,
        standard_deviation: float,
    ) -> float:
        return 0.5 * erfc(
            (value - mean) / (standard_deviation * _SQRT_TWO)
        )

    def _diameter_impact_rows(self, analyses: list[dict]) -> list[dict]:
        rows = []
        for analysis in analyses:
            ordered_details = sorted(
                analysis["details"].values(),
                key=lambda detail: detail["removal_rank"],
            )
            for detail in ordered_details:
                taxon = detail["taxon"]
                p_value = analysis["p_values"].get(taxon, 1.0)
                if detail["signature"] <= 0.0 or p_value >= self.alpha:
                    continue
                rows.append(
                    {
                        "tree": analysis["tree"],
                        "taxon": taxon,
                        "signature": round(detail["signature"], 6),
                        "p_value": round(p_value, 6),
                        "removal_rank": detail["removal_rank"],
                        "terminal_branch_length": round(
                            detail["terminal_branch_length"],
                            6,
                        ),
                        "diameter_before": round(detail["diameter_before"], 6),
                        "diameter_after": round(detail["diameter_after"], 6),
                    }
                )
        return rows

    @staticmethod
    def _format_diameter_impact_row(row: dict, include_tree: bool) -> str:
        fields = []
        if include_tree:
            fields.append(row["tree"])
        fields.extend(
            [
                row["taxon"],
                str(row["signature"]),
                str(row["p_value"]),
                str(row["removal_rank"]),
                str(row["terminal_branch_length"]),
                str(row["diameter_before"]),
                str(row["diameter_after"]),
            ]
        )
        return "\t".join(fields)

    def identify_spurious_sequence(
        self,
        tree: Newick.Tree,
        factor: float,
    ) -> tuple[
        dict[str, float],
        float,
        float,
    ]:
        branch_lengths, name_and_branch_len = \
            self.get_branch_lengths_and_their_names(tree)

        median = self._median_branch_length(branch_lengths)

        threshold = median * factor

        return name_and_branch_len, threshold, median

    @staticmethod
    def _has_spurious_sequence(name_and_branch_len, threshold) -> bool:
        if not name_and_branch_len:
            return False

        longest_branch_length = max(name_and_branch_len.values())
        if longest_branch_length >= threshold:
            return True

        if longest_branch_length != longest_branch_length:
            for length in name_and_branch_len.values():
                if length >= threshold:
                    return True

        return False

    @staticmethod
    def _median_branch_length(branch_lengths: list[float]) -> float:
        count = len(branch_lengths)
        if count == 0:
            raise ValueError("no median for empty data")

        middle = count // 2
        if count >= _MEDIAN_NUMPY_THRESHOLD:
            values = np.asarray(branch_lengths, dtype=float)
            if count % 2:
                values.partition(middle)
                return float(values[middle])
            values.partition((middle - 1, middle))
            return float((values[middle - 1] + values[middle]) / 2)

        values = sorted(branch_lengths)
        if count % 2:
            return values[middle]
        return (values[middle - 1] + values[middle]) / 2

    def get_branch_lengths_and_their_names(
        self,
        tree: Newick.Tree,
    ) -> tuple[
        list[float],
        dict[str, float],
    ]:
        direct_result = self._get_branch_lengths_and_names_standard_tree(tree)
        if direct_result is not None:
            return direct_result

        branch_lengths = []
        name_and_branch_len = {}

        # collect terminal branch lengths only for spurious sequence detection
        # (internal branches are not considered for spurious sequence detection)
        terminals = self._iter_terminal_clades(tree)
        for terminal in terminals:
            if terminal.branch_length is not None:
                branch_lengths.append(terminal.branch_length)
                name_and_branch_len[terminal.name] = terminal.branch_length

        return branch_lengths, name_and_branch_len

    @staticmethod
    def _get_branch_lengths_and_names_standard_tree(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        branch_lengths = []
        name_and_branch_len = {}
        stack = [root]
        pop = stack.pop
        append = stack.append
        append_length = branch_lengths.append
        try:
            while stack:
                clade = pop()
                children = clade.clades
                child_count = len(children)
                if child_count == 2:
                    append(children[1])
                    append(children[0])
                elif child_count:
                    for idx in range(child_count - 1, -1, -1):
                        append(children[idx])
                else:
                    branch_length = clade.branch_length
                    if branch_length is not None:
                        append_length(branch_length)
                        name_and_branch_len[clade.name] = branch_length
        except AttributeError:
            return None

        return branch_lengths, name_and_branch_len

    @staticmethod
    def _iter_terminal_clades(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return tree.get_terminals()

        terminals = []
        stack = [root]
        pop = stack.pop
        append = stack.append
        append_terminal = terminals.append
        while stack:
            clade = pop()
            children = clade.clades
            child_count = len(children)
            if child_count == 2:
                append(children[1])
                append(children[0])
            elif child_count:
                for idx in range(child_count - 1, -1, -1):
                    append(children[idx])
            else:
                append_terminal(clade)
        return terminals
