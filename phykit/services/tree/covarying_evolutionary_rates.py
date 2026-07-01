import math

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyPickle:
    def __getattr__(self, name):
        import pickle as _pickle

        return getattr(_pickle, name)

    def dumps(self, *args, **kwargs):
        import pickle as _pickle

        return _pickle.dumps(*args, **kwargs)

    def loads(self, *args, **kwargs):
        import pickle as _pickle

        return _pickle.loads(*args, **kwargs)


pickle = _LazyPickle()


def ProcessPoolExecutor(*args, **kwargs):
    from concurrent.futures import ProcessPoolExecutor as _ProcessPoolExecutor

    return _ProcessPoolExecutor(*args, **kwargs)


def as_completed(*args, **kwargs):
    from concurrent.futures import as_completed as _as_completed

    return _as_completed(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()
_ZSCORE_FAST_MAX_SIZE = 50_000


def _zscore(values):
    arr = np.asarray(values, dtype=float)
    if 0 < arr.size <= _ZSCORE_FAST_MAX_SIZE:
        centered = arr - arr.mean()
        variance = (centered * centered).sum() / centered.size
        return centered / math.sqrt(float(variance))
    return (arr - np.mean(arr)) / np.std(arr)


def _pearsonr(x_values, y_values):
    x = np.asarray(x_values, dtype=float)
    y = np.asarray(y_values, dtype=float)
    if x.size != y.size:
        raise ValueError("x and y must have the same length")
    if x.size < 2:
        raise ValueError("x and y must have length at least 2")

    x_centered = x - x.mean()
    y_centered = y - y.mean()
    x_ss = float(np.dot(x_centered, x_centered))
    y_ss = float(np.dot(y_centered, y_centered))
    if x_ss == 0.0 or y_ss == 0.0:
        return float("nan"), float("nan")

    r_value = float(np.dot(x_centered, y_centered) / math.sqrt(x_ss * y_ss))
    r_value = max(min(r_value, 1.0), -1.0)
    if x.size == 2:
        return r_value, 1.0

    from scipy.special import betainc

    shape = x.size / 2.0 - 1.0
    p_value = 2.0 * float(betainc(shape, shape, 0.5 * (1.0 - abs(r_value))))
    return r_value, min(max(p_value, 0.0), 1.0)


class CovaryingEvolutionaryRates(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            tree1_file_path=parsed["tree1_file_path"],
            reference=parsed["reference"],
            verbose=parsed["verbose"],
        )
        self.json_output = parsed["json_output"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def run(self):
        tree_zero = self.read_tree_file()
        tree_one = self.read_tree1_file()
        tree_ref = self.read_reference_tree_file()

        # - Calculate correlation between two gene trees
        # and save results to an array, corrArr.
        # - Branch lengths will also be part of output

        # get tree tip names
        tree_zero_tips = self.get_tip_names_from_tree(tree_zero)
        tree_one_tips = self.get_tip_names_from_tree(tree_one)
        tree_ref_tips = self.get_tip_names_from_tree(tree_ref)

        # get shared tips between the two trees
        shared_tree_tips = self.shared_tips(tree_zero_tips, tree_one_tips)

        # find differences between tree tips and shared tips
        # to determine what tips to prune
        (
            tree_zero_tips_to_prune,
            tree_one_tips_to_prune,
            tree_ref_tips_to_prune,
        ) = self._tips_to_prune_for_shared(
            tree_zero_tips, tree_one_tips, tree_ref_tips, shared_tree_tips
        )

        # get a set of pruned trees
        tree_zero = self.prune_tips(tree_zero, tree_zero_tips_to_prune)
        tree_one = self.prune_tips(tree_one, tree_one_tips_to_prune)
        tree_ref = self.prune_tips(tree_ref, tree_ref_tips_to_prune)

        # obtain corrected branch lengths where branch lengths
        # are corrected by the species tree branch length
        (
            tree_zero_corr_branch_lengths,
            tree_one_corr_branch_lengths,
            tip_names,
        ) = self.correct_branch_lengths(tree_zero, tree_one, tree_ref)

        # remove corrected BLs greater than 5
        outlier_indices = []
        outlier_indices = self.get_indices_of_outlier_branch_lengths(
            tree_zero_corr_branch_lengths, outlier_indices
        )
        outlier_indices = self.get_indices_of_outlier_branch_lengths(
            tree_one_corr_branch_lengths, outlier_indices
        )

        tree_zero_corr_branch_lengths = self.remove_outliers_based_on_indices(
            tree_zero_corr_branch_lengths, outlier_indices
        )
        tree_one_corr_branch_lengths = self.remove_outliers_based_on_indices(
            tree_one_corr_branch_lengths, outlier_indices
        )
        tip_names = self.remove_outliers_based_on_indices(tip_names, outlier_indices)

        # standardize values for final correction
        tree_zero_corr_branch_lengths = _zscore(tree_zero_corr_branch_lengths)
        tree_one_corr_branch_lengths = _zscore(tree_one_corr_branch_lengths)

        # Calculate correlation and append to results array
        # also keep a list of p values
        corr = list(
            _pearsonr(tree_zero_corr_branch_lengths, tree_one_corr_branch_lengths)
        )

        if self.plot:
            self._plot_covarying_rates_scatter(
                tree_zero_corr_branch_lengths,
                tree_one_corr_branch_lengths,
                corr,
            )

        try:
            if self.json_output:
                if self.verbose:
                    rows = [
                        dict(
                            tree_zero_rate=round(float(val_zero), 4),
                            tree_one_rate=round(float(val_one), 4),
                            branch=";".join(tip_name),
                        )
                        for val_zero, val_one, tip_name in zip(
                            tree_zero_corr_branch_lengths,
                            tree_one_corr_branch_lengths,
                            tip_names,
                        )
                    ]
                    payload = dict(verbose=True, rows=rows, branches=rows)
                    if self.plot:
                        payload["plot_output"] = self.plot_output
                    print_json(payload)
                else:
                    payload = dict(
                        verbose=False,
                        correlation=round(float(corr[0]), 4),
                        p_value=round(float(corr[1]), 6),
                    )
                    if self.plot:
                        payload["plot_output"] = self.plot_output
                    print_json(payload)
                return

            if self.verbose:
                lines = [
                    f"{round(val_zero, 4)}\t{round(val_one, 4)}\t{';'.join(tip_name)}"
                    for val_zero, val_one, tip_name in zip(
                        tree_zero_corr_branch_lengths,
                        tree_one_corr_branch_lengths,
                        tip_names,
                    )
                ]
                if lines:
                    print("\n".join(lines))
            else:
                print(f"{round(corr[0], 4)}\t{round(corr[1], 6)}")
            if self.plot:
                print(f"Saved covarying rates plot: {self.plot_output}")
        except BrokenPipeError:
            pass

    def process_args(self, args):
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree_file_path=args.tree_zero,
            tree1_file_path=args.tree_one,
            reference=args.reference,
            verbose=args.verbose,
            json_output=getattr(args, "json", False),
            plot=getattr(args, "plot", False),
            plot_output=getattr(args, "plot_output", "covarying_rates_plot.png"),
            plot_config=PlotConfig.from_args(args),
        )

    @staticmethod
    def _tips_to_prune_for_shared(
        tree_zero_tips, tree_one_tips, tree_ref_tips, shared_tree_tips
    ):
        shared_tip_set = set(shared_tree_tips)
        return (
            [tip for tip in tree_zero_tips if tip not in shared_tip_set],
            [tip for tip in tree_one_tips if tip not in shared_tip_set],
            [tip for tip in tree_ref_tips if tip not in shared_tip_set],
        )

    def _plot_covarying_rates_scatter(self, x_vals, y_vals, corr):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print(
                "matplotlib is required for --plot in covarying_evolutionary_rates. Install matplotlib and retry."
            )
            raise SystemExit(2)

        x = np.asarray(x_vals, dtype=float)
        y = np.asarray(y_vals, dtype=float)
        finite_mask = np.isfinite(x) & np.isfinite(y)
        x = x[finite_mask]
        y = y[finite_mask]

        config = self.plot_config
        config.resolve(n_rows=len(x), n_cols=None)
        default_colors = ["#2b8cbe", "#000000"]
        colors = config.merge_colors(default_colors)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        ax.scatter(
            x,
            y,
            s=14,
            alpha=0.6,
            color=colors[0],
            edgecolors="none",
        )

        if x.size >= 2:
            slope, intercept = np.polyfit(x, y, 1)
            x_line = np.linspace(float(np.min(x)), float(np.max(x)), 200)
            y_line = slope * x_line + intercept
            ax.plot(
                x_line,
                y_line,
                color=colors[1],
                linestyle="--",
                linewidth=2.0,
                label=f"Best fit (slope={slope:.4f})",
            )
            legend_loc = config.legend_position or "best"
            if legend_loc != "none":
                ax.legend(loc=legend_loc, frameon=False)

        if config.show_title:
            ax.set_title(config.title or "Covarying Evolutionary Rates", fontsize=config.title_fontsize)
        ax.set_xlabel("Tree zero relative rate (z-score)")
        ax.set_ylabel("Tree one relative rate (z-score)")
        ax.text(
            0.02,
            0.98,
            f"r={float(corr[0]):.4f}\np={float(corr[1]):.6g}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.8, edgecolor="none"),
        )
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)
        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def get_indices_of_outlier_branch_lengths(
        self, corr_branch_lengths, outlier_indices
    ):
        """
        create index for branch lengths that
        have an absolute value greater than 5
        """
        # Convert to numpy array for vectorized operations
        arr = np.array(corr_branch_lengths, dtype=float)

        # Find outliers using vectorized operations
        new_outliers = np.flatnonzero((np.abs(arr) > 5) | np.isnan(arr))
        if not outlier_indices:
            return new_outliers.tolist()

        # Combine with existing outliers
        all_outliers = set(outlier_indices)
        all_outliers.update(new_outliers.tolist())

        return list(all_outliers)

    def remove_outliers_based_on_indices(self, corr_branch_lengths, outlier_indices):
        """
        remove value if the value is an outlier according
        to the outlier indices list
        """
        if not outlier_indices:
            return corr_branch_lengths

        value_type = corr_branch_lengths.__class__
        if value_type.__module__ == "numpy" and value_type.__name__ == "ndarray":
            mask = np.ones(corr_branch_lengths.shape[0], dtype=bool)
            mask[list(outlier_indices)] = False
            return corr_branch_lengths[mask].tolist()

        length = len(corr_branch_lengths)
        normalized_indices = []
        append = normalized_indices.append
        for index in outlier_indices:
            normalized_index = index if index >= 0 else length + index
            if normalized_index < 0 or normalized_index >= length:
                raise IndexError("outlier index out of range")
            append(normalized_index)

        result = []
        extend = result.extend
        start = 0
        for index in sorted(set(normalized_indices)):
            extend(corr_branch_lengths[start:index])
            start = index + 1
        extend(corr_branch_lengths[start:])
        return result

    def prune_tips(self, tree, tips):
        """
        prune tips from trees
        """

        for tip in tips:
            tree.prune(tip)

        return tree

    @staticmethod
    def _process_terminal_batch(tree0_pickle, tree1_pickle, terminals_data):
        """Process a batch of terminals in parallel."""
        t0 = pickle.loads(tree0_pickle)
        t1 = pickle.loads(tree1_pickle)

        results = []
        for terminal_name, terminal_bl, sp_tips in terminals_data:
            try:
                newtree = t0.common_ancestor(terminal_name)
                newtree1 = t1.common_ancestor(terminal_name)

                bl0 = round(newtree.branch_length / terminal_bl, 6) if newtree.branch_length else None
                bl1 = round(newtree1.branch_length / terminal_bl, 6) if newtree1.branch_length else None

                if bl0 is not None and bl1 is not None:
                    results.append((bl0, bl1, sp_tips))
            except Exception as err:
                if isinstance(err, (SystemExit, KeyboardInterrupt)):
                    raise
                continue
        return results

    @staticmethod
    def _process_nonterminal_batch(tree0_pickle, tree1_pickle, nonterminals_data):
        """Process a batch of nonterminals in parallel."""
        t0 = pickle.loads(tree0_pickle)
        t1 = pickle.loads(tree1_pickle)

        results = []
        for sp_tips, nonterminal_bl in nonterminals_data:
            try:
                newtree = t0.common_ancestor(sp_tips)
                newtree1 = t1.common_ancestor(sp_tips)

                if newtree.branch_length and newtree1.branch_length and nonterminal_bl:
                    bl0 = round(newtree.branch_length / nonterminal_bl, 6)
                    bl1 = round(newtree1.branch_length / nonterminal_bl, 6)
                    results.append((bl0, bl1, sp_tips))
            except Exception as err:
                if isinstance(err, (SystemExit, KeyboardInterrupt)):
                    raise
                continue
        return results

    def _process_terminals_sequential(self, terminals, t0, t1, l0, l1, tip_names):
        for i in terminals:
            sp_tips = self.get_tip_names_from_tree(i)
            tip_names.append(sp_tips)
            try:
                newtree = t0.common_ancestor(i.name)
                newtree1 = t1.common_ancestor(i.name)
                if newtree.branch_length and i.branch_length:
                    l0.append(round(newtree.branch_length / i.branch_length, 6))
                    l1.append(round(newtree1.branch_length / i.branch_length, 6))
            except Exception as err:
                if isinstance(err, (SystemExit, KeyboardInterrupt)):
                    raise
                continue

    def _process_nonterminals_sequential(self, nonterminals, t0, t1, l0, l1, tip_names):
        for i in nonterminals:
            sp_tips = self.get_tip_names_from_tree(i)
            try:
                newtree = t0.common_ancestor(sp_tips)
                newtree1 = t1.common_ancestor(sp_tips)
                if newtree.branch_length and newtree1.branch_length and i.branch_length:
                    l0.append(round(newtree.branch_length / i.branch_length, 6))
                    l1.append(round(newtree1.branch_length / i.branch_length, 6))
                    tip_names.append(sp_tips)
            except Exception as err:
                if isinstance(err, (SystemExit, KeyboardInterrupt)):
                    raise
                continue

    @staticmethod
    def _branch_lengths_by_tipset(tree):
        direct_result = CovaryingEvolutionaryRates._branch_lengths_by_tipset_direct(
            tree
        )
        if direct_result is not None:
            return direct_result

        clade_tips = {}
        branch_lengths = {}
        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                tips = frozenset([clade.name])
            else:
                tips = frozenset().union(
                    *(clade_tips[id(child)] for child in clade.clades)
                )
            clade_tips[id(clade)] = tips
            branch_lengths[tips] = clade.branch_length
        return branch_lengths

    @staticmethod
    def _branch_lengths_by_tipset_direct(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        clade_tips = {}
        branch_lengths = {}
        clades = []
        stack = [root]
        try:
            while stack:
                clade = stack.pop()
                clades.append(clade)
                children = clade.clades
                if children:
                    stack.extend(children)

            for clade in reversed(clades):
                children = clade.clades
                if children:
                    if len(children) == 2:
                        tips = (
                            clade_tips[id(children[0])]
                            | clade_tips[id(children[1])]
                        )
                    else:
                        tips = frozenset().union(
                            *(clade_tips[id(child)] for child in children)
                        )
                else:
                    tips = frozenset([clade.name])
                clade_tips[id(clade)] = tips
                branch_lengths[tips] = clade.branch_length
        except (AttributeError, KeyError, TypeError):
            return None
        return branch_lengths

    @staticmethod
    def _reference_tipsets_and_order_direct(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        sp_tips_by_id = {}
        sp_tip_names_by_id = {}
        try:
            preorder = []
            terminals = []
            nonterminals = []
            stack = [root]
            pop = stack.pop
            append = stack.append
            preorder_append = preorder.append
            terminals_append = terminals.append
            nonterminals_append = nonterminals.append
            while stack:
                clade = pop()
                preorder_append(clade)
                children = clade.clades
                if children:
                    nonterminals_append(clade)
                    child_count = len(children)
                    if child_count == 2:
                        append(children[1])
                        append(children[0])
                    else:
                        for index in range(child_count - 1, -1, -1):
                            append(children[index])
                else:
                    terminals_append(clade)

            for clade in reversed(preorder):
                children = clade.clades
                if children:
                    if len(children) == 2:
                        tips = (
                            sp_tips_by_id[id(children[0])]
                            | sp_tips_by_id[id(children[1])]
                        )
                    else:
                        tips = frozenset().union(
                            *(sp_tips_by_id[id(child)] for child in children)
                        )
                    tip_names = tuple(
                        name
                        for child in children
                        for name in sp_tip_names_by_id[id(child)]
                    )
                else:
                    tips = frozenset([clade.name])
                    tip_names = (clade.name,)
                sp_tips_by_id[id(clade)] = tips
                sp_tip_names_by_id[id(clade)] = tip_names
        except (AttributeError, KeyError, TypeError):
            return None

        return sp_tips_by_id, sp_tip_names_by_id, terminals, nonterminals

    def _correct_branch_lengths_from_same_tree(self, tree):
        reference_data = self._reference_tipsets_and_order_direct(tree)
        if reference_data is None:
            return None

        _, sp_tip_names_by_id, terminals, nonterminals = reference_data
        l0 = []
        l1 = []
        tip_names = []
        l0_append = l0.append
        l1_append = l1.append
        tip_names_append = tip_names.append

        for clade_group in (terminals, nonterminals):
            for clade in clade_group:
                if clade.branch_length:
                    l0_append(1.0)
                    l1_append(1.0)
                    tip_names_append(list(sp_tip_names_by_id[id(clade)]))

        return l0, l1, tip_names

    def _correct_branch_lengths_from_exact_splits(self, t0, t1, sp):
        if t0 is t1 and t0 is sp:
            same_tree_result = self._correct_branch_lengths_from_same_tree(sp)
            if same_tree_result is not None:
                return same_tree_result

        try:
            t0_lengths = self._branch_lengths_by_tipset(t0)
            t1_lengths = t0_lengths if t0 is t1 else self._branch_lengths_by_tipset(t1)
            reference_data = self._reference_tipsets_and_order_direct(sp)
            if reference_data is None:
                sp_tips_by_id = {}
                sp_tip_names_by_id = {}
                for clade in sp.find_clades(order="postorder"):
                    if clade.is_terminal():
                        tips = frozenset([clade.name])
                        tip_names = (clade.name,)
                    else:
                        tips = frozenset().union(
                            *(sp_tips_by_id[id(child)] for child in clade.clades)
                        )
                        tip_names = tuple(
                            name
                            for child in clade.clades
                            for name in sp_tip_names_by_id[id(child)]
                        )
                    sp_tips_by_id[id(clade)] = tips
                    sp_tip_names_by_id[id(clade)] = tip_names

                terminals = []
                nonterminals = []
                stack = [sp.root]
                while stack:
                    clade = stack.pop()
                    children = clade.clades
                    if children:
                        nonterminals.append(clade)
                        stack.extend(reversed(children))
                    else:
                        terminals.append(clade)
            else:
                sp_tips_by_id, sp_tip_names_by_id, terminals, nonterminals = (
                    reference_data
                )
        except (AttributeError, TypeError, KeyError):
            return None

        l0 = []
        l1 = []
        tip_names = []

        for clade_group in (terminals, nonterminals):
            for clade in clade_group:
                tips = sp_tips_by_id.get(id(clade))
                if not tips:
                    return None
                if tips not in t0_lengths or tips not in t1_lengths:
                    return None

                sp_bl = clade.branch_length
                bl0 = t0_lengths[tips]
                bl1 = t1_lengths[tips]
                if bl0 and bl1 and sp_bl:
                    l0.append(round(bl0 / sp_bl, 6))
                    l1.append(round(bl1 / sp_bl, 6))
                    tip_names.append(list(sp_tip_names_by_id[id(clade)]))

        return l0, l1, tip_names

    def correct_branch_lengths(self, t0, t1, sp):
        """
        obtain a list of corrected branch lengths with parallel processing
        """
        fast_result = self._correct_branch_lengths_from_exact_splits(t0, t1, sp)
        if fast_result is not None:
            return fast_result

        l0 = []
        l1 = []
        tip_names = []

        # Collect terminal data
        terminals = sp.get_terminals()
        nonterminals = sp.get_nonterminals()

        # Process sequentially if small dataset or use parallel processing
        if len(terminals) + len(nonterminals) < 50:
            # Original sequential processing for small datasets
            self._process_terminals_sequential(terminals, t0, t1, l0, l1, tip_names)
            self._process_nonterminals_sequential(nonterminals, t0, t1, l0, l1, tip_names)
        else:
            # Parallel processing for large datasets
            tree0_pickle = pickle.dumps(t0)
            tree1_pickle = pickle.dumps(t1)

            # Prepare terminal data
            terminals_data = []
            for i in terminals:
                sp_tips = self.get_tip_names_from_tree(i)
                if i.branch_length:
                    terminals_data.append((i.name, i.branch_length, sp_tips))

            # Prepare nonterminal data
            nonterminals_data = []
            for i in nonterminals:
                if i.branch_length:
                    sp_tips = self.get_tip_names_from_tree(i)
                    nonterminals_data.append((sp_tips, i.branch_length))

            # Process in batches
            batch_size = max(10, (len(terminals_data) + len(nonterminals_data)) // 4)

            try:
                with ProcessPoolExecutor(max_workers=min(4, len(terminals_data) + len(nonterminals_data) // 10)) as executor:
                    futures = []

                    # Submit terminal batches
                    for i in range(0, len(terminals_data), batch_size):
                        batch = terminals_data[i:i+batch_size]
                        futures.append(executor.submit(self._process_terminal_batch, tree0_pickle, tree1_pickle, batch))

                    # Submit nonterminal batches
                    for i in range(0, len(nonterminals_data), batch_size):
                        batch = nonterminals_data[i:i+batch_size]
                        futures.append(executor.submit(self._process_nonterminal_batch, tree0_pickle, tree1_pickle, batch))

                    for future in as_completed(futures):
                        batch_results = future.result()
                        for bl0, bl1, sp_tips in batch_results:
                            l0.append(bl0)
                            l1.append(bl1)
                            tip_names.append(sp_tips)
            except (OSError, ValueError, RuntimeError):
                self._process_terminals_sequential(terminals, t0, t1, l0, l1, tip_names)
                self._process_nonterminals_sequential(nonterminals, t0, t1, l0, l1, tip_names)

        return (l0, l1, tip_names)
