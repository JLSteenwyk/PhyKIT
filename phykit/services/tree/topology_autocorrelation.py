"""Genomic distance-decay analysis of strict focal topology classifications."""

import csv
import json
import math
import sys
import textwrap
from pathlib import Path

import numpy as np

from ...helpers.genomic_coordinates import select_rows
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig
from ...helpers.topology_autocorrelation import (
    METRICS,
    analyze,
    read_classified,
    uncertainty,
    validate_edges,
)
from ...helpers.topology_landscape import CLASSES, fail
from .topology_landscape import COLORS, TopologyLandscape


class TopologyAutocorrelation:
    def __init__(self, args):
        self.args = args

    def prepare(self):
        args = self.args
        if args.classified_genes:
            raw_options = (args.manifest, args.coordinates, args.groups, args.reference_tree,
                           args.branch_taxa, args.outgroup, args.min_support, args.support_scale,
                           args.missing_support, args.unmapped)
            if any(value is not None for value in raw_options):
                fail("--classified-genes cannot be combined with raw-tree or support-filtering options.")
            rows = select_rows(read_classified(args.classified_genes), args.chromosome, args.interval)
            if not rows:
                fail("No classified genes remain in the selected chromosomes/interval.")
            labels = dict(zip(CLASSES[:3], args.labels or CLASSES[:3]))
            return {"genes": rows, "labels": labels, "classification_source": args.classified_genes,
                    "diagnostics": {"classification_policy": "retained_from_input_not_recomputed"}}, {
                        Path(args.classified_genes).resolve()}
        if not args.manifest or not args.coordinates:
            fail("Raw input requires --manifest and --coordinates.")
        args.support_scale = args.support_scale if args.support_scale is not None else 100
        args.missing_support = args.missing_support or "collapse"
        args.unmapped = args.unmapped or "error"
        return TopologyLandscape(args).prepare()

    def run(self):
        args = self.args
        if args.distance_edges is not None:
            if args.bin_width is not None or args.max_distance is not None:
                fail("Use --distance-edges or --bin-width/--max-distance, not both.")
            edges = validate_edges(args.distance_edges).tolist()
        else:
            width = args.bin_width if args.bin_width is not None else 10000
            maximum = args.max_distance if args.max_distance is not None else 1000000
            if width <= 0 or maximum <= 0 or maximum > 2**62:
                fail("Distance width and maximum must be positive; maximum must be <= 2**62.")
            if (maximum + width - 1) // width > 1000:
                fail("Use at most 1000 distance bins.")
            edges = [*range(0, maximum, width), maximum]
        if args.seed < 0:
            fail("--seed must be nonnegative.")
        source, inputs = self.prepare()
        rows, labels = source["genes"], source["labels"]
        result, _ = analyze(rows, edges)
        intervals, blocks = uncertainty(rows, edges, args.block_sizes, args.replicates, args.seed) \
            if args.block_sizes else ([], [])
        warnings = [
            "Chromosome-wide baselines do not control spatially varying topology frequencies.",
            "Unresolved genes are excluded from pair denominators; nonrandom missingness can bias results.",
            "Clustering does not establish introgression, recombination breakpoints, or independent-locus counts.",
        ]
        if intervals:
            warnings.append("Intervals are pointwise, not simultaneous bands or random-label null intervals.")
            if any(r["status"] == "withheld" for r in intervals):
                warnings.append("Some intervals were withheld; inspect interval reasons and block diagnostics.")
            if any(r["block_size_sensitive"] for r in intervals):
                warnings.append("Interval widths are sensitive to block size.")
        else:
            warnings.append("No block sizes supplied: descriptive estimates only, without uncertainty intervals.")
        payload = {
            "command": "topology_autocorrelation", "coordinate_system": "0-based-half-open",
            "parameters": {"distance_edges": edges, "maximum_distance_exclusive": edges[-1],
                           "block_sizes": sorted(set(args.block_sizes or [])), "seed": args.seed,
                           "replicates": args.replicates if args.block_sizes else 0,
                           "chromosome": args.chromosome, "interval": args.interval},
            "labels": labels, "input_metadata": {k: v for k, v in source.items()
                                                  if k not in ("genes", "classifications", "labels", "command")},
            "genes": rows, **result, "uncertainty": intervals, "block_diagnostics": blocks,
            "clustering_range": {"estimate": None, "status": "not_estimable",
                                 "reason": "No calibrated range estimator is implemented."},
            "warnings": warnings,
        }
        prefix = str(args.output_prefix)
        outputs = {key: prefix + suffix for key, suffix in (
            ("reference_estimates", ".autocorrelation.tsv"), ("chromosome_estimates", ".chromosomes.tsv"),
            ("genes", ".genes.tsv"), ("uncertainty", ".uncertainty.tsv"),
            ("block_diagnostics", ".blocks.tsv"), ("json", ".json"),
        )}
        plots = self.plot_paths(rows) if args.plot else []
        paths = [Path(p).resolve() for p in [*outputs.values(), *(p for _, p in plots)]]
        if len(set(paths)) != len(paths):
            fail("Output paths must be distinct.")
        if inputs.intersection(paths):
            fail("Output would overwrite an input file.")
        if any(not path.parent.is_dir() for path in paths):
            fail("Output parent directories must exist.")
        payload["output_files"] = outputs
        try:
            if plots:
                payload["plots"] = self.plot(payload, plots)
            for key, path in outputs.items():
                if key != "json":
                    self.write_tsv(path, payload[key], key)
            with open(outputs["json"], "w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2, allow_nan=False)
                handle.write("\n")
        except (OSError, ValueError) as exc:
            fail(f"Cannot write topology autocorrelation outputs: {exc}")
        for warning in warnings:
            print(f"Warning: {warning}", file=sys.stderr)
        if args.json:
            print_json(payload)
        else:
            print("reference\tdistance_start\tdistance_end\tpair_count\texcess_agreement")
            for row in result["reference_estimates"]:
                if row["metric"] == "agreement":
                    print("\t".join("NA" if row[k] is None else str(row[k]) for k in
                                    ("reference", "distance_start", "distance_end", "pair_count", "excess")))

    @staticmethod
    def write_tsv(path, rows, kind):
        if rows:
            TopologyLandscape.write_tsv(path, rows)
            return
        columns = {
            "uncertainty": ["reference", "distance_start", "distance_end", "metric", "requested_block_bp",
                            "replicates", "valid_replicates", "seed", "confidence_level", "interval_type",
                            "lower", "upper", "status", "reasons", "block_size_sensitive"],
            "block_diagnostics": ["reference", "chromosome", "requested_block_bp", "block_index", "start",
                                  "end", "total_genes", "resolved_genes", *CLASSES],
        }
        with open(path, "w", newline="", encoding="utf-8") as handle:
            csv.writer(handle, delimiter="\t").writerow(columns[kind])

    def plot_paths(self, rows):
        path = Path(self.args.plot_output or (str(self.args.output_prefix) + ".png"))
        if path.suffix.lower() not in (".png", ".pdf", ".svg", ".jpg", ".jpeg", ".tif", ".tiff"):
            fail("Plot output must use PNG, PDF, SVG, JPEG, or TIFF format.")
        references = sorted({r["reference"] for r in rows})
        return [(ref, str(path if len(references) == 1 else
                          path.with_name(f"{path.stem}.{i}{path.suffix}")))
                for i, ref in enumerate(references, 1)]

    def plot(self, payload, paths):
        import matplotlib.pyplot as plt
        from matplotlib.colors import is_color_like

        config = PlotConfig.from_args(self.args)
        config.fig_width = config.fig_width if config.fig_width is not None else 11
        config.fig_height = config.fig_height if config.fig_height is not None else 9
        config.resolve()
        if any(not math.isfinite(v) or v <= 0 for v in (config.fig_width, config.fig_height, config.dpi)):
            fail("Plot dimensions and DPI must be finite and positive.")
        colors = config.merge_colors(COLORS[:3])
        if any(not is_color_like(color) for color in colors):
            fail("Invalid plot color.")
        for reference, path in paths:
            estimates = [r for r in payload["reference_estimates"] if r["reference"] == reference]
            intervals = [r for r in payload["uncertainty"] if r["reference"] == reference]
            with plt.rc_context({"font.size": 9}):
                fig, axes = plt.subplots(3, 1, figsize=(config.fig_width, config.fig_height),
                                         sharex=True, layout="constrained")
                try:
                    for k, metric in enumerate(METRICS):
                        ax = axes[0] if k == 0 else axes[1]
                        data = [r for r in estimates if r["metric"] == metric]
                        x = [(r["distance_start"] + r["distance_end"]) / 2 for r in data]
                        y = [r["excess"] if r["excess"] is not None else np.nan for r in data]
                        color = "#333333" if k == 0 else colors[k - 1]
                        label = "All topologies" if k == 0 else payload["labels"][metric]
                        ax.plot(x, y, marker="o", markersize=3, color=color,
                                label=textwrap.fill(label, 30))
                        for size in payload["parameters"]["block_sizes"]:
                            local = [r for r in intervals if r["metric"] == metric and
                                     r["requested_block_bp"] == size]
                            low = [r["lower"] if r["lower"] is not None else np.nan for r in local]
                            high = [r["upper"] if r["upper"] is not None else np.nan for r in local]
                            ax.fill_between(x, low, high, color=color, alpha=0.12)
                    for ax in axes[:2]:
                        ax.axhline(0, color="#777777", linestyle="--", linewidth=0.8)
                    axes[0].set_ylabel("Excess same-topology agreement")
                    axes[1].set_ylabel("Excess topology-specific joint support")
                    aggregate = [r for r in estimates if r["metric"] == "agreement"]
                    axes[2].bar([r["distance_start"] for r in aggregate], [r["pair_count"] for r in aggregate],
                                width=[r["distance_end"] - r["distance_start"] for r in aggregate],
                                align="edge", color="#777777")
                    axes[2].set_ylabel("Resolved gene pairs")
                    axes[2].set_xlabel("Genomic distance between gene anchors (bp)")
                    axes[2].set_xlim(0, payload["parameters"]["maximum_distance_exclusive"])
                    if config.show_title:
                        axes[0].set_title(textwrap.fill(config.title or reference, 75))
                    if config.legend_position != "none":
                        placement = {"loc": config.legend_position} if config.legend_position else {
                            "loc": "upper left", "bbox_to_anchor": (1.01, 1)}
                        axes[1].legend(**placement, fontsize=8)
                    for ax in axes:
                        ax.spines[["top", "right"]].set_visible(False)
                        ax.tick_params(axis="x", labelsize=config.xlabel_fontsize)
                        ax.tick_params(axis="y", labelsize=config.ylabel_fontsize)
                        ax.xaxis.label.set_size(config.axis_fontsize)
                        ax.yaxis.label.set_size(config.axis_fontsize)
                        ax.title.set_size(config.title_fontsize)
                    if intervals:
                        fig.supxlabel("Shading: 95% pointwise intervals across requested block sizes", fontsize=8)
                    fig.savefig(path, dpi=config.dpi)
                finally:
                    plt.close(fig)
        return [path for _, path in paths]
