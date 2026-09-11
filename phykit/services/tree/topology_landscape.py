"""Map strict four-group gene-tree concordance to reference coordinates."""

import csv
import json
import math
import sys
import textwrap
from pathlib import Path

from ...helpers.genomic_coordinates import (
    neighborhoods,
    read_coordinates,
    read_manifest,
    select_rows,
    summarize,
)
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig
from ...helpers.topology_landscape import (
    CLASSES,
    classify_tree,
    fail,
    groups_from_branch,
    read_tree,
    validate_groups,
)

COLORS = ["#0072B2", "#D55E00", "#009E73", "#999999", "#E6C84F", "#CC79A7"]


class TopologyLandscape:
    def __init__(self, args):
        self.args = args

    def prepare(self):
        """Classify and select mapped genes without writing files or neighborhoods."""
        args = self.args
        if bool(args.groups) == bool(args.reference_tree):
            fail("Provide either --groups or --reference-tree with --branch-taxa.")
        if bool(args.reference_tree) != bool(args.branch_taxa):
            fail("--reference-tree and --branch-taxa must be provided together.")
        if args.groups:
            try:
                with open(args.groups, encoding="utf-8") as handle:
                    def unique_object(pairs):
                        obj = {}
                        for key, value in pairs:
                            if key in obj:
                                fail(f"Duplicate group name: {key}.")
                            obj[key] = value
                        return obj
                    groups = json.load(handle, object_pairs_hook=unique_object)
            except (OSError, ValueError, UnicodeError) as exc:
                fail(f"Cannot read group JSON: {exc}")
            if not isinstance(groups, dict):
                fail("Group JSON must be an object mapping four names to taxon lists.")
            groups = validate_groups(groups)
        else:
            groups = groups_from_branch(read_tree(args.reference_tree), args.branch_taxa)
        if args.outgroup:
            if args.outgroup not in groups:
                fail("--outgroup must name one of the four groups.")
            groups = {args.outgroup: groups[args.outgroup], **{
                name: taxa for name, taxa in groups.items() if name != args.outgroup
            }}
        names = list(groups)
        labels = {}
        for index, partner in enumerate((1, 2, 3), 1):
            other = [names[j] for j in range(1, 4) if j != partner]
            labels[f"topology_{index}"] = (
                f"{names[partner]}-sister" if args.outgroup
                else f"{names[0]} + {names[partner]} | {' + '.join(other)}"
            )
        if args.labels:
            labels.update(zip(CLASSES[:3], args.labels))
        labels.update(unresolved="Unresolved", insufficient_sampling="Insufficient sampling",
                      incompatible_groups="Incompatible groups")
        manifest = read_manifest(args.manifest)
        coordinates, diagnostics = read_coordinates(args.coordinates)
        located = {r["gene_id"] for r in coordinates}
        diagnostics["unmapped_genes"] = sorted(set(manifest) - located)
        diagnostics["coordinate_genes_without_trees"] = sorted(located - set(manifest))
        diagnostics["unmapped_by_reference"] = {
            ref: sorted(set(manifest) - {r["gene_id"] for r in coordinates if r["reference"] == ref})
            for ref in diagnostics["references"]
        }
        if diagnostics["unmapped_genes"] and args.unmapped == "error":
            fail("Genes lack coordinates: " + ", ".join(diagnostics["unmapped_genes"]) +
                 ". Use --unmapped skip to report and omit them from maps.")
        for field in ("unmapped_genes", "coordinate_genes_without_trees"):
            if diagnostics[field]:
                print(f"Warning: {field}: {len(diagnostics[field])}; see diagnostics JSON.", file=sys.stderr)
        if diagnostics["merged_records"]:
            print(f"Warning: merged repeated gene records: {diagnostics['merged_records']}", file=sys.stderr)
        classifications = {}
        for gene, path in manifest.items():
            classifications[gene] = classify_tree(
                read_tree(path), groups, args.min_support, args.support_scale, args.missing_support,
            )
        rows = [dict(r, **classifications[r["gene_id"]]) for r in coordinates if r["gene_id"] in manifest]
        rows = select_rows(rows, args.chromosome, args.interval)
        if not rows:
            fail("No mapped genes remain in the selected chromosomes/interval.")
        inputs = {Path(p).resolve() for p in [args.manifest, *manifest.values(),
                  *(s[2] for s in args.coordinates), args.groups or args.reference_tree]}
        return {
            "command": "topology_landscape", "mode": "strict", "coordinate_system": "0-based-half-open",
            "groups": {name: sorted(taxa) for name, taxa in groups.items()}, "labels": labels,
            "outgroup": args.outgroup, "min_support": args.min_support, "support_scale": args.support_scale,
            "missing_support": args.missing_support, "genes": rows, "classifications": classifications,
            "diagnostics": diagnostics,
        }, inputs

    def run(self):
        args = self.args
        if args.window_bp is not None and args.window_bp <= 0:
            fail("--window-bp must be positive.")
        if args.window_genes is not None and args.window_genes <= 0:
            fail("--window-genes must be positive.")
        payload, inputs = self.prepare()
        rows, labels, diagnostics = payload["genes"], payload["labels"], payload["diagnostics"]
        windows = neighborhoods(
            rows, window_bp=args.window_bp if args.window_genes is None else None,
            window_genes=args.window_genes, interval=args.interval,
        )
        overall = [dict(reference=ref, **summarize([r for r in rows if r["reference"] == ref]))
                   for ref in sorted({r["reference"] for r in rows})]
        payload.update(neighborhoods=windows, overall=overall)
        prefix = Path(args.output_prefix)
        outputs = {
            key: str(prefix) + suffix for key, suffix in (
                ("genes", ".genes.tsv"), ("neighborhoods", ".neighborhoods.tsv"),
                ("overall", ".overall.tsv"), ("diagnostics", ".diagnostics.json"),
            )
        }
        # Do not let a mistyped prefix or plot filename overwrite any source.
        plot_paths = self.plot_paths(rows) if args.plot else []
        for path in [*outputs.values(), *(p[2] for p in plot_paths)]:
            if Path(path).resolve() in inputs:
                fail(f"Output would overwrite an input file: {path}")
        try:
            if args.plot:
                payload["plots"] = self.plot(rows, windows, labels, plot_paths)
            self.write_tsv(outputs["genes"], rows)
            self.write_tsv(outputs["neighborhoods"], windows)
            self.write_tsv(outputs["overall"], overall)
            with open(outputs["diagnostics"], "w", encoding="utf-8") as handle:
                json.dump(diagnostics, handle, indent=2)
                handle.write("\n")
        except OSError as exc:
            fail(f"Cannot write topology landscape outputs: {exc}")
        payload["output_files"] = outputs
        if args.json:
            print_json(payload)
        else:
            print("reference\ttotal_genes\tinformative_genes\ttopology_1\ttopology_2\ttopology_3")
            for record in overall:
                print("\t".join(str(record[k]) for k in (
                    "reference", "total_genes", "informative_genes", *CLASSES[:3],
                )))

    @staticmethod
    def write_tsv(path, rows):
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow({key: json.dumps(value, sort_keys=True) if isinstance(value, (dict, list))
                                 else value for key, value in row.items()})

    def plot_paths(self, rows):
        path = Path(self.args.plot_output or (self.args.output_prefix + ".png"))
        if path.suffix.lower() not in (".png", ".pdf", ".svg", ".jpg", ".jpeg", ".tif", ".tiff"):
            fail("Plot output must use PNG, PDF, SVG, JPEG, or TIFF format.")
        panels = sorted({(r["reference"], r["chromosome"]) for r in rows})
        return [(ref, chrom, str(path if len(panels) == 1 else
                path.with_name(f"{path.stem}.{index}{path.suffix}")))
                for index, (ref, chrom) in enumerate(panels, 1)]

    def plot(self, rows, windows, labels, paths):
        import matplotlib.pyplot as plt
        from matplotlib.colors import is_color_like
        from matplotlib.lines import Line2D

        config = PlotConfig.from_args(self.args)
        config.fig_width = config.fig_width or 13
        config.fig_height = config.fig_height or 8
        config.resolve()
        if any(not math.isfinite(value) or value <= 0 for value in
               (config.fig_width, config.fig_height, config.dpi)):
            fail("Plot dimensions and DPI must be finite and positive.")
        colors = config.merge_colors(COLORS)
        if any(not is_color_like(c) for c in colors):
            fail("Invalid plot color.")
        outputs = []
        for ref, chromosome, path in paths:
            genes = [r for r in rows if r["reference"] == ref and r["chromosome"] == chromosome]
            local = [w for w in windows if w["reference"] == ref and w["chromosome"] == chromosome]
            with plt.rc_context({"font.size": 9}):
                fig, axes = plt.subplots(3, 1, figsize=(config.fig_width, config.fig_height),
                                         gridspec_kw={"height_ratios": [1.5, 1, 1]},
                                         layout="constrained")
                try:
                    track, proportions, counts = axes
                    for index, category in enumerate(CLASSES):
                        members = [r for r in genes if r["classification"] == category]
                        track.scatter([r["anchor"] for r in members], [index] * len(members),
                                      marker="|", s=100, color=colors[index])
                    track.set_yticks(range(6), [textwrap.fill(labels[c], 28) for c in CLASSES])
                    track.invert_yaxis()
                    track.set_xlabel("Reference position (bp); gene-span midpoint")
                    if config.show_title:
                        track.set_title(textwrap.fill(config.title or f"{ref}: {chromosome}", 85),
                                        fontsize=config.title_fontsize)
                    gene_mode = local[0]["window_mode"] == "genes"
                    x = [w["window_index"] if gene_mode else w["start"] for w in local]
                    width = [0.85 if gene_mode else w["end"] - w["start"] for w in local]
                    bottom = [0.] * len(local)
                    align = "center" if gene_mode else "edge"
                    for index, category in enumerate(CLASSES[:3]):
                        heights = [w[f"{category}_proportion"] or 0 for w in local]
                        proportions.bar(x, heights, width=width, bottom=bottom, align=align,
                                        color=colors[index], linewidth=0)
                        bottom = [b + h for b, h in zip(bottom, heights)]
                    proportions.set_ylim(0, 1)
                    proportions.set_ylabel("Fraction of resolved genes")
                    counts.bar(x, [w["total_genes"] for w in local], width=width,
                               align=align, color="#BBBBBB", label="All mapped genes")
                    counts.bar(x, [w["informative_genes"] for w in local], width=width,
                               align=align, color="#333333", label="Resolved genes")
                    counts.set_ylabel("Gene count")
                    counts.set_xlabel("Gene-count window index" if gene_mode else "Reference position (bp)")
                    if config.legend_position != "none":
                        placement = ({"loc": config.legend_position} if config.legend_position
                                     else {"loc": "upper left", "bbox_to_anchor": (1.01, 1)})
                        counts.legend(**placement, fontsize=8)
                        proportions.legend(handles=[Line2D([], [], color=colors[i], lw=3,
                                           label=textwrap.fill(labels[c], 28)) for i, c in enumerate(CLASSES[:3])],
                                           **placement, fontsize=8)
                    for ax in axes:
                        ax.spines[["top", "right"]].set_visible(False)
                        ax.tick_params(axis="x", labelsize=config.xlabel_fontsize)
                        ax.tick_params(axis="y", labelsize=config.ylabel_fontsize)
                        ax.xaxis.label.set_size(config.axis_fontsize)
                        ax.yaxis.label.set_size(config.axis_fontsize)
                    bounds = self.args.interval or (0, max(r["end"] for r in genes))
                    track.set_xlim(*bounds)
                    if not gene_mode:
                        proportions.set_xlim(*bounds)
                        counts.set_xlim(*bounds)
                    fig.savefig(path, dpi=config.dpi)
                finally:
                    plt.close(fig)
            outputs.append({"reference": ref, "chromosome": chromosome, "path": path})
        return outputs
