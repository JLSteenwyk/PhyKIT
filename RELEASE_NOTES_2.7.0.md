# PhyKIT 2.7.0

## Genomic Topology Landscapes

New command: `topology_landscape` (alias `topomap`), also available as
`pk_topology_landscape` and `pk_topomap`.

- Classify gene trees into three focal quartet resolutions, unresolved,
  insufficiently sampled, or incompatible group relationships.
- Define four taxon groups explicitly or select an internal reference-tree edge.
  Classification uses every sampled group representative and unrooted splits.
- Map gene identifiers to BED or TSV coordinates across separate reference
  genomes, with duplicate-record and missing-mapping diagnostics.
- Summarize physical-distance or fixed-gene-count neighborhoods with explicit
  total and resolved-gene denominators.
- Export per-gene, neighborhood, and overall TSVs, diagnostics JSON, complete
  JSON results, and chromosome figures in raster or vector formats.
- Configure support thresholds, missing-support policies, genomic selections,
  hypothesis labels, colors, and plot dimensions.

Tutorial 22 includes downloadable synthetic animal-root data demonstrating
ctenophore-sister, sponge-sister, and ctenophore-plus-sponge classifications.

This is descriptive gene-tree concordance mapping, not a topology significance
test, introgression test, or recombination-breakpoint caller. It does not
implement Twisst's fractional subtree weighting. GFF3 and automatic
transcript-to-gene identifier conversion are not included in this release.

Related method: Martin and Van Belleghem (2017), Genetics,
doi:10.1534/genetics.116.194720.

Documentation deployment now follows `main`, the branch containing this release.
