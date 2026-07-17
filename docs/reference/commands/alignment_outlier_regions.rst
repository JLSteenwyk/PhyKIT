.. _cmd-alignment_outlier_regions:
.. _command-alignment_outlier_regions:

Alignment outlier regions
=========================

Detect localized outlier regions in multiple sequence alignments

Command identity
----------------

:Canonical command: ``alignment_outlier_regions``
:Handler: ``alignment_outlier_regions``
:Aliases: aor, outlier_regions
:Standalone executables: pk_alignment_outlier_regions, pk_aor, pk_outlier_regions
:Categories: Alignment quality & statistics, Homology assessment

Runtime interface
-----------------

.. include:: /_generated/commands/alignment_outlier_regions.inc

Guidance, interpretation, and examples
--------------------------------------

``alignment_outlier_regions`` identifies localized stretches that are unusually
divergent in both dimensions of a multiple sequence alignment: relative to the
other characters in the same alignment columns and relative to other regions in
the same sequence. These stretches can result from sequencing, assembly, gene
prediction, contamination, or alignment errors.

The command reports individual taxon-region combinations. It does not remove an
entire sequence or alignment column. With ``--mask-output``, only flagged
characters are replaced, while gaps, existing ambiguous characters, unaffected
taxa, and unaffected columns are retained.

Method
######

PhyKIT uses a two-dimensional procedure inspired by TAPER:

1. Each valid character ``x`` in column ``i`` receives the divergence score
   ``1 / (u_i * p_i,x)``, where ``u_i`` is the number of unique valid character
   states and ``p_i,x`` is the frequency of ``x`` among valid characters in the
   column. Rare characters in otherwise conserved columns therefore receive
   high scores.
2. Scores are summarized within overlapping windows along each ungapped
   sequence.
3. A two-class Jenks natural break separates ordinary and high-scoring windows
   within each sequence.
4. Across-taxon, within-sequence, and absolute score thresholds prevent every
   sequence from being forced to contain an outlier.
5. Dynamic programming smooths window classifications into contiguous regions.

The analysis combines the published TAPER scale settings: window sizes 5, 9,
and 17; taxon-tail proportions 0.25, 0.25, and 0.10; and sequence-tail
proportions 0.10, 0.25, and 0.50. Calls from the three scales are combined.
Regions longer than 30 or 54 valid characters are ignored at window sizes 5 or
9, respectively; the window-size-17 analysis has no upper limit.

This is an independent implementation of the method described in the TAPER
paper, not a wrapper around the TAPER software. Exact calls can differ between
the programs. Cite TAPER when using this command:

   Zhang C, Zhao Y, Braun EL, and Mirarab S. 2021. TAPER: Pinpointing errors in
   multiple sequence alignments despite varying rates of evolution. *Methods in
   Ecology and Evolution* 12:2145-2158.
   https://doi.org/10.1111/2041-210X.13696

Input and missing data
######################

The input must be a multiple sequence alignment. PhyKIT accepts the alignment
formats supported elsewhere in the toolkit, including FASTA, PHYLIP, Clustal,
MAF, and Stockholm.

For nucleotide alignments, ``-``, ``?``, ``*``, ``X``, and ``N`` are treated as
missing. For protein alignments, ``-``, ``?``, ``*``, and ``X`` are treated as
missing. Missing characters do not contribute to column frequencies, window
scores, or sequence coordinates.

Output
######

The default tab-separated report contains one row per detected region:

``taxon``
   Sequence identifier.

``alignment_start`` and ``alignment_end``
   One-based, closed coordinates in the original alignment.

``sequence_start`` and ``sequence_end``
   One-based, closed coordinates among valid characters after gaps and ambiguous
   characters are excluded.

``length``
   Number of valid characters in the region.

``mean_divergence_score`` and ``max_divergence_score``
   Character-level divergence summaries for the region.

``window_sizes``
   TAPER-inspired scales that supported at least part of the reported region.

If no region is detected, the tabular output contains only its header. JSON
output additionally records the citation, parameters, alignment dimensions,
number of affected taxa, and number of masked characters.

Examples
########

Print a tab-separated region report:

.. code-block:: shell

   phykit alignment_outlier_regions alignment.fa

Write the report and a masked FASTA alignment:

.. code-block:: shell

   phykit alignment_outlier_regions alignment.fa \
     --report suspicious_regions.tsv \
     --mask-output alignment.masked.fa

Use ``N`` for masked nucleotide characters and write JSON metadata:

.. code-block:: shell

   phykit alignment_outlier_regions alignment.fa \
     --mask-output alignment.masked.fa \
     --mask-character N \
     --json

Increase sensitivity by lowering the absolute cutoff while keeping it greater
than 1:

.. code-block:: shell

   phykit alignment_outlier_regions alignment.fa --cutoff 2.5

Interpretation and limitations
##############################

A reported region is a statistical outlier, not proof that the underlying
characters are nonhomologous. Inspect important calls and compare downstream
results with and without masking. Highly divergent but valid biological regions
can be flagged, especially when taxon sampling is sparse.

As discussed for TAPER, two-dimensional detection has limited power for very
short errors. Very long errors, or the same error occurring in many taxa, can
look like normal biological variation rather than an outlier. The method also
assumes that the alignment mostly represents the same homologous locus; resolve
clear paralogy before interpreting localized calls.

Lower ``--cutoff`` values are more aggressive and can increase false positives.
The default value of 3.0 follows the TAPER recommendation. Masking should be
treated as a sensitivity-analysis step rather than an automatic guarantee of a
better phylogenetic estimate.

Related commands
################

- :doc:`alignment_outlier_taxa` evaluates whole-sequence quality features.
- :doc:`spurious_sequence` detects tree-based long-branch and diameter-impact
  outliers.
- :doc:`hidden_paralogy_check` tests expected monophyletic groups for evidence
  of hidden paralogy.
- :doc:`mask_alignment` removes complete alignment columns using occupancy,
  gap, or entropy thresholds.

