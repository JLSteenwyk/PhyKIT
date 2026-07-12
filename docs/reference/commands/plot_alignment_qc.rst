.. _cmd-plot_alignment_qc:
.. _command-plot_alignment_qc:

Plot alignment QC
=================

Visual quality control plots for alignments

Command identity
----------------

:Canonical command: ``plot_alignment_qc``
:Handler: ``plot_alignment_qc``
:Aliases: paqc, plot_qc
:Standalone executables: pk_plot_alignment_qc, pk_paqc, pk_plot_qc
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/plot_alignment_qc.inc

Guidance, interpretation, and examples
--------------------------------------

Generate a multi-panel alignment quality-control plot.

The figure includes:
1) occupancy per taxon |br|
2) gap rate per taxon |br|
3) composition distance vs long-branch proxy scatter |br|
4) count of flagged outliers by feature

Outlier evaluation uses the same features as ``alignment_outlier_taxa``:
``gap_rate``, ``occupancy``, ``composition_distance``, ``long_branch_proxy``,
``rcvt``, and ``entropy_burden``.

.. code-block:: shell

	phykit plot_alignment_qc <alignment> [-o/--output <path>] [--width <float>] [--height <float>] [--dpi <int>] [--gap-z <float>] [--composition-z <float>] [--distance-z <float>] [--rcvt-z <float>] [--occupancy-z <float>] [--entropy-z <float>] [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*-o/--output*: output image path (default: ``alignment_qc.png``) |br|
*--width*: figure width in inches (default: ``14.0``) |br|
*--height*: figure height in inches (default: ``10.0``) |br|
*--dpi*: output image DPI (default: ``300``) |br|
*--gap-z*: z-threshold for gap-rate outliers (default: ``3.0``) |br|
*--composition-z*: z-threshold for composition-distance outliers (default: ``3.0``) |br|
*--distance-z*: z-threshold for long-branch-proxy outliers (default: ``3.0``) |br|
*--rcvt-z*: z-threshold for RCVT outliers (default: ``3.0``) |br|
*--occupancy-z*: z-threshold for low-occupancy outliers (default: ``3.0``) |br|
*--entropy-z*: z-threshold for entropy-burden outliers (default: ``3.0``) |br|
*--json*: optional argument to print plot metadata and outlier summary as JSON
