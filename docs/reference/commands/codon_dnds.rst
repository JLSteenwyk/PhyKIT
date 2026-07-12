.. _cmd-codon_dnds:
.. _command-codon_dnds:

Codon dN/dS
===========

Pairwise synonymous and nonsynonymous substitution-rate estimates

Command identity
----------------

:Canonical command: ``codon_dnds``
:Handler: ``codon_dnds``
:Aliases: dn_ds, dnds, kaks
:Standalone executables: pk_codon_dnds, pk_dn_ds, pk_dnds, pk_kaks
:Categories: Alignment quality & statistics, Evolutionary rate analysis

Runtime interface
-----------------

.. include:: /_generated/commands/codon_dnds.inc

Guidance, interpretation, and examples
--------------------------------------

Estimate pairwise nonsynonymous substitutions per nonsynonymous site (dN),
synonymous substitutions per synonymous site (dS), and their ratio
``omega = dN/dS`` from an in-frame nucleotide alignment.

Four estimators are available:

* ``NG86``: Nei and Gojobori (1986), the default counting method.
* ``LWL85``: Li, Wu, and Luo (1985).
* ``YN00``: Yang and Nielsen (2000), accounting for transition and codon-frequency biases.
* ``ML``: Goldman and Yang (1994) codon-model maximum likelihood. This can be substantially slower than the counting methods.

By default, PhyKIT analyzes every unique sequence pair and reports aggregate
statistics. ``--verbose`` prints a TSV row for each pair. ``--reference``
restricts analysis to pairs containing one exact taxon ID, which is useful
for comparisons against a designated reference and for limiting expensive
ML calculations.

**Input requirements and exclusions**

The input must be an aligned nucleotide coding sequence that starts at codon
position one. Alignment length must be divisible by three, all sequences must
have equal length, and gaps must occupy complete codons (``---``). RNA ``U``
is treated as ``T`` and ``.`` is treated as a gap.

Codons containing a gap or an ambiguous nucleotide are excluded independently
for each sequence pair. Terminal stop codons are also excluded. Internal stop
codons fail by default; ``--stop-policy skip`` explicitly excludes them
pairwise. ``codons_used`` and ``codons_skipped`` in verbose output disclose
the effective data for each estimate.

Estimates that are undefined, saturated, or have ``dS = 0`` are emitted as
``NA`` with an explanatory status instead of producing an infinite or invalid
omega value.

.. important::

   Pairwise omega is a descriptive average across the analyzed sequence and
   pair. An omega greater than one is not, by itself, a statistical test for
   positive selection. This command does not fit branch-specific or
   site-specific selection models.

.. code-block:: shell

   # NG86 summary across all sequence pairs
   phykit codon_dnds codon_alignment.fa

   # Pairwise TSV using YN00
   phykit dnds codon_alignment.fa --method YN00 --verbose

   # Compare every sequence with one reference using genetic code 4
   phykit dnds codon_alignment.fa --reference taxon_a --genetic-code 4 --json

   # Goldman-Yang ML with F61 codon frequencies
   phykit dnds codon_alignment.fa --method ML --codon-frequency F61 --verbose

Verbose columns: |br|
``taxon_a``; ``taxon_b``; ``dN``; ``dS``; ``omega``; ``codons_used``;
``codons_skipped``; ``status``

Options: |br|
*<alignment>*: in-frame codon alignment (required) |br|
*--method*: ``NG86``, ``LWL85``, ``YN00``, or ``ML`` (default: ``NG86``) |br|
*--genetic-code*: NCBI genetic-code ID (default: ``1``) |br|
*--kappa*: transition/transversion ratio used by NG86 (default: ``1.0``) |br|
*--codon-frequency*: ``F1x4``, ``F3x4``, or ``F61`` for ML (default: ``F3x4``) |br|
*--stop-policy*: ``error`` or ``skip`` for internal stops (default: ``error``) |br|
*--reference*: analyze only pairs containing this exact taxon ID |br|
*-v/--verbose*: print one TSV row per sequence pair |br|
*--json*: output structured JSON

References: |br|
Nei M, Gojobori T. 1986. Molecular Biology and Evolution 3:418-426.
doi: `10.1093/oxfordjournals.molbev.a040410 <https://doi.org/10.1093/oxfordjournals.molbev.a040410>`_. |br|
Li WH, Wu CI, Luo CC. 1985. Molecular Biology and Evolution 2:150-174.
doi: `10.1093/oxfordjournals.molbev.a040343 <https://doi.org/10.1093/oxfordjournals.molbev.a040343>`_. |br|
Yang Z, Nielsen R. 2000. Molecular Biology and Evolution 17:32-43.
doi: `10.1093/oxfordjournals.molbev.a026236 <https://doi.org/10.1093/oxfordjournals.molbev.a026236>`_. |br|
Goldman N, Yang Z. 1994. Molecular Biology and Evolution 11:725-736.
doi: `10.1093/oxfordjournals.molbev.a040153 <https://doi.org/10.1093/oxfordjournals.molbev.a040153>`_.
