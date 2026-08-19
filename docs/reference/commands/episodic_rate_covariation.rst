.. _cmd-episodic_rate_covariation:
.. _command-episodic_rate_covariation:

Episodic rate covariation
==========================

Localize coordinated gene-rate shifts to clades on a reference tree

Command identity
----------------

:Canonical command: ``episodic_rate_covariation``
:Handler: ``episodic_rate_covariation``
:Aliases: erc_scan
:Standalone executables: pk_episodic_rate_covariation, pk_erc_scan
:Categories: Evolutionary rate analysis
:Status: Experimental

Runtime interface
-----------------

.. include:: /_generated/commands/episodic_rate_covariation.inc

Guidance, interpretation, and examples
--------------------------------------

Use ``episodic_rate_covariation`` to ask *where* on a reference tree two genes
show coordinated evolutionary-rate changes. The command complements a single
whole-tree correlation by scanning nonterminal rooted clades for localized
concordant or antagonistic rate shifts.

This is a novel experimental PhyKIT method. It is not the published CovER
branch-ratio statistic, a branch-site selection test, or evidence by itself of
a direct molecular interaction. Treat discoveries as hypotheses for
independent biological and statistical validation.

Inputs and preprocessing
^^^^^^^^^^^^^^^^^^^^^^^^

Provide two gene trees and one rooted reference tree in Newick format. All
trees must have named tips and finite, nonnegative branch lengths. The
reference is typically a species tree whose root defines the clades being
tested.

PhyKIT takes the intersection of the taxa in all three trees and prunes each
tree to that shared set. It then applies the same reference-edge projection,
rate filtering, and standardization used by
:doc:`projected_covarying_rates`. Gene trees may therefore disagree in
topology and may contain taxa absent from the other inputs. Inspect both
projection NRMSE values before interpreting any local scan.

Clade statistic
^^^^^^^^^^^^^^^

Let :math:`x_i` and :math:`y_i` be the standardized projected rates of the two
genes on retained reference edge :math:`i`. Its local correlation
contribution is

.. math::

   c_i = x_i y_i.

For a candidate clade :math:`C` containing :math:`k` of the :math:`n`
localizable edges, PhyKIT calculates

.. math::

   T_C = \frac{\sum_{i \in C}(c_i-\bar{c})}
   {\sqrt{k(1-k/n)}}.

This size-normalized statistic contrasts the clade's contributions with the
background tree. The output also reports the mean contribution inside the
clade, the mean outside it, and their difference. ``coupling`` is
``concordant`` when the local mean is nonnegative and ``antagonistic`` when it
is negative. The sign of ``score`` instead indicates whether the local mean is
greater or less than the background mean; these two signs answer different
questions.

A candidate contains its stem edge and all descendant edges. Nonterminal
clades are tested only when at least ``--min-edges`` localizable edges occur
both inside and outside the clade. Terminal tips and the complete tree are not
tested.

Pairwise distances identify only the sum of the two complementary branches
below a degree-two root. That merged unrooted edge contributes to the global
correlation but cannot be assigned to either rooted child clade, so it is
excluded from the local scan. ``ambiguous_root_edge_count`` records this
exclusion.

Permutation test and FWER control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The null test permutes the second gene's standardized rates among edges with
similar topological root depth while holding the first gene fixed. Equal
depths are never split across bins. When there are more unique depths than
``--depth-bins``, quantile boundaries combine adjacent depths into at most the
requested number of strata.

For :math:`B` permutations, the unadjusted two-sided empirical p-value for a
clade is

.. math::

   p_C = \frac{1 + \#\{|T^{(b)}_C| \geq |T_C|\}}{B+1}.

PhyKIT controls the family-wise error rate across all candidate clades with a
maximum absolute statistic:

.. math::

   p_{C,\mathrm{FWER}} =
   \frac{1 + \#\{\max_D |T^{(b)}_D| \geq |T_C|\}}{B+1}.

``significant`` is based on this FWER-adjusted value and ``--alpha``. The
smallest possible p-value is ``1 / (--permutations + 1)``. Use more than the
default 999 permutations when fine tail resolution is important. ``--seed``
makes both distance-pair subsampling and permutation results reproducible.

The null assumes that edges are exchangeable within each depth stratum. This
constraint reduces broad depth effects but cannot account for every source of
branch dependence. Nested clades also share edges. FWER correction protects
the tested family under the implemented permutation scheme; it does not make
the biological hypotheses independent.

Global and edge diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The reported ``correlation`` and ``p_value`` describe the complete retained
projected-rate vectors. The Pearson p-value remains descriptive because
reference edges are not independent. A local clade can be notable even when
the global coefficient is modest, and a strong global coefficient need not be
concentrated in one clade.

Every retained branch receives:

- ``local_contribution``: :math:`x_i y_i`;
- ``leave_one_out_correlation``: the global correlation after removing that
  edge;
- ``correlation_delta_without_edge``: leave-one-out correlation minus the
  full correlation.

``most_influential_edge`` identifies the retained edge with the largest
absolute leave-one-out change. These are sensitivity diagnostics, not
independent branch-level tests.

Output
^^^^^^

Default text output summarizes the global correlation, projection and scan
edge counts, depth bins, permutation count, significant-clade count, and most
influential edge. It then reports FWER-significant clades. If none are
significant, all candidates are shown so the null result remains inspectable;
``--verbose`` always shows every candidate.

``--output`` writes one TSV row per candidate with the score, local and
background means, contribution contrast, coupling direction, unadjusted
p-value, FWER p-value, and significance call. ``--annotated-tree`` writes a
copy of the reference tree whose significant clades carry Newick comments
named ``erc_score``, ``erc_fwer_p``, and ``erc_coupling``.

``--json`` returns complete structured diagnostics. Important fields include:

- ``clade_scans``: every candidate and its local statistics;
- ``candidate_count`` and ``significant_clade_count``;
- ``scan_edge_count``, ``permutable_edge_count``, and
  ``ambiguous_root_edge_count``;
- ``requested_depth_bins`` and ``effective_depth_bins``;
- ``branches`` and ``most_influential_edge``;
- both projection NRMSE records and the exact shared-taxon set.

Examples
^^^^^^^^

Run a reproducible exploratory scan and save the full result:

.. code-block:: shell

   phykit episodic_rate_covariation gene_a.tre gene_b.tre \
       --reference species.tre --permutations 9999 --seed 814 \
       --output gene_a_gene_b_clades.tsv --json > gene_a_gene_b_scan.json

Write significant discoveries onto the reference tree:

.. code-block:: shell

   phykit erc_scan gene_a.tre gene_b.tre -r species.tre \
       --min-edges 4 --alpha 0.05 \
       --annotated-tree gene_a_gene_b_episodic.tre

Use :doc:`projected_covarying_rates` when only a topology-robust whole-tree
coefficient is needed, or :doc:`covarying_evolutionary_rates` when the gene
trees were constrained to the same rooted topology and the published CovER
branch-ratio calculation is required.

Citations and validation
^^^^^^^^^^^^^^^^^^^^^^^^

The evolutionary-rate covariation framework follows `Clark et al. (2012)
<https://doi.org/10.1101/gr.132647.111>`_ and `Steenwyk et al. (2022)
<https://doi.org/10.1126/sciadv.abn0105>`_. Reference-edge projection and the
depth-matched maximum-clade scan are experimental PhyKIT methodology and
should be described as such rather than attributed to those publications.

PhyKIT's validation suite checks exact projection on additive synthetic trees,
recovery of a planted clade-level covariation event after FWER correction, and
the absence of a discovery in a fixed independent-rate null simulation. These
tests verify implementation behavior; they do not establish power or false
positive rates for every tree shape, taxon count, or evolutionary process.
