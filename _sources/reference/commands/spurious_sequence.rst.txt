.. _cmd-spurious_sequence:
.. _command-spurious_sequence:

Spurious homolog identification
===============================

Identify potentially spurious sequences from phylogenetic trees

Command identity
----------------

:Canonical command: ``spurious_sequence``
:Handler: ``spurious_sequence``
:Aliases: spurious_seq, ss
:Standalone executables: pk_spurious_sequence, pk_spurious_seq, pk_ss
:Categories: Homology assessment

Runtime interface
-----------------

.. include:: /_generated/commands/spurious_sequence.inc

Guidance, interpretation, and examples
--------------------------------------

``spurious_sequence`` provides two methods for identifying tips whose branch
lengths are inconsistent with the rest of a gene tree. Both methods report
candidates for review; neither edits an alignment or tree.

Median-factor method
^^^^^^^^^^^^^^^^^^^^

``--method median-factor`` is the default and preserves the original command
behavior. It reports terminal branches that are at least ``--factor`` times
the median terminal-branch length. The default factor is 20.

Text output has four tab-delimited columns:

1. taxon name
2. terminal branch length
3. long-branch threshold
4. median terminal-branch length

If there are no putatively spurious sequences, "None" is reported.

Using this method to identify potentially spurious sequences
was, to my knowledge, first introduced by Shen et al., (2018)
Cell doi: 10.1016/j.cell.2018.10.023. 

.. code-block:: shell

   phykit spurious_sequence gene.tre
   phykit spurious_sequence gene.tre --factor 10 --json

Diameter-impact method
^^^^^^^^^^^^^^^^^^^^^^

``--method diameter-impact`` detects tips that disproportionately increase a
tree's weighted diameter, the greatest branch-length distance between any two
active tips. At each iteration, PhyKIT evaluates both endpoints of the current
diameter and removes the endpoint whose removal produces the smaller next
diameter. The tip signature is
``log(diameter_before / diameter_after)``. This repeats for
``--max-remove`` tips; the default is
``min(floor(n / 4), floor(5 * sqrt(n)))``, capped so that two tips remain.

For one tree, PhyKIT fits a one-sided normal model to all log-diameter
signatures, including zero-impact tips. Equivalently, this models the
multiplicative diameter ratios with a log-normal distribution. A tip is
reported when its upper-tail p-value is less than ``--alpha`` (default 0.05).

This is an independently implemented, deterministic TreeShrink-inspired
endpoint-removal heuristic. It uses TreeShrink's diameter-reduction concept
and statistical calibration, but it does not reproduce TreeShrink's full
optimal removal-set search. See Mai and Mirarab (2018), BMC Genomics 19:272,
doi: 10.1186/s12864-018-4620-2.

The input must contain at least four uniquely named tips. Every non-root
branch must have a finite, nonnegative length. Rooting does not affect the
weighted tip-to-tip distances.

Single-tree text output has seven tab-delimited columns:

1. taxon name
2. log-diameter reduction signature
3. upper-tail p-value
4. removal rank
5. terminal branch length
6. diameter before removal
7. diameter after removal

.. code-block:: shell

   phykit spurious_sequence gene.tre \
       --method diameter-impact --alpha 0.05

Tree collections
^^^^^^^^^^^^^^^^

Use ``--tree-list`` when the positional input is a file containing one tree
path per line. Blank lines and lines beginning with ``#`` are ignored. Relative
tree paths are resolved from the list file's directory.

By default, signatures from all listed trees are pooled. PhyKIT estimates
upper-tail probabilities with a Gaussian kernel density estimate using
Silverman's bandwidth and leave-one-out evaluation. ``--per-species`` instead
builds a separate distribution for each taxon across trees in which it occurs.
Groups with fewer than four observations are not called as outliers.

Collection text output prepends the tree path to the seven columns described
above.

.. code-block:: shell

   phykit spurious_sequence trees.txt \
       --method diameter-impact --tree-list --per-species --json

JSON output
^^^^^^^^^^^

The median-factor JSON contract is unchanged: ``rows`` and
``spurious_sequences`` contain identical candidate objects with ``taxon``,
``branch_length``, ``threshold``, and ``median`` fields.

Diameter-impact JSON adds ``method``, ``scope``, ``alpha``, and ``tree_count``
metadata. ``rows`` and ``spurious_sequences`` contain identical objects with
the tree path, taxon, signature, p-value, removal rank, terminal branch length,
and before/after diameters.

.. code-block:: json

   {
     "alpha": 0.05,
     "method": "diameter-impact",
     "rows": [
       {
         "diameter_after": 2.0,
         "diameter_before": 101.0,
         "p_value": 0.000007,
         "removal_rank": 1,
         "signature": 3.921973,
         "taxon": "outlier",
         "terminal_branch_length": 100.0,
         "tree": "gene.tre"
       }
     ],
     "scope": "per-gene",
     "spurious_sequences": [
       {
         "diameter_after": 2.0,
         "diameter_before": 101.0,
         "p_value": 0.000007,
         "removal_rank": 1,
         "signature": 3.921973,
         "taxon": "outlier",
         "terminal_branch_length": 100.0,
         "tree": "gene.tre"
       }
     ],
     "tree_count": 1
   }

.. code-block:: shell

   phykit spurious_seq gene.tre --factor 20 --json

Options: |br|
*<file>*: first argument after function name should be a tree file |br|
*-f/--factor*: factor to multiply median branch length by to calculate
the threshold of long branches. (Default: 20) |br|
*--method*: ``median-factor`` or ``diameter-impact`` (default:
``median-factor``) |br|
*--alpha*: upper-tail significance level for ``diameter-impact`` (default:
``0.05``) |br|
*--max-remove*: maximum diameter endpoints evaluated; by default follows the
TreeShrink-inspired rule described above |br|
*--tree-list*: interpret the positional input as a tree-path list and pool
signatures across trees |br|
*--per-species*: with ``--tree-list``, calibrate signatures independently for
each taxon |br|
*--json*: optional argument to print results as JSON
