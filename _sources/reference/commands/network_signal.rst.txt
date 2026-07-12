.. _cmd-network_signal:
.. _command-network_signal:

Network signal
==============

Phylogenetic signal on networks

Command identity
----------------

:Canonical command: ``network_signal``
:Handler: ``network_signal``
:Aliases: net_signal, netsig
:Standalone executables: pk_network_signal, pk_net_signal, pk_netsig
:Categories: Phylogenetic signal

Runtime interface
-----------------

.. include:: /_generated/commands/network_signal.inc

Guidance, interpretation, and examples
--------------------------------------

Compute phylogenetic signal (Blomberg's K and/or Pagel's lambda)
on a **phylogenetic network** rather than a tree. This accounts for
hybridization and introgression when estimating how strongly a
continuous trait tracks evolutionary history.

Standard phylogenetic signal methods assume a strictly bifurcating tree.
When the true history includes reticulation, the tree-based
variance-covariance (VCV) matrix is incorrect and signal estimates
may be biased. ``network_signal`` replaces the tree VCV with a
**network VCV** computed using the recursive algorithm of
Bastide et al. (*Systematic Biology*, 2018), which properly
weights shared ancestry through both tree-like and hybrid lineages.

Polytomies (collapsed branches) in the input tree are represented as star
topologies in the network VCV, which correctly models unresolved
relationships as equal covariance among all children.

Two signal metrics are available (same as ``phylogenetic_signal``):

- **Blomberg's K** (Blomberg et al. 2003): K = 1 under Brownian motion;
  K < 1 = less signal than expected; K > 1 = more. P-value via
  permutation test. Computing K on a network is a **novel capability**
  not available in any other tool.
- **Pagel's lambda** (Pagel 1999): lambda = 0 = no signal; lambda = 1 =
  full BM signal. P-value via likelihood ratio test.

**Network specification** — two options:

1. **Explicit hybrid edges** (``--hybrid``): specify one or more
   reticulation events as ``donor:recipient:gamma`` where gamma is
   the inheritance probability from the donor lineage (0 < gamma < 0.5).
2. **From quartet_network JSON** (``--quartet-json``): auto-infer
   hybrid edges from the output of ``phykit quartet_network --json``.
   The command identifies taxon pairs that swap across hybrid quartets
   and estimates gamma from concordance factor ratios.

.. code-block:: shell

   # With explicit hybrid edges
   phykit network_signal -t <tree> -d <trait_data> --hybrid <donor:recipient:gamma> [--method both|blombergs_k|lambda] [--permutations 1000] [--json]

   # With quartet_network JSON output
   phykit network_signal -t <tree> -d <trait_data> --quartet-json <quartets.json> [--method both|blombergs_k|lambda] [--permutations 1000] [--json]

Options: |br|
*-t/--tree*: a rooted species tree in Newick format (with branch lengths) |br|
*-d/--trait-data*: tab-delimited trait file (taxon_name<tab>trait_value) |br|
*--hybrid*: one or more hybrid edge specifications (donor:recipient:gamma); |br|
donor is the source lineage, recipient receives gene flow, gamma is the |br|
inheritance proportion from the donor (e.g., ``B:C:0.3``) |br|
*--quartet-json*: path to JSON output from ``phykit quartet_network --json`` |br|
*--method*: ``both`` (default), ``blombergs_k``, or ``lambda`` |br|
*--permutations*: number of permutations for K p-value (default: 1000) |br|
*-v/--verbose*: print network VCV matrix details |br|
*--json*: optional argument to print results as JSON

Output for default (both) mode: |br|
``Hybrid edge: B -> C (gamma=0.3000)`` |br|
``Network taxa: 5`` |br|
``---`` |br|
``Blomberg's K: 0.8234    p-value: 0.0320`` |br|
``Pagel's lambda: 0.7651    log-likelihood: -12.3456    p-value: 0.0012``

|

**Tutorial: Wing pattern evolution in** *Heliconius* **butterflies**

This example shows a realistic workflow for computing phylogenetic
signal on a network, starting from gene tree discordance analysis
through to signal estimation. The scenario is motivated by the
*Heliconius* butterfly system, where *H. melpomene* and *H. cydno*
are sister species that hybridize with *H. heurippa*, producing
introgression of wing pattern genes across species boundaries
(Mavárez et al., *Nature*, 2006).

*Step 1: Identify hybridization from gene trees.*

You have gene trees from 200 loci across 6 *Heliconius* species.
First, use ``quartet_network`` to test whether gene tree discordance
is due to ILS alone or also involves hybridization:

.. code-block:: shell

   phykit quartet_network -t gene_trees.nwk --json > quartets.json

Examine the output to see which quartets are classified as hybrid:

.. code-block:: shell

   # Quick summary
   python -c "
   import json
   data = json.load(open('quartets.json'))
   print(f'Tree-like: {data[\"tree_count\"]}')
   print(f'Hybrid: {data[\"hybrid_count\"]}')
   print(f'Unresolved: {data[\"unresolved_count\"]}')
   for q in data['quartets']:
       if q['classification'] == 'hybrid':
           print(f'  {q[\"dominant_topology\"]}  CFs: {q[\"cfs\"]}')
   "

Suppose the output shows that quartets involving *H. melpomene*
and *H. heurippa* are consistently classified as hybrid, with
asymmetric minor concordance factors — evidence of gene flow
between these lineages.

*Step 2a: Compute signal using quartet_network output directly.*

Feed the quartet JSON into ``network_signal`` along with a rooted
species tree and wing pattern measurements (e.g., forewing red band
area, log-transformed):

.. code-block:: shell

   phykit network_signal \
       -t species_tree.nwk \
       -d wing_pattern.tsv \
       --quartet-json quartets.json

The command automatically identifies the strongest hybrid signal
from the quartet classifications and estimates the inheritance
probability (gamma).

*Step 2b: Alternatively, specify hybrid edges explicitly.*

If you know the donor and recipient lineages (e.g., from prior
knowledge or external network inference), you can specify the
hybrid edge directly. Here, *H. melpomene* is the donor of wing
pattern alleles to *H. heurippa* with an estimated 25% introgression:

.. code-block:: shell

   phykit network_signal \
       -t species_tree.nwk \
       -d wing_pattern.tsv \
       --hybrid H_melpomene:H_heurippa:0.25

*Step 3: Interpret the results.*

Example output:

.. code-block:: text

   Hybrid edge: H_melpomene -> H_heurippa (gamma=0.2500)
   Network taxa: 6
   ---
   Blomberg's K: 0.6821    p-value: 0.0410
   Pagel's lambda: 0.5934    log-likelihood: -8.7231    p-value: 0.0085

Interpretation:

- **K = 0.68** (p = 0.04): significant phylogenetic signal, but less
  than expected under Brownian motion (K < 1). This is consistent with
  the wing pattern being phylogenetically conserved in most lineages but
  displaced in *H. heurippa* due to introgression from *H. melpomene*.
- **Lambda = 0.59** (p = 0.009): moderate phylogenetic signal. The
  trait is not evolving independently of the network (lambda > 0), but
  the fit is better with reduced covariance (lambda < 1).
- **Comparison with tree-based signal**: running ``phylogenetic_signal``
  on the species tree alone would likely produce a lower K value because
  the tree VCV does not account for the shared ancestry introduced by
  introgression. The network-based K is a more accurate estimate of how
  much evolutionary history explains trait variation.

*Why this matters*: Without accounting for the network, the tree
treats *H. heurippa*'s wing pattern as an independent observation. In
reality, its wing pattern was partly inherited from *H. melpomene*
through hybridization — the network VCV correctly reflects this shared
ancestry, producing a more accurate phylogenetic signal estimate.

|

**Algorithm — Bastide et al. 2018 network VCV**

For Brownian motion on a network, the trait at a hybrid node *h*
with parents *p1* and *p2* is:

``X_h = gamma * (X_p1 + noise_1) + (1-gamma) * (X_p2 + noise_2)``

The network VCV is computed recursively in topological order:

- Tree node *c* with parent *p*, edge length *l*: |br|
  ``V[c,c] = V[p,p] + l`` |br|
  ``V[c,j] = V[p,j]`` for all other nodes *j*
- Hybrid node *h* with parents *p1* (weight gamma) and *p2* (weight 1-gamma): |br|
  ``V[h,h] = gamma^2 * (V[p1,p1] + l1) + (1-gamma)^2 * (V[p2,p2] + l2) + 2*gamma*(1-gamma)*V[p1,p2]`` |br|
  ``V[h,j] = gamma * V[p1,j] + (1-gamma) * V[p2,j]``

The tip-by-tip submatrix is the VCV used for K and lambda.
