.. _cmd-kf_distance:
.. _command-kf_distance:

Kuhner-Felsenstein distance
===========================

Branch score distance between trees (topology + branch lengths)

Command identity
----------------

:Canonical command: ``kf_distance``
:Handler: ``kf_distance``
:Aliases: kf, kf_dist, kuhner_felsenstein_distance
:Standalone executables: pk_kf_distance, pk_kf, pk_kf_dist, pk_kuhner_felsenstein_distance
:Categories: Tree comparison & consensus

Runtime interface
-----------------

.. include:: /_generated/commands/kf_distance.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate the Kuhner-Felsenstein (KF) branch score distance between two trees.

Unlike Robinson-Foulds distance which only considers topology, KF distance
incorporates both topology and branch length differences. The KF distance
is calculated as KF = sqrt(sum over all splits of (b1_i - b2_i)^2), where
b1_i and b2_i are branch lengths for split i in each tree. Splits absent
from one tree use branch length 0. Both internal and terminal branches are
included.

PhyKIT will print out
col 1: the plain KF distance and
col 2: the normalized KF distance.

KF distances are calculated following Kuhner & Felsenstein, Journal of
Computational Biology (1994), doi: 10.1089/cmb.1994.1.183.

Cross-validated against R's phangorn::KF.dist().

.. code-block:: shell

	phykit kf_distance <tree_file_zero> <tree_file_one> [--json]

Options: |br|
*<tree_file_zero>*: first argument after function name should be a tree file |br|
*<tree_file_one>*: second argument after function name should be a tree file |br|
*--json*: optional argument to print results as JSON

**R validation:** Validated against ``phangorn`` in R
(see ``tests/r_validation/validate_kf_distance.R``).
