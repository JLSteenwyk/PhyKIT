.. _cmd-collapse_branches:
.. _command-collapse_branches:

Collapse bipartitions
=====================

Collapse low-support bipartitions

Command identity
----------------

:Canonical command: ``collapse_branches``
:Handler: ``collapse_branches``
:Aliases: cb, collapse
:Standalone executables: pk_collapse_branches, pk_cb, pk_collapse
:Categories: Tree manipulation & utilities

Runtime interface
-----------------

.. include:: /_generated/commands/collapse_branches.inc

Guidance, interpretation, and examples
--------------------------------------

Collapse branches on a phylogeny according to bipartition support.

Bipartitions will be collapsed if they are less than the user specified
value.    

.. code-block:: shell

   phykit collapse_branches <tree> -s/--support n [-o/--output <output_file>] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-s/--support*: bipartitions with support less than this value will be 
collapsed |br|
*-o/--output*: optional argument to name the outputted tree file. Default 
output will have the same name as the input file but with the suffix 
".collapsed_(support).tre" |br|
*--json*: optional argument to print results as JSON
