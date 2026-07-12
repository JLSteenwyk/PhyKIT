.. _cmd-internode_labeler:
.. _command-internode_labeler:

Internode labeler
=================

Label internal nodes of a tree

Command identity
----------------

:Canonical command: ``internode_labeler``
:Handler: ``internode_labeler``
:Aliases: il
:Standalone executables: pk_internode_labeler, pk_il
:Categories: Tree manipulation & utilities

Runtime interface
-----------------

.. include:: /_generated/commands/internode_labeler.inc

Guidance, interpretation, and examples
--------------------------------------

Appends numerical identifiers to bipartitions in place of support values.
This is helpful for pointing to specific internodes in supplementary files
or otherwise.  

.. code-block:: shell

   phykit internode_labeler <tree> [-o/--output <file>] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-o/--output*: optional argument to name the outputted tree file |br|
*--json*: optional argument to print results as JSON
