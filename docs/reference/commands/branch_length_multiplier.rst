.. _cmd-branch_length_multiplier:
.. _command-branch_length_multiplier:

Branch length multiplier
========================

Multiply branch lengths by a factor

Command identity
----------------

:Canonical command: ``branch_length_multiplier``
:Handler: ``branch_length_multiplier``
:Aliases: blm
:Standalone executables: pk_branch_length_multiplier, pk_blm
:Categories: Tree manipulation & utilities

Runtime interface
-----------------

.. include:: /_generated/commands/branch_length_multiplier.inc

Guidance, interpretation, and examples
--------------------------------------

Multiply branch lengths in a phylogeny by a given factor.
                
This can help modify reference trees when conducting simulations
or other analyses.  

.. code-block:: shell

   phykit branch_length_multiplier <tree> -f n [-o/--output <output_file>] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-f/--factor*: factor to multiply branch lengths by |br|
*-o/--output*: optional argument to name the outputted tree file. Default 
output will have the same name as the input file but with the suffix ".factor_(n).tre" |br|
*--json*: optional argument to print results as JSON
