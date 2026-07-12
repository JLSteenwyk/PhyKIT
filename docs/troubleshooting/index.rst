.. _troubleshooting:

Troubleshooting
===============

Start with the command's canonical reference page and its live interface::

   phykit <command> --help

Installation and command discovery
----------------------------------

Confirm that PhyKIT and the active Python interpreter come from the same
environment::

   python -m pip show phykit
   command -v phykit
   phykit version

If ``phykit`` is not found, activate the environment used for installation or
run ``python -m pip install --upgrade phykit`` in the intended environment.
PhyKIT requires Python 3.10 or newer.

Input and parsing errors
------------------------

- Resolve relative paths from the current working directory. Use absolute paths
  when an input list refers to files outside that directory.
- Ensure FASTA records have unique identifiers and aligned sequences have equal
  lengths.
- Ensure Newick trees end with a semicolon and every required tip has a name.
- Match taxon names exactly across files. Capitalization, spaces, and underscores
  are significant.
- Use tab characters, rather than runs of spaces, in documented TSV inputs.
- Preserve the documented header or no-header convention for each command.

Unexpected or empty results
---------------------------

Check filters, thresholds, missing-data rules, and taxon intersections before
changing the analysis. An empty result can be valid when no site, split, taxon,
or model satisfies the selected criteria. Use ``--json`` where available to
distinguish an empty structured result from formatting intended for people.

Reproducibility
---------------

Record the output of ``phykit version``, the complete command, input checksums,
and random seed. Commands that expose ``--seed`` should use an explicit value in
reproducible workflows. Do not assume defaults are unchanged across releases;
consult the versioned command catalog.

Reporting a problem
-------------------

Open a `GitHub issue <https://github.com/JLSteenwyk/PhyKIT/issues>`_ with the
PhyKIT and Python versions, operating system, complete command, full error text,
and a minimal input that reproduces the problem. Remove sensitive data before
attaching files.
