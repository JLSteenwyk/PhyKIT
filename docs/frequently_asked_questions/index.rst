.. _faq:

Frequently asked questions
==========================

How do I choose a command?
--------------------------

Start with the :doc:`task-oriented command reference </reference/index>`. Each
canonical page lists aliases, accepted arguments, defaults, output behavior,
scientific interpretation, examples, limitations, validation, and citations.

How do I see the options for my installed version?
--------------------------------------------------

Run ``phykit <command> --help``. The website documents the current release, while
the local help reflects the version active in your environment.

Can PhyKIT produce machine-readable output?
-------------------------------------------

Most analytical commands expose ``--json``. Check the generated runtime table on
the canonical command page. The full command and argument catalog is available as
:download:`commands.json </_static/commands.json>`.

Why are taxa being dropped or reported as missing?
---------------------------------------------------

Taxon names must match exactly across inputs. Commands differ in whether they
reject mismatches, use the shared intersection, or allow partial taxon sets. Read
the command's missing-taxa option and the :doc:`format conventions </formats/index>`.

What alignment and tree formats are supported?
----------------------------------------------

See :doc:`Input formats and shared conventions </formats/index>`. Alignment format
is detected from content; tree inputs are normally Newick.

How do I make an analysis reproducible?
---------------------------------------

Record the PhyKIT version, complete command, input checksums, and explicit random
seed. Keep structured JSON output when available and cite the methods identified
on the command page.

How do I report a bug or request a feature?
-------------------------------------------

Open a `GitHub issue <https://github.com/JLSteenwyk/PhyKIT/issues>`_. For bugs,
include a minimal reproducible example, complete error text, PhyKIT and Python
versions, and operating system. See :doc:`Troubleshooting </troubleshooting/index>`
before posting installation or input-format problems.
