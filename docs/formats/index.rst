.. _input-formats:

Input formats and shared conventions
====================================

The canonical command page is authoritative when a command has stricter
requirements than the shared conventions below.

Multiple sequence alignments
----------------------------

PhyKIT detects alignment format from file content. Supported formats are FASTA,
Clustal, MAF, Mauve, PHYLIP (interleaved, sequential, and relaxed), and
Stockholm. Aligned sequences must have equal lengths.

FASTA records use the first whitespace-delimited token after ``>`` as the taxon
identifier::

   >species_A optional description
   ACGTACGT
   >species_B
   ACGT-CGT

Use unique identifiers. Gap characters are normally ``-``; the handling of
ambiguous characters depends on the statistic being calculated.

Phylogenetic trees
------------------

Tree inputs use Newick unless a command page states otherwise::

   ((species_A:0.1,species_B:0.1):0.2,species_C:0.3);

Terminate each tree with a semicolon. Branch lengths and internal support values
are required only by commands that use them. Some multi-tree commands accept
either one Newick tree per line or one tree-file path per line.

Taxon names
-----------

Taxon matching is exact and case-sensitive. ``species_A``, ``Species_A``, and
``species A`` are different identifiers. Use the same spelling in alignments,
tree tips, trait tables, grouping files, and mapping files. Avoid whitespace in
identifiers when possible.

Single-trait tables
-------------------

Commands such as ``phylogenetic_signal`` accept two tab-separated columns with
no header. Blank lines and lines beginning with ``#`` are ignored::

   # taxon<TAB>continuous_value
   species_A	1.25
   species_B	0.80

Multi-trait tables
------------------

Regression, ordination, and other multi-trait commands use a header. The first
column contains taxa and subsequent columns contain named traits::

   taxon	body_mass	brain_size
   species_A	1.25	0.90
   species_B	0.80	0.72

Use literal tab characters, not aligned runs of spaces. Numeric analyses reject
non-numeric trait values unless the command explicitly supports categorical data.

Discrete traits and grouping files
-----------------------------------

Discrete-trait commands commonly use a header so ``--trait`` can identify the
selected column::

   taxon	diet
   species_A	carnivore
   species_B	herbivore

Other two-column files, such as rename maps, usually contain an existing name and
its replacement. Consult the command page for column order and header rules.

Missing values and missing taxa
-------------------------------

PhyKIT does not define one missing-value token for every command. Do not assume
that ``NA``, ``NaN``, an empty field, or ``?`` is accepted. Check the command's
input section and runtime interface. Missing taxa are also command-specific:
some commands use shared taxa, some expose ``--missing-taxa``, and others reject
nonmatching sets.

Paths and list files
--------------------

Paths passed directly on the command line are resolved from the current working
directory. For files containing lists of other file paths, consult the command
page to determine whether relative paths are resolved from the working directory
or the list file's directory. Absolute paths remove that ambiguity.
