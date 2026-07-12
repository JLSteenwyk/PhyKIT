.. _glossary:

Glossary
========

Alias
   An alternative, shorter command name. For example, ``aln_len`` is an alias
   for ``alignment_length`` and dispatches to the same implementation.

Bipartition
   A split of the taxa induced by removing an edge from a phylogenetic tree.

Canonical command
   The stable, descriptive command name used for its reference page and primary
   machine-readable identity.

Cladogram
   A tree drawing that emphasizes topology and does not interpret branch lengths
   as evolutionary distance or time.

Gene tree discordance
   Differences among gene-tree topologies or between gene trees and a species tree.

Handler
   The Python function that implements a canonical command and its aliases.

Newick
   A parenthetical plain-text representation of a phylogenetic tree.

Patristic distance
   The sum of branch lengths along the path between two tree tips.

Phylogenetic signal
   Statistical dependence between trait similarity and phylogenetic relatedness.
   Individual statistics such as Blomberg's K and Pagel's lambda have distinct
   definitions and interpretations.

Phylogram
   A tree drawing in which branch lengths represent an estimated amount of change.

Polytomy
   An internal node with more than two immediate descendants, often representing
   unresolved relationships.

Standalone executable
   A ``pk_*`` command installed for a canonical command or alias, such as
   ``pk_alignment_length`` or ``pk_aln_len``.

Taxon
   A named unit represented by an alignment record, tree tip, or trait-table row.

TSV
   A tab-separated values file. Tabs delimit columns; spaces are not equivalent.

VCV matrix
   A phylogenetic variance-covariance matrix describing expected covariance among
   taxa under an evolutionary model.
