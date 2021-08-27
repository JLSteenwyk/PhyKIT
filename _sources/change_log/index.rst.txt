.. _change_log:


Change log
==========

^^^^^

Major changes to PhyKIT are summarized here.

**1.2.2**:a dded function to get the subtree of the last common ancestor among a set of taxa

**1.2.0**: added command line interfaces for all functions so that each command 
can easily be executed. For example, 'phykit aln_len -h' can now be
called using 'pk_aln_len -h'

**1.1.0**: added faidx (alias: get_entry; ge) function to extract fasta entries from a
multi-fasta file

**1.0.3**: added rooting procedure before calculating RF to handle comparing unrooted
and rooted trees

**1.0.2**: function that calculates Robinson Foulds distance (robinson_foulds_distance;
rf_distance; rf_dist; rf) now can take trees that differ in topology. PhyKIT
will first determine shared tips between the two trees and prune both trees
to a common set of tips. Next, PhyKIT will calculate the Robinson Foulds 
distance.

**0.1.3**: Added function (column_score; cs) to calculate the quality of
an alignment given an input query alignment and a reference
alignment to compare it to

**0.1.2**: Added function (sum_of_pairs_score; sops; sop) to calculate
the quality of an alignment given an input query alignment
and a reference alignment to compare it to

**0.0.9**: PhyKIT now handles error stemming from piping output