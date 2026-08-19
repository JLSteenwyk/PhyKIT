.. _tutorial-02:
.. _tutorial-evaluating-gene-gene-covariation:

Tutorial 2: Evaluating gene-gene covariation
============================================

Objectives
----------

- Complete the evaluating gene-gene covariation workflow.
- Interpret the reported values and generated artifacts in their scientific context.
- Identify the canonical command references for each analysis step.

Prerequisites and working directory
-----------------------------------

Install the current PhyKIT release and create a dedicated working directory.
Download the data linked in this tutorial into that directory before running
the commands. All paths below are relative to this directory.

.. code-block:: shell

   mkdir phykit-tutorial-02
   cd phykit-tutorial-02

Related command references
--------------------------

- :doc:`Covarying Evolutionary Rates </reference/commands/covarying_evolutionary_rates>`
- :doc:`Episodic Rate Covariation </reference/commands/episodic_rate_covariation>`
- :doc:`Projected Covarying Rates </reference/commands/projected_covarying_rates>`
- :doc:`Root Tree </reference/commands/root_tree>`

Workflow
--------

Identifying genes that significantly covary (or coevolve) with one another is known to accurately and sensitively 
identify genes that have shared functions, are coexpressed, and/or are part of the same multimeric complexes 
(`Sato et al. 2005 <https://academic.oup.com/bioinformatics/article/21/17/3482/212654>`_; 
`Clark et al. 2012 <https://genome.cshlp.org/content/22/4/714.full>`_).
Furthermore, gene-gene covariation serves as a powerful evolution-based genetic screen for predicting gene function
(`Brunette et al. 2019 <https://www.pnas.org/content/116/39/19593>`_).

PhyKIT implements the CovER branch-ratio method to identify genes that covary with one another. For every branch,
the gene-tree length is divided by the corresponding reference-tree length, and the resulting branch-rate vectors
are correlated. Thus, the two gene trees and reference tree must have the same rooted topology after all three are
pruned to their shared taxa. Gene trees should be inferred with the reference topology constrained rather than
estimated as unconstrained topologies. See the :doc:`command reference
</reference/commands/covarying_evolutionary_rates>` for filtering and statistical details.

To provide a comprehensive tutorial, we will start with the sequence alignments for three genes and their constrained 
tree topologies that match the putative species tree from `Shen et al. 2020
<https://www.biorxiv.org/content/10.1101/2020.05.11.088658v1.abstract>`_. 

Choosing a method when topologies differ
-----------------------------------------

The published CovER implementation, ``cover``, is appropriate for the
constrained and topology-matched trees used in this tutorial. If gene trees
were inferred without a topological constraint and disagree with one another
or with the species tree, the experimental ``pcover`` command provides a
different analysis:

.. code-block:: shell

   phykit pcover unconstrained_gene_a.tre unconstrained_gene_b.tre \
       -r species_tree.tre --json > projected_rates.json

``pcover`` projects each gene's pairwise patristic distances onto the species
tree's unrooted edge basis before correlating relative rates. It reports a
projection NRMSE for each gene; zero is an exact fit to the reference basis,
whereas larger values indicate unexplained topological distance signal.
Interpret the correlation only alongside both NRMSE values. The ``cover`` and
``pcover`` coefficients estimate different quantities and should not be mixed
within the same gene-pair network. See :doc:`Projected Covarying Rates
</reference/commands/projected_covarying_rates>` for assumptions and output
details.

Localizing covariation after projection
---------------------------------------

A whole-tree ``pcover`` coefficient asks whether projected rates covary across
all retained reference edges. To investigate whether that signal is
concentrated in a particular part of the tree, run the experimental episodic
scan on the same inputs and settings:

.. code-block:: shell

   phykit erc_scan unconstrained_gene_a.tre unconstrained_gene_b.tre \
       -r species_tree.tre --permutations 9999 --seed 814 \
       --output episodic_clades.tsv \
       --annotated-tree episodic_reference.tre --json > episodic_scan.json

For each eligible rooted clade, the scan contrasts the products of the two
standardized projected rates inside the clade with those elsewhere in the
tree. It permutes rates among similarly deep edges and reports both
unadjusted and maximum-statistic family-wise error rate (FWER) p-values. Base
discovery claims on ``fwer_p_value``, inspect both projection NRMSE values, and
record the permutation count and seed. The scan is exploratory and does not
replace independent validation of a proposed episodic interaction. See
:doc:`Episodic Rate Covariation
</reference/commands/episodic_rate_covariation>` for the statistic, edge
eligibility rules, and interpretation of concordant versus antagonistic
coupling.

.. centered::
   Download test data:
   :download:`gene_gene_covariation_tutorial.tar.gz </data/gene_gene_covariation_tutorial.tar.gz>`

|

Step 0: Prepare data
********************
The mirror tree method for determining significant gene-gene covariation requires that both input phylogenies have the same topology.
As a result, gene trees must be constrained to the species tree, which is typically inferred from whole genome or proteome data.
In the present tutorial, the species tree has already been inferred. Additionally, the guide trees used to constrain the gene trees
have been generated. These trees were generated by pruning the species tree to match the taxon representation of the sequences in
the multiple sequence alignment.

Step 1: Estimate gene tree branch lengths
*****************************************
To infer the constrained tree, we will use `IQ-TREE2 <http://www.iqtree.org/>`_. The species tree (or guide tree) is specified
with the *-g* argument. Lastly, the best-fitting substitution model was specified according to what was reported in 
`Shen et al. 2020 supplementary data <https://jlsteenwyk.com/publication_pdfs/2020_Shen_etal_Science_Advances.pdf>`_; however,
if the best-fitting model is unknown, this will have to be determined prior to estimating gene tree branch lengths.

Estimate the gene tree branch lengths using the following commands:

.. code-block:: shell

   # infer constrained trees
   iqtree2 -s Shen_etal_SciAdv_2020_NDC80.fa -te Shen_etal_SciAdv_2020_NDC80.constrained.tre -pre Shen_etal_SciAdv_2020_NDC80 -m JTT+G4+F -keep-ident
   iqtree2 -s Shen_etal_SciAdv_2020_NUF2.fa -te Shen_etal_SciAdv_2020_NUF2.constrained.tre -pre Shen_etal_SciAdv_2020_NUF2 -m LG+G4 -keep-ident
   iqtree2 -s Shen_etal_SciAdv_2020_SEC7.fa -te Shen_etal_SciAdv_2020_SEC7.constrained.tre -pre Shen_etal_SciAdv_2020_SEC7 -m LG+G4 -keep-ident

Step 2: Root constrained trees
******************************
To ensure PhyKIT traverses each tree the same, root each tree using the outgroup taxa. PhyKIT has
a function for rooting and takes as input a single column file with the names of the outgroup taxa.
For sake of simplicity, I have provided the necessary input files.

Root the trees with the inferred branch lengths using the following commands:

.. code-block:: shell

   # root trees
   pk_root Shen_etal_SciAdv_2020_NUF2.treefile -r Shen_etal_SciAdv_2020_NUF2_taxa_for_rooting.txt -o Shen_etal_SciAdv_2020_NUF2.treefile.rooted
   pk_root Shen_etal_SciAdv_2020_SEC7.treefile -r Shen_etal_SciAdv_2020_SEC7_taxa_for_rooting.txt -o Shen_etal_SciAdv_2020_SEC7.treefile.rooted
   pk_root Shen_etal_SciAdv_2020_NDC80.treefile -r Shen_etal_SciAdv_2020_NDC80_taxa_for_rooting.txt -o Shen_etal_SciAdv_2020_NDC80.treefile.rooted

Step 3: Evaluate gene-gene covariation
**************************************
When determining gene-gene covariation, it is best to use a high significance threshold for coevolutionary coefficients.
Here, we will use a threshold of 0.5; however, I recommend users explore their own data and distribution of coevolutionary 
coefficients.

To evaluate gene-gene covariation, execute the following commands:

.. code-block:: shell

   # Evaluate gene-gene covariation between NUF2 and SEC7
   phykit cover Shen_etal_SciAdv_2020_NUF2.treefile.rooted Shen_etal_SciAdv_2020_SEC7.treefile.rooted -r Shen_etal_SciAdv_2020_species_tree.tre

   # Evaluate gene-gene covariation between NDC80 and SEC7
   phykit cover Shen_etal_SciAdv_2020_NDC80.treefile.rooted Shen_etal_SciAdv_2020_SEC7.treefile.rooted -r Shen_etal_SciAdv_2020_species_tree.tre

.. code-block:: text

   0.1885  0.0
   0.2105   0.0

Given our thresholds, neither *NUF2* nor *NDC80* significantly covary with *SEC7*. Next, evaluate gene-gene covariation between
*NUF2* and *NDC80*.

.. code-block:: shell

   # Evaluate gene-gene covariation between NUF2 and NDC80
   phykit cover Shen_etal_SciAdv_2020_NUF2.treefile.rooted Shen_etal_SciAdv_2020_NDC80.treefile.rooted -r Shen_etal_SciAdv_2020_species_tree.tre

.. code-block:: text

   0.6693  0.0

These two genes significantly covary with one another. This raises the hypothesis that these two genes have shared function. A literature-
based examination of these genes reveals the encoded proteins are part of the same kinetochore-associated complex termed the 
`NDC80 complex <https://www.yeastgenome.org/complex/CPX-548>`_. Thus, PhyKIT is useful for determining gene-gene covariation, which can be 
driven by shared function, coexpression, and/or membership in the same multimeric complexes.

|

Expected artifacts
------------------

Each step identifies its expected terminal output or generated files. Confirm
that those artifacts exist before continuing to the next step; filenames are
relative to the tutorial working directory unless an absolute path is shown.

Troubleshooting
---------------

- Run ``phykit <command> --help`` to compare an invocation with the live interface.
- Confirm that downloaded files are in the current working directory and retain
  the filenames shown in the tutorial.
- For parsing errors, compare taxon names exactly across alignments, trees, and
  trait tables, including capitalization and underscores.
- See :doc:`Troubleshooting </troubleshooting/index>` for installation, format,
  and error-reporting guidance.
