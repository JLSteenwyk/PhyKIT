.. _cmd-alignment_recoding:
.. _command-alignment_recoding:

Alignment recoding
==================

Recode alignment into reduced alphabets

Command identity
----------------

:Canonical command: ``alignment_recoding``
:Handler: ``alignment_recoding``
:Aliases: aln_recoding, recode
:Standalone executables: pk_alignment_recoding, pk_aln_recoding, pk_recode
:Categories: Alignment & dataset utilities

Runtime interface
-----------------

.. include:: /_generated/commands/alignment_recoding.inc

Guidance, interpretation, and examples
--------------------------------------

Recode alignments using reduced character states.

Alignments can be recoded using established or
custom recoding schemes. Recoding schemes are
specified using the -c/--code argument. Custom
recoding schemes can be used and should be formatted
as a two column file wherein the first column is the
recoded character and the second column is the character
in the alignment.

.. code-block:: shell

	phykit alignment_recoding <fasta> [-c/--code <code>] [--json]

Codes for which recoding scheme to use: |br|
**RY-nucleotide** |br|
R = purines (i.e., A and G) |br|
Y = pyrimidines (i.e., T and C) |br|

**SandR-6** |br|
0 = A, P, S, and T |br|
1 = D, E, N, and G |br|
2 = Q, K, and R |br|
3 = M, I, V, and L |br|
4 = W and C |br|
5 = F, Y, and H |br|

**KGB-6** |br|
0 = A, G, P, and S |br|
1 = D, E, N, Q, H, K, R, and T |br|
2 = M, I, and L |br|
3 = W |br|
4 = F and Y |br|
5 = C and V |br|

**Dayhoff-6** |br|
0 = A, G, P, S, and T |br|
1 = D, E, N, and Q |br|
2 = H, K, and R |br|
3 = I, L, M, and V |br|
4 = F, W, and Y |br|
5 = C |br|

**Dayhoff-9** |br|
0 = D, E, H, N, and Q |br|
1 = I, L, M, and V |br|
2 = F and Y |br|
3 = A, S, and T |br|
4 = K and R |br|
5 = G |br|
6 = P |br|
7 = C |br|
8 = W |br|

**Dayhoff-12** |br|
0 = D, E, and Q |br|
1 = M, L, I, and V |br|
2 = F and Y |br|
3 = K, H, and R |br|
4 = G |br|
5 = A |br|
6 = P |br|
7 = S |br|
8 = T |br|
9 = N |br|
A = W |br|
B = C |br|

**Dayhoff-15** |br|
0 = D, E, and Q |br|
1 = M and L |br|
2 = I and V |br|
3 = F and Y |br|
4 = G |br|
5 = A |br|
6 = P |br|
7 = S |br|
8 = T |br|
9 = N |br|
A = K |br|
B = H |br|
C = R |br|
D = W |br|
E = C |br|

**Dayhoff-18** |br|
0 = F and Y |br|
1 = M and L |br|
2 = I |br|
3 = V |br|
4 = G |br|
5 = A |br|
6 = P |br|
7 = S |br|
8 = T |br|
9 = D |br|
A = E |br|
B = Q |br|
C = N |br|
D = K |br|
E = H |br|
F = R |br|
G = W |br|
H = C |br|

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*-c/--code*: argument to specify the recoding scheme to use |br|
*--json*: optional argument to print results as JSON
