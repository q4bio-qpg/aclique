AClique
=======

Aclique.c is the main product of interest here.

Build with

    gcc -g -O3 aclique.c -lm -o aclique

It takes an unweighted or weighted clique in DIMACS format as input,
but works best on unweighted cliques. Use `aclique -h` for help.

An example usage:

    ./aclique gen400_p0.9_75.clq -s 1 -t 10

There are more clique examples in datasets-*.tar.gz.
These are weighted cliques, but they can be converted to unweighted
using the DIMACS-w2u.pl script.

Evaluation tools
================

These are in the evaluation subdirectory and primarily exist as
supplementary material for the paper.  They are not suitable for
production work.

The preprocessing directory also contains tools for some experiments
used in the paper, such as adding noise to explore the impact on
parsimony and maximum likelhood methods vs maximum compatibility.

Covid pipeline
--------------

The covid data was locally compiled from GISAID and consisted of
1,549,263 genomes.

Our file format was one line per genome, white space separated, with fields
Name Lineage Date Accession-number md5sum MD5Sum Sequence

Producing a subset
------------------

We can select random subsets using "shuf | head". Eg:

    xz -d < sars_1549263_sequences_assigned.tsv.xz \
      | shuf | head -100 \
      | awk '{printf(">%04d#%s\t%s\n%s\n",NR,$2,$2,$NF)}' > shuf100.fa

A version of this standard utility using a fixed seed as a shell function is:

    shuf() {
        perl -e '$count=shift(@ARGV);
                 srand(shift(@ARGV));
                 while (<>) {push(@lines, $_)}
                 for ($i=0;$i<$count;) {
                     $l = int(rand $#lines+.999);
                     next unless defined($lines[$l]);
                     print $lines[$l];
                     undef $lines[$l];
                     $i++;
                 }' $1 $2
    }


Creating a true tree to compare against
---------------------------------------

It has sequence names along with Pangolin assigned LINEAGE.  The
lineage is of the style A.B.C is a child of A.B which is a child of A.
This is essentially a tree.

We can convert a fasta file with name + lineage to a Newick tree using
"lineage_fa2tree.pl".

    evaluation/lineage_fa2tree.pl shuf100.fa > shuf100.nwk


Align sequences
---------------

We could use MAFFT here to do a better job, but a crude strategy of
getting all the data aligned is to align each to a common reference.

    covid=/nfs/srpipe_references/references/SARS-CoV-2/default/all/fasta/MN908947.3.fa
    preprocessing/align_seqs.pl $covid shuf100.fa > shuf100.aln

This aligns using minimap2 and then processes the SAM output POS,
CIGAR and SEQ fields to create reference-coordinate sequences, padded
with Ns.

It's not perfect and no filtering is done currently.

A MAFFT example:

mafft --auto unaligned.fa | \
  perl -lane 'if (/^>/) {print;next} tr/acgt-/ACGTN/;tr/ACGTN/N/c;print' \
  > mafft.fa

Also see preprocessing/align_seqs_mafft.pl.


Build trees
-----------

Two methods are compat (Maximum Compatibility) and iqtree3's own
method (Maximum likelihood?).

Use compat.sh which sets PYTHONPATH and runs the compiled compat.

    evaluation/compat.sh shuf100.aln shuf100.aln.compat.tree

(Also see compat options such as reducing number of trees to iterate over.)

Use ~jkb/bin/iqtree3 to build its own tree.

    iqtree3 --redo -s shuf100.aln

This produces shuf100.aln.treefile


Comparing trees
---------------

iqtree3 -rf computes a Robinson-Foulds edit distance between a
pair of trees.  It doesn't use branch lengths, but does tell us about
topology changes.

    iqtree3 -rf shuf100.nwk shuf100.aln.compat.tree
    cat shuf100.aln.compat.tree.rfdist

    iqtree3 -rf shuf100.nwk shuf100.aln.treefile
    cat shuf100.aln.treefile.rfdist

The evaluation/treedist.r also uses the TreeDist R package to compare
trees.

Note that tree comparison scores are problematic and depend heavily on
the N-way nature of the tree construction.  A binary tree may look bad
compared to a system that puts many children into a single node.
These tree comparisons work on shape and do not take distances into account.


Scoring trees
-------------

Instead of comparing trees to a known lineage, we can also score trees
based on the lineage clustering.  This is more flexible and can take
distance into account.

A variety of scripts are available in evaluation/tree_score*.pl.

They work by sorting the child nodes by the dominant lineage present
in the children, and then performing a left-to-right tree walk.  We
can then count lineage switches as we walk.  A perfectly clustered
tree should then have N-1 switches for N lineages.  This is less
susceptible to differences in binary vs non-binary tree construction.

Variants of these scripts focus purely on lineage sorting (*lin.pl)
and also by branch lengths (*len.pl).


Links
=====

- https://academic.oup.com/ve/article/7/2/veab064/6315289?login=false

  Pangolin paper.  This is the tool that assigned the covid lineages
  we are taking as the truth set.


- https://pmc.ncbi.nlm.nih.gov/articles/PMC12380668

  A paper that uses iqtree3 on covid data.
  Has useful hints on alignment methods and options to iqtree3.  Also
  gives credence to our methodology as we can cite a previous paper.

  Contains in supplementary material:

  #!/bin/bash
  # Declare variables
  INPUT_FASTA="gisaid_hcov-19_2025_04_10_16_Goias.fasta"
  ALN_FASTA="aln_sarscov2_go_2025_04_16.fasta"
  MODEL="MFP"
  BOOTSTRAP=1000
  ALRT=1000
  THREADS=100

  # Multiple sequence alignment
  mafft --auto "$INPUT_FASTA" > "$ALN_FASTA"

  # Phylogenetic tree analyses
  iqtree3 -s "$ALN_FASTA" -m "$MODEL" -bb "$BOOTSTRAP" -alrt "$ALRT" -nt "$THREADS"
  