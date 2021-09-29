#!/bin/sh

cd /media/ink/EukMetaSanityAnnot/orthologs/

fasttree busco_orthologs.fasta.trim  > busco_tree_pretrim.newick
