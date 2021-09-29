#!/bin/sh

cd /media/ink/EukMetaSanityAnnot/orthologs/

/media/ink/catfasta2phyml/catfasta2phyml.pl -f -c *.faa.aln.trim.clean > busco_orthologs.fasta
trimal -gappyout -in busco_orthologs.fasta -out busco_orthologs.fasta.trim
