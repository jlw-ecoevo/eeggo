#!/bin/sh

cd /media/ink/EukMetaSanityAnnot/delmont/

find -name "*.1.Tier.faa" | awk '{gsub(".1.Tier.faa",""); print}' > eukms.txt

while read NAME;do
  NAMEBASE=`basename $NAME` 
  blastp -db /media/ink/riboblastdb/riboprot -query ${NAME}.1.Tier.faa -outfmt 6 -out ${NAME}.riboblast -num_threads 10
  awk '$11<1e-10' ${NAME}.riboblast | awk '{print $1}' | sort | uniq > ${NAME}.riboprot 
  gffread -x ${NAME}.cds -g Contigs/${NAMEBASE}.fa ${NAME}.1.Tier.gff3 
done <eukms.txt
