
#!/bin/sh

cd /media/ink/EukMetaSanityAnnot/alexander/

find -name "*.faa.gz" | awk '{gsub(".faa.gz",""); print}' > eukms.txt

while read NAME;do
  gunzip ${NAME}.faa.gz
  gunzip ${NAME}.fna.gz
  gunzip ${NAME}.nr.gff3.gz
  blastp -db /media/ink/riboblastdb/riboprot -query ${NAME}.faa -outfmt 6 -out ${NAME}.riboblast -num_threads 10
  awk '$11<1e-10' ${NAME}.riboblast | awk '{print $1}' | sort | uniq > ${NAME}.riboprot 
  gffread -x ${NAME}.cds -g ${NAME}.fna ${NAME}.nr.gff3 
done <eukms.txt
