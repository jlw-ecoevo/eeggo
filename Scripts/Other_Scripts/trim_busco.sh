#!/bin/sh

cd /media/ink/EukMetaSanityAnnot/orthologs/

ls *.faa > faa_files.txt

while read FILE; do
  echo $FILE
  trimal -automated1 -in ${FILE}.aln -out ${FILE}.aln.trim
done <faa_files.txt

while read FILE; do
  awk '{gsub(" .*",""); print}' ${FILE}.aln.trim > ${FILE}.aln.trim.clean
done <faa_files.txt
