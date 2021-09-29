#!/bin/sh

cd /media/ink/EukMetaSanityAnnot/orthologs/

ls *.faa > faa_files.txt

while read FILE; do
  echo $FILE
  muscle -in $FILE -out ${FILE}.aln &
done <faa_files.txt
