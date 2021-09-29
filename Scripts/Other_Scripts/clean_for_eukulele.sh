#!/bin/sh

cd /media/ink/EukMetaSanityAnnot/delmont/faa_for_eukulele

ls *.faa > faa.list

mkdir clean_faa
while read FILE; do
  awk '{gsub("\\.","*"); print}' $FILE > clean_faa/$FILE
done <faa.list
