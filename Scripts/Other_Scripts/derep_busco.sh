#!/bin/sh

cd /media/ink/EukMetaSanityAnnot/busco_euks/

find -name "*.faa" | grep "multi_copy_busco_sequences" > multi_files.txt
awk 'gsub("multi_copy","derep")' multi_files.txt > derep.txt
paste multi_files.txt derep.txt > to_derep.txt

while read FILE OUTPUT; do
  echo $FILE
  echo $OUTPUT
  mkdir `dirname $OUTPUT`
  awk '/^>/{if(N)exit;++N;} {print;}' $FILE > $OUTPUT
done <to_derep.txt

