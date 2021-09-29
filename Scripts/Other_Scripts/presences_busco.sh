#!/bin/sh

cd /media/ink/EukMetaSanityAnnot/busco_euks/

find -name "*.faa" | grep -E "single_copy_busco_sequences|derep_busco_sequences" | awk 'gsub("/","\t")' | awk '{print $2"\t"$6}' > ../busco_presences.txt
