#!/bin/sh

cd /media/ink/EukMetaSanityAnnot/

busco -c 60 -i pep/ -l eukaryota_odb10 -o busco_euks -m proteins
