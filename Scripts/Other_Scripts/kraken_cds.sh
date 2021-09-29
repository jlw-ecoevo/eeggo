#!/bin/sh

cd /media/ink/MMETSP/cds/

find -type f -name "*.cds" | awk 'gsub("./","")' > cds.list

while read FILE; do
  NAME=${FILE%.cds}
  echo $NAME
  kraken2 --db /media/ink/kraken2/ntfast --use-names --threads 50 --output ${NAME}.kraken --report ${NAME}.krakenreport $FILE
done <cds.list
