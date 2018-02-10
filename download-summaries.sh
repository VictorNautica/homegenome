#!/bin/bash

set -e

awk -F '\t' 'NR>6{print $4}' phenolist_final.tsv > wget_commands.txt
i=0

while read l; do
	i=$((i+1))
	echo "command: $l"
	eval $l
	f=${l##* }
	nf=${f%.gz}
	echo "unzipping file $f, number $i"
	gunzip -c $f | awk '$9<1e-5' > $nf
	rm $f
done < wget_commands.txt

