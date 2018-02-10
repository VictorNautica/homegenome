#!/bin/bash

awk '$2!~/Treatment/ && $2!~/Medication/' UKBB_full_manifest.tsv > phenolist_nomedications.tsv
grep -v -E '(Leg|Arm) fat' phenolist_nomedications.tsv |\
	grep -v 'Impedance' |\
	grep -v 'Home area' |\
	grep -v 'reported: fracture' |\
	grep -v -E 'Illnesses of' > phenolist_final.tsv

echo "Number of lines, initial:"
wc -l UKBB_full_manifest.tsv
echo "Number of lines, remaining:"
wc -l phenolist_final.tsv
