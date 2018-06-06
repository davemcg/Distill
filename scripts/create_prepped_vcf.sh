#!/bin/bash

input=$1
output=$2

# takes unsorted, blank vcf
# and outputs
# a vcf that is sorted and has a fake sample
/home/mcgaugheyd/git/eye_var_Pathogenicity/scripts/vcf_add_fake_sample.py $input | \
	awk '$1 ~ /^#/ {print $0;next} {print $0 | "LC_ALL=C sort -k1,1 -k2,2n"}' | bgzip > $output
tabix -p vcf $output
