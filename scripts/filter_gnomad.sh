#!/bin/bash

module load samtools
module load GATK

bcftools view -f PASS -Q 0.01 /fdb/gnomad/release-171003/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz | \
	/home/mcgaugheyd/git/eye_var_Pathogenicity/scripts/vcf_add_fake_sample.py | bgzip > /data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/temp.vcf.bgz
tabix -p vcf /data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/temp.vcf.bgz
GATK -m 8g SelectVariants \
	-R /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta \
	-V /data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/temp.vcf.bgz \
	-o /data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/cgnomad.exomes.r2.0.2.sites.maxAF_0.01_7percent.vcf.gz \
	-L /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/converted_exome_bait_beds/TruSeq1.0.b37.merge.bed \
	-fraction 0.07 \
	--interval_padding 20
