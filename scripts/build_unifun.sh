#!/bin/bash

cd /data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/unifun
wget https://ndownloader.figshare.com/articles/4704955/versions/1
unzip 1

/home/mcgaugheyd/git/eye_var_Pathogenicity/scripts/vcf_add_fake_sample.py UniFun_deleterious.vcf | bgzip > UniFun_deleterious.fakeS.vcf.gz
/home/mcgaugheyd/git/eye_var_Pathogenicity/scripts/vcf_add_fake_sample.py UniFun_benign.vcf | bgzip > UniFun_benign.fakeS.vcf.gz

tabix -p vcf UniFun_deleterious.fakeS.vcf.gz
tabix -p vcf UniFun_benign.fakeS.vcf.gz

bash ~/git/variant_prioritization/Snakemake.wrapper.sh ~/git/eye_var_Pathogenicity/config_variant_prioritization__unifunD.yaml
bash ~/git/variant_prioritization/Snakemake.wrapper.sh ~/git/eye_var_Pathogenicity/config_variant_prioritization__unifunB.yaml


module load gemini
time gemini query --header -q "SELECT * from variants" UniFun_deleterious.fakeS.PED_faux.gemini.db | bgzip > unifun_deleterious.gemini.tsv.gz
time gemini query --header -q "SELECT * from variants" UniFun_benign.fakeS.PED_faux.gemini.db | bgzip > unifun_benign.gemini.tsv.gz