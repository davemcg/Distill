# run in /data/mcgaugheyd/projects/nei/hufnagel/colombia_cohort
# module load gemini
# gemini query --header --show-samples --format sampledetail -q "select * from variants where aaf < 0.05 and filter is NULL AND (aaf_esp_all < 0.01 AND aaf_1kg_all < 0.01 AND af_exac_all < 0.01  AND (is_coding=1 OR is_splicing=1)) OR impact_severity='HIGH' OR clinvar_sig LIKE '%patho%' " vcfs.GATK.PED_columbia_master_pedigree.gemini.db | grep -v "FAIL_McGaughey" | bgzip -f >  columbia.gemini.tsv.gz

library(tidyverse)
library(data.table)
colombia_raw <- fread('gzcat /Volumes/data/projects/nei/hufnagel/colombia_cohort/columbia.gemini.tsv.gz')
colombia_raw %>% head()

