# run in /data/mcgaugheyd/projects/nei/hufnagel/colombia_cohort
# module load gemini
# gemini query --header --show-samples --format sampledetail -q "select * from variants where aaf < 0.05 and filter is NULL AND (aaf_esp_all < 0.01 AND aaf_1kg_all < 0.01 AND af_exac_all < 0.01  AND (is_coding=1 OR is_splicing=1)) OR impact_severity='HIGH' OR clinvar_sig LIKE '%patho%' "  --gt-filter "(gt_types).(phenotype==2).(==HET).(any)" vcfs.GATK.PED_columbia_master_pedigree.gemini.db | grep -v "FAIL_McGaughey" | bgzip -f >  columbia.gemini.het.tsv.gz &
# gemini query --header --show-samples --format sampledetail -q "select * from variants where aaf < 0.05 and filter is NULL AND (aaf_esp_all < 0.01 AND aaf_1kg_all < 0.01 AND af_exac_all < 0.01  AND (is_coding=1 OR is_splicing=1)) OR impact_severity='HIGH' OR clinvar_sig LIKE '%patho%' "  --gt-filter "(gt_types).(phenotype==2).(==HOM).(any)" vcfs.GATK.PED_columbia_master_pedigree.gemini.db | grep -v "FAIL_McGaughey" | bgzip -f >  columbia.gemini.hom.tsv.gz &

library(tidyverse)
library(data.table)
colombia_raw <- rbind(fread('gzcat /Volumes/data/projects/nei/hufnagel/colombia_cohort/columbia.gemini.het.tsv.gz') %>% mutate(Allele = 'Het'),
                      fread('gzcat /Volumes/data/projects/nei/hufnagel/colombia_cohort/columbia.gemini.hom.tsv.gz') %>% mutate(Allele = 'Hom'))
colombia_meta <- bind_rows(readxl::read_xlsx('/Volumes/data/projects/nei/hufnagel/colombia_cohort/summary_variants_solved.xlsx', sheet = 'Annovar output'),
                           readxl::read_xlsx('/Volumes/data/projects/nei/hufnagel/colombia_cohort/summary_variants_solved.xlsx', sheet = 'IVA output') %>% mutate(End = `End Position`, Start = Position))


# ID path
path <- bind_rows(colombia_raw %>% right_join(., colombia_meta, by = c("end" = "End")) %>% select(chrom, start, end, ref, alt, Allele, name, `Study ID`, `Gene.refGene`),           colombia_raw %>% right_join(., colombia_meta, by = c("end" = "Start")) %>% 
                    select(variant_id, chrom, start, end, ref, alt, Allele, name, `Study ID`, `Gene.refGene`)) %>% 
  unique() %>% 
  filter(!is.na(start)) %>% 
  filter(start != 11200234)

#  samples with path vars
path$`Study ID` %>% unique()

# remove variants from unsolved
colombia_out <- colombia_raw %>% 
  filter(name %in% gsub('-', '_', path$`Study ID` %>% unique())) %>% 
  select(-name, -family_id, -sex, -phenotype, -Allele) %>% 
  unique() %>% 
  mutate(Status = case_when(variant_id %in% path$variant_id ~ 'Pathogenic',
         TRUE ~ 'NotPathogenic')) 

save(colombia_out, file = '/Volumes/data/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/colombia.Rdata')
