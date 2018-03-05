#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)


# Load data from ~/git/EGA_EGAD00001002656_NGS_reanalyze/scripts/stats.Rmd
load('uk10k_gemini_rare_variants.Rdata')
## Set up data for modeling
library(tidyverse)

all_processed <- uk10k_gemini_rare_variants %>% 
  separate(gene_eyediseaseclass, c('RDGene','DiseaseClass'), sep='_') %>%  #split off RD disease type
  select(-RDGene) %>% 
  mutate(impact_severity = case_when(impact_severity == 'HIGH' ~ 3, # convert to integer 
                                     impact_severity == 'MED' ~ 2, 
                                     TRUE ~ 1),
         Status = case_when(Status=='Pathogenic' ~ 'Pathogenic',
                            TRUE ~ 'NotPathogenic')) %>% 
  mutate(Status = factor(Status, levels=c('Pathogenic','NotPathogenic'))) %>% 
  mutate_at(vars(matches('ac_|an_|^n_')), funs(as.integer(.))) %>% # convert columns with ac_|whatever to integer (ac is allele count)
  mutate_at(vars(matches('af_|dann|revel|mpc|gerp|polyphen_score|sift_score|fitcons_float|gerp_elements|^adj|_z$|^pli$|^pnull$|precessive|^phylop')), funs(as.numeric(.))) %>%  # af is allele frequency
  select(variant_id, Status, Complicated_Status, is_exonic, is_coding, is_lof, is_splicing, impact_severity, polyphen_score, sift_score, dann, gerp_elements, DiseaseClass, mpc, revel, aaf_1kg_afr_float:an_exac_sas, fitcons_float, gno_ac_afr:gno_an_popmax, lof_z:precessive, phylop_100way, grantham, maxentscan, cadd_phred, fathmm_mkl_coding_score, genesplicer, spliceregion) %>% 
  filter(max_aaf_all < 0.01) %>% 
  unique() # remove any common variants

all_processed[is.na(all_processed)] <- -1

all_PATH <- all_processed %>% 
  filter(Status=='Pathogenic') %>% 
  unique()
all_NOT_PATH <- all_processed %>% 
  filter(Status=='NotPathogenic') %>% 
  unique()

all_set__uk10k <- all_processed
all_set__uk10k$Source <- 'UK10K'

all_PATH <- all_set__uk10k %>% filter(Status=='Pathogenic')
# cut down pathogenic 
set.seed(115470)
all_NOT_PATH__CUT <- all_set__uk10k %>% filter(Status=='NotPathogenic') #%>% sample_n(20000)

ML_set__UK10K <- rbind(all_PATH, all_NOT_PATH__CUT)

# Load clinvar 
library(data.table)
clinvar <- fread('zcat ~/git/eye_var_Pathogenicity/processed_data/clinvar.gemini.tsv.gz')


## Prep data for modeling
clinvar_processed <- clinvar %>% 
  separate(gene_eyediseaseclass, c('RDGene','DiseaseClass'), sep='_') %>%  #split off RD disease type
  select(-RDGene) %>% 
  mutate(impact_severity = case_when(impact_severity == 'HIGH' ~ 3, # convert to integer 
                                     impact_severity == 'MED' ~ 2, 
                                     TRUE ~ 1),
         Status = case_when(status=='PATHOGENIC_RD' ~ 'Pathogenic',
                            TRUE ~ 'NotPathogenic'),
         genesplicer = case_when(genesplicer == "" ~ 'No',
                                 grepl('^gain', genesplicer) ~ 'Gain',
                                 grepl('^loss', genesplicer) ~ 'Loss',
                                 grepl('^diff', genesplicer) ~ 'Diff',
                                 TRUE ~ 'Else')) %>% 
  mutate(Status = factor(Status, levels=c('Pathogenic','NotPathogenic'))) %>% 
  mutate_at(vars(matches('ac_|an_|^n_')), funs(as.integer(.))) %>% # convert columns with ac_|whatever to integer (ac is allele count)
  mutate_at(vars(matches('af_|dann|revel|mpc|gerp|polyphen_score|sift_score|fitcons_float|gerp_elements|^adj|_z$|^pli$|^pnull$|precessive|^phylop_100')), funs(as.numeric(.))) %>%  # af is allele frequency
  select(variant_id, Status, is_exonic, is_coding, is_lof, is_splicing, impact_severity, polyphen_score, sift_score, dann, gerp_elements, DiseaseClass, mpc, revel, aaf_1kg_afr_float:an_exac_sas, fitcons_float, gno_ac_afr:gno_an_popmax, lof_z:precessive, phylop_100way, grantham, maxentscan, cadd_phred, fathmm_mkl_coding_score, genesplicer, spliceregion) %>% 
  filter(max_aaf_all < 0.01) # remove any common variants

clinvar_processed[is.na(clinvar_processed)] <- -1

all_PATH <- clinvar_processed %>% 
  filter(Status == 'Pathogenic') %>% 
  unique()
all_NOT_PATH <- clinvar_processed %>% 
  filter(Status != 'Pathogenic') %>% 
  unique()

all_set__clinvar <- clinvar_processed
all_set__clinvar$Source <- 'ClinVar'

all_PATH <- all_set__clinvar %>% filter(Status=='Pathogenic')
# cut down pathogenic to 5x of path
set.seed(115470)
all_NOT_PATH__CUT <- all_set__clinvar %>% filter(Status=='NotPathogenic')

ML_set__clinvar <- rbind(all_PATH, all_NOT_PATH__CUT)

# Combine UK10K and ClinVar data
ML_set__all <- bind_rows(ML_set__clinvar %>% select_(.dots = colnames(ML_set__UK10K %>% select(-Complicated_Status))), ML_set__UK10K %>% select(-Complicated_Status))


set.seed(115470)
train_set <- ML_set__all %>% 
  group_by(Status, Source) %>% 
  # filter(!Complicated_Status=='Comp_Het') %>% # remove comp hets for now
  sample_frac(0.33) %>% ungroup()

set.seed(115470)
validate_set <- ML_set__all %>% 
  filter(!variant_id %in% train_set$variant_id) %>% 
  group_by(Status, Source) %>% 
  sample_frac(0.5) %>% ungroup()

library(caret)
library(mlbench)
library(parallel)
library(doParallel)
library(MLmetrics)
#library(plotROC)

# CV on rf seems to overfit
fitControl_RF <- trainControl(
  classProbs=T,
  savePredictions = T,
  allowParallel = T,
  summaryFunction = prSummary,
  returnData = T)


cluster <- makeCluster(24) 
registerDoParallel(cluster)

rfFit <- caret::train(Status ~ ., data=train_set %>% select(-variant_id, -Source), 
                      preProcess = c("scale", "center"),
                      method = "rf", metric='F',
                      trControl = fitControl_RF)
save(rfFit, file='rfFit_2018_03_03.Rdata')

rfFit_noDC <- caret::train(Status ~ ., data=train_set %>% select(-variant_id, -Source, -DiseaseClass), 
                      preProcess = c("scale", "center"),
                      method = "rf", metric='F',
                      trControl = fitControl_RF)
save(rfFit_noDC, file='rfFit_noDC_2018_03_03.Rdata')