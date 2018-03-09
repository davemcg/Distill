#!/usr/bin/env Rscript

# give this 30 cores and 200GB of memory on biowulf

# args = commandArgs(trailingOnly=TRUE)
# uk10k_data <- args[1]
# clinvar_file <- args[2]
# gnomad_file <- args[3]

uk10k_data <- '/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/uk10k_gemini_rare_variants.Rdata'
clinvar_file <- '/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/clinvar.gemini.tsv.gz'
gnomad_file <- '/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/gnomad.gemini.tsv.gz'

library(tidyverse)
library(data.table)
library(dummies)
library(caret)
library(mlbench)
library(parallel)
library(doParallel)
library(MLmetrics)


###############################
# UK10K processing
###############################
# Load data from ~/git/EGA_EGAD00001002656_NGS_reanalyze/scripts/stats.Rmd
#load('uk10k_gemini_rare_variants.Rdata')
load(uk10k_data)

## Set up data for modeling
all_processed <- uk10k_gemini_rare_variants %>% 
  separate(gene_eyediseaseclass, c('RDGene','DiseaseClass'), sep='_') %>%  #split off RD disease type
  select(-RDGene) %>% 
  mutate(impact_severity = case_when(impact_severity == 'HIGH' ~ 3, # convert to integer 
                                     impact_severity == 'MED' ~ 2, 
                                     TRUE ~ 1),
         Status = case_when(Status=='Pathogenic' ~ 'Pathogenic',
                            TRUE ~ 'NotPathogenic'),
         genesplicer = case_when(genesplicer == "" ~ 'No',
                                 grepl('^gain', genesplicer) ~ 'Gain',
                                 grepl('^loss', genesplicer) ~ 'Loss',
                                 grepl('^diff', genesplicer) ~ 'Diff',
                                 TRUE ~ 'Else')) %>% 
  mutate(Status = factor(Status, levels=c('Pathogenic','NotPathogenic'))) %>% 
  mutate_at(vars(matches('ac_|an_|^n_')), funs(as.integer(.))) %>% # convert columns with ac_|whatever to integer (ac is allele count)
  mutate_at(vars(matches('af_|dann|revel|mpc|gerp|polyphen_score|sift_score|fitcons_float|gerp_elements|^adj|_z$|^pli$|^pnull$|precessive|^phylop')), funs(as.numeric(.))) %>%  # af is allele frequency
  select(variant_id, Status, Complicated_Status, is_exonic, is_coding, is_lof, is_splicing, impact_severity, polyphen_score, sift_score, dann, gerp_elements, DiseaseClass, mpc, revel, aaf_1kg_afr_float:an_exac_sas, fitcons_float, gno_ac_afr:gno_an_popmax, lof_z:precessive, phylop_100way, grantham, cadd_phred, fathmm_mkl_coding_score, genesplicer, spliceregion) %>% 
  filter(max_aaf_all < 0.01) %>% 
  unique() # remove any common variants

# fill missing with -1
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
# chance to cut down non pathogenic
# i'm just keeping all right now
set.seed(115470)
all_NOT_PATH__CUT <- all_set__uk10k %>% filter(Status=='NotPathogenic') #%>% sample_n(20000)

ML_set__UK10K <- rbind(all_PATH, all_NOT_PATH__CUT)

###########################################
# ClinVar  Processing
###########################################
#clinvar <- fread('zcat ~/git/eye_var_Pathogenicity/processed_data/clinvar.gemini.tsv.gz')
clinvar <- fread(paste0('zcat ', clinvar_file))

## Prep data for modeling
clinvar_processed <- clinvar %>% 
  separate(gene_eyediseaseclass, c('RDGene','DiseaseClass'), sep='_') %>%  #split off RD disease type
  select(-RDGene) %>% 
  mutate(impact_severity = case_when(impact_severity == 'HIGH' ~ 3, # convert to integer 
                                     impact_severity == 'MED' ~ 2, 
                                     TRUE ~ 1),
         Status = case_when(status=='PATHOGENIC_EYE' ~ 'Pathogenic',
                            TRUE ~ 'NotPathogenic'),
         genesplicer = case_when(genesplicer == "" ~ 'No',
                                 grepl('^gain', genesplicer) ~ 'Gain',
                                 grepl('^loss', genesplicer) ~ 'Loss',
                                 grepl('^diff', genesplicer) ~ 'Diff',
                                 TRUE ~ 'Else')) %>% 
  mutate(Status = factor(Status, levels=c('Pathogenic','NotPathogenic'))) %>% 
  mutate_at(vars(matches('ac_|an_|^n_')), funs(as.integer(.))) %>% # convert columns with ac_|whatever to integer (ac is allele count)
  mutate_at(vars(matches('af_|dann|revel|mpc|gerp|polyphen_score|sift_score|fitcons_float|gerp_elements|^adj|_z$|^pli$|^pnull$|precessive|^phylop_100')), funs(as.numeric(.))) %>%  # af is allele frequency
  select(variant_id, Status, is_exonic, is_coding, is_lof, is_splicing, impact_severity, polyphen_score, sift_score, dann, gerp_elements, DiseaseClass, mpc, revel, aaf_1kg_afr_float:an_exac_sas, fitcons_float, gno_ac_afr:gno_an_popmax, lof_z:precessive, phylop_100way, grantham, cadd_phred, fathmm_mkl_coding_score, genesplicer, spliceregion) %>% 
  filter(max_aaf_all < 0.01) # remove any common variants

# fill missing with -1
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


###############################################
# gnomAD benign? processing
##############################################
gnomad <- fread(paste0('zcat ', gnomad_file))
## Prep data for modeling
gnomad_processed <- gnomad %>% 
  separate(gene_eyediseaseclass, c('RDGene','DiseaseClass'), sep='_') %>%  #split off RD disease type
  select(-RDGene) %>% 
  mutate(impact_severity = case_when(impact_severity == 'HIGH' ~ 3, # convert to integer 
                                     impact_severity == 'MED' ~ 2, 
                                     TRUE ~ 1),
         Status = 'NotPathogenic',
         genesplicer = case_when(genesplicer == "" ~ 'No',
                                 grepl('^gain', genesplicer) ~ 'Gain',
                                 grepl('^loss', genesplicer) ~ 'Loss',
                                 grepl('^diff', genesplicer) ~ 'Diff',
                                 TRUE ~ 'Else')) %>% 
  mutate(Status = factor(Status, levels=c('Pathogenic','NotPathogenic'))) %>% 
  filter(hgmd_overlap=='None' & clinvar_pathogenic == 'None') %>% # remove possible pathogenic by checking against hgmd or clinvar presence
  mutate_at(vars(matches('ac_|an_|^n_')), funs(as.integer(.))) %>% # convert columns with ac_|whatever to integer (ac is allele count)
  mutate_at(vars(matches('af_|dann|revel|mpc|gerp|polyphen_score|sift_score|fitcons_float|gerp_elements|^adj|_z$|^pli$|^pnull$|precessive|^phylop_100')), funs(as.numeric(.))) %>%  # af is allele frequency
  select(variant_id, Status, is_exonic, is_coding, is_lof, is_splicing, impact_severity, polyphen_score, sift_score, dann, gerp_elements, DiseaseClass, mpc, revel, aaf_1kg_afr_float:an_exac_sas, fitcons_float, gno_ac_afr:gno_an_popmax, lof_z:precessive, phylop_100way, grantham, cadd_phred, fathmm_mkl_coding_score, genesplicer, spliceregion) %>% 
  filter(max_aaf_all < 0.01) # remove any common variants

# fill missing with -1
gnomad_processed[is.na(gnomad_processed)] <- -1

# add 250X the number of clinvar pathogenic variants 
set.seed(13457)
gnomad_processed_sub <- gnomad_processed %>% sample_n((ML_set__clinvar %>% filter(Status=='Pathogenic') %>% nrow()) * 250)
gnomad_processed_sub$Source <- 'gnomAD'

################################################
# Combine UK10K and ClinVar and gnomAD data
################################################

ML_set__all <- bind_rows(ML_set__clinvar %>% select_(.dots = colnames(ML_set__UK10K %>% select(-Complicated_Status))), 
                         ML_set__UK10K %>% select(-Complicated_Status),
                         gnomad_processed_sub %>% select_(.dots = colnames(ML_set__UK10K %>% select(-Complicated_Status))))

################################
# one hot encode and center scale
##################################

temp <- ML_set__all %>% dplyr::select(-Status, -Source, -variant_id)
temp <- dummy.data.frame(temp, sep='_')
ML_set_dummy <- temp %>% mutate(variant_id = ML_set__all$variant_id, Status = ML_set__all$Status, Source = ML_set__all$Source)
# center scale
ML_set_dummy_CS <- preProcess(ML_set_dummy, method = c('center','scale')) %>% predict(., ML_set_dummy)

set.seed(115470)
train_set <- ML_set_dummy_CS %>% 
  group_by(Status, Source) %>% 
  # filter(!Complicated_Status=='Comp_Het') %>% # remove comp hets for now
  sample_frac(0.33) %>% ungroup()

set.seed(115470)
validate_set <- ML_set_dummy_CS %>% 
  filter(!variant_id %in% train_set$variant_id) %>% 
  group_by(Status, Source) %>% 
  sample_frac(0.5) %>% ungroup()

test_set <- all_processed %>% 
  filter(!variant_id %in% c(train_set$variant_id, validate_set$variant_id))


########################################
# Set up fitcontrols for caret modeling
#########################################
# CV on rf seems to overfit
fitControl_RF <- trainControl(
  classProbs=T,
  savePredictions = T,
  allowParallel = T,
  summaryFunction = prSummary,
  returnData = T)

fitControl <- trainControl(## 5-fold CV
  method = "cv",
  number = 5,
  classProbs=T,
  savePredictions = T,
  allowParallel = T,
  summaryFunction = prSummary,
  returnData = T)

fitControl_min <- trainControl(
  classProbs=T,
  savePredictions = T,
  allowParallel = T,
  summaryFunction = prSummary,
  returnData = T)


###########################################
# multi processing
##########################################
cluster <- makeCluster(30) 
registerDoParallel(cluster)

##############################################
# BUILD MODELS!!!!!!!!!!!
#############################################
rfFit <- caret::train(Status ~ ., data=train_set %>% select(-variant_id, -Source), 
                      method = "rf", metric='F',
                      trControl = fitControl_RF)
# use the first rf model to pick the useful predictors and limit the models to these
most_imp_predictors <- varImp(rfFit)$importance  %>% rownames_to_column('Predictors') %>% arrange(-Overall) %>% filter(Overall > 2) %>% pull(Predictors)
# variant with no disease class predictos
most_imp_predictors_no_disease_class <- most_imp_predictors[!grepl('DiseaseClass', most_imp_predictors)]


rfFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                      method = "rf", metric='F',
                      trControl = fitControl_RF)

rfFit_noDC <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors_no_disease_class)), 
                           method = "rf", metric='F',
                           trControl = fitControl_RF)

bglmFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                        method = "bayesglm", metric='F',
                        trControl = fitControl)

bglmFit_noDC <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors_no_disease_class)), 
                             method = "bayesglm", metric='F',
                             trControl = fitControl)

# check again, had errors when running
glmboostFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                            method = "glmboost", metric='F',
                            trControl = fitControl)

LogitBoostFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                              method = "LogitBoost", metric='F',
                              trControl = fitControl)

avNNetFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                          method = "avNNet", metric='F',
                          trControl = fitControl)

avNNetFit_noDC <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors_no_disease_class)), 
                               method = "avNNet", metric='F',
                               trControl = fitControl)

# tossed, terrible performance
# xgbTreeFit <- caret::train(Status ~ ., data=train_set %>% select(-variant_id, -Source), 
#                      
#                       method = "xgbTree",  metric='AUC',
#                       trControl = fitControl)

stepLDAFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                           method = "stepLDA",  metric='AUC',
                           trControl = fitControl)

# tossed, terrible performance
# naive_bayesFit <- caret::train(Status ~ ., data=train_set %>% select(-variant_id, -Source), 
#                      
#                       method = "naive_bayes",  metric='AUC',
#                       trControl = fitControl)

dnnFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                       method = "dnn",  metric='AUC',
                       trControl = fitControl)


monmlpFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                          method = "monmlp",  metric='AUC',
                          trControl = fitControl_min)


# monmlpFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
#                           method = "monmlp",  metric='F',
#                           trControl = fitControl_min)

# mlpKerasDropoutFit <- caret::train(x = train_set %>% select(-variant_id, -Source, -Status),
#                                  y= train_set$Status,
#                      
#                       method = "mlpKerasDropout",  metric='F',
#                       trControl = fitControl)


svmLinearWeightsFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                                    method = "svmLinearWeights",  metric='Precision',
                                    trControl = fitControl)

caddFit <- caret::train(Status ~ ., data=train_set %>% select(cadd_phred, Status), 
                        method = "glm",  metric='Precision',
                        trControl = fitControl)

cadd_plus_DiseaseClassFit <- caret::train(Status ~ ., data=train_set %>% select(cadd_phred, Status, `DiseaseClass_-1`), 
                                          method = "glm",  metric='Precision',
                                          trControl = fitControl)

revelFit <- caret::train(Status ~ ., data=train_set %>% select(revel, Status), 
                         method = "glm",
                         trControl = fitControl)

revel_plus_DiseaseClassFit <- caret::train(Status ~ ., data=train_set %>% select(revel, Status, `DiseaseClass_-1`), 
                                           method = "glm",
                                           trControl = fitControl)

dannFit <- caret::train(Status ~ ., data=train_set %>% select(dann, Status), 
                        method = "glm",
                        trControl = fitControl)

dann_plus_DiseaseClassFit <- caret::train(Status ~ ., data=train_set %>% select(dann, Status, `DiseaseClass_-1`), 
                                          method = "glm",
                                          trControl = fitControl)

##############################
# SAVE DATA AND MODELS
###############################

my_models <-list()
for (i in ls()[grepl('Fit',ls())]) {my_models[[i]] <- get(i)}
my_models$most_imp_predictors <- most_imp_predictors
my_models$most_imp_predictors_no_disease_class <- most_imp_predictors_no_disease_class
my_models$train_set <- train_set 
my_models$validate_set <- validate_set
my_models$test_set <- test_set
save(my_models, file='eye_var_path_models__2018_03_09.Rdata')
