################################################
# give this 30 cores and 256GB of memory on biowulf
################################################


# args = commandArgs(trailingOnly=TRUE)
# uk10k_data <- args[1]
# clinvar_file <- args[2]
# gnomad_file <- args[3]

# biowulf paths
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
# Highly correlated predictors removed by inspecting correlation heatmap in remove_correlated_predictors.Rmd
# 
# or rather unique(ish) predictors kept
# reflected in the `select` commands below
###############################

###############################
# UK10K processing
###############################
load(uk10k_data)
## Set up data for modeling
all_processed <- uk10k_gemini_rare_variants %>% 
  filter(!(af_exac_all > 0.0001 & Genotype=='Het' & Status=='NotPathogenic')) %>%  # remove more common het mutations
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
  mutate_at(vars(matches('af_|dann|revel|mpc|gerp|polyphen_score|sift_score|fitcons_float|gerp_elements|^adj|_z$|^pli$|^pnull$|precessive|^phylop|linsight|metalr_score|lrt_converted_rankscore|metasvm|eigen|mutationtaster|provean')), funs(as.numeric(.))) %>%  # af is allele frequency
  select(variant_id, Status, Complicated_Status, is_exonic, is_coding, is_lof, is_splicing, impact_severity, polyphen_score, sift_score, dann, gerp_elements, DiseaseClass, mpc, revel, max_aaf_all, gno_ac_afr, gno_ac_eas, gno_ac_all, gno_ac_popmax, ac_exac_sas, ac_exac_fin, aaf_1kg_all_float, aaf_esp_all, ac_exac_all, ac_exac_amr, ac_exac_oth, gno_af_all, gno_an_popmax, an_exac_all, af_exac_all, fitcons_float, lof_z:precessive, phylop_100way, grantham, cadd_phred, fathmm_mkl_coding_score, linsight, metalr_score, lrt_converted_rankscore, metasvm_rankscore, eigen_phred, mutationtaster_score, provean_converted_rankscore, genesplicer, spliceregion) %>% 
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
#clinvar <- fread('gzcat ~/git/eye_var_Pathogenicity/processed_data/clinvar.gemini.tsv.gz')
clinvar <- fread(paste0('gzcat ', clinvar_file))

## Prep data for modeling
clinvar_processed <- clinvar %>% 
  #filter(status!='PATHOGENIC_OTHER') %>% # drop non eye pathogenic variants for the model learning 
  separate(gene_eyediseaseclass, c('RDGene','DiseaseClass'), sep='_') %>%  #split off RD disease type
  select(-RDGene) %>% 
  mutate(impact_severity = case_when(impact_severity == 'HIGH' ~ 3, # convert to integer 
                                     impact_severity == 'MED' ~ 2, 
                                     TRUE ~ 1),
         Status = case_when(status=='PATHOGENIC_EYE' ~ 'Pathogenic',
                            status=='PATHOGENIC_OTHER' ~ 'Pathogenic_NOTEYE',
                            TRUE ~ 'NotPathogenic'),
         genesplicer = case_when(genesplicer == "" ~ 'No',
                                 grepl('^gain', genesplicer) ~ 'Gain',
                                 grepl('^loss', genesplicer) ~ 'Loss',
                                 grepl('^diff', genesplicer) ~ 'Diff',
                                 TRUE ~ 'Else')) %>% 
  mutate_at(vars(matches('ac_|an_|^n_')), funs(as.integer(.))) %>% # convert columns with ac_|whatever to integer (ac is allele count)
  mutate_at(vars(matches('af_|dann|revel|mpc|gerp|polyphen_score|sift_score|fitcons_float|gerp_elements|^adj|_z$|^pli$|^pnull$|precessive|^phylop_100|linsight|metalr_score|lrt_converted_rankscore|metasvm|eigen|mutationtaster|provean')), funs(as.numeric(.))) %>%  # af is allele frequency
  select(variant_id, Status, is_exonic, is_coding, is_lof, is_splicing, impact_severity, polyphen_score, sift_score, dann, gerp_elements, DiseaseClass, mpc, revel, max_aaf_all, gno_ac_afr, gno_ac_eas, gno_ac_all, gno_ac_popmax, ac_exac_sas, ac_exac_fin, aaf_1kg_all_float, aaf_esp_all, ac_exac_all, ac_exac_amr, ac_exac_oth, gno_af_all, gno_an_popmax, an_exac_all, af_exac_all, fitcons_float, lof_z:precessive, phylop_100way, grantham, cadd_phred, fathmm_mkl_coding_score, linsight, metalr_score, lrt_converted_rankscore, metasvm_rankscore, eigen_phred, mutationtaster_score, provean_converted_rankscore, genesplicer, spliceregion) %>% 
  filter(max_aaf_all < 0.01) # remove any common variants

# fill missing with -1
clinvar_processed[is.na(clinvar_processed)] <- -1

ML_set__clinvar <- clinvar_processed %>% 
  filter(Status != 'Pathogenic_NOTEYE') %>% 
  mutate(Status = factor(Status, levels=c('Pathogenic','NotPathogenic'))) 
ML_set__clinvar$Source <- 'ClinVar'

ML_set__clinvar__otherPath <- clinvar_processed %>% 
  filter(Status == 'Pathogenic_NOTEYE') %>% 
  mutate(Status = gsub('Pathogenic_NOTEYE','Pathogenic',Status)) %>% 
  mutate(Status = factor(Status, levels=c('Pathogenic','NotPathogenic'))) 
ML_set__clinvar__otherPath$Source <- 'ClinVar'

###############################################
# gnomAD benign? processing
##############################################
gnomad <- fread(paste0('gzcat ', gnomad_file))
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
  mutate_at(vars(matches('af_|dann|revel|mpc|gerp|polyphen_score|sift_score|fitcons_float|gerp_elements|^adj|_z$|^pli$|^pnull$|precessive|^phylop_100|linsight|metalr_score|lrt_converted_rankscore|metasvm|eigen|mutationtaster|provean')), funs(as.numeric(.))) %>%  # af is allele frequency
  select(variant_id, Status, is_exonic, is_coding, is_lof, is_splicing, impact_severity, polyphen_score, sift_score, dann, gerp_elements, DiseaseClass, mpc, revel, max_aaf_all, gno_ac_afr, gno_ac_eas, gno_ac_all, gno_ac_popmax, ac_exac_sas, ac_exac_fin, aaf_1kg_all_float, aaf_esp_all, ac_exac_all, ac_exac_amr, ac_exac_oth, gno_af_all, gno_an_popmax, an_exac_all, af_exac_all, fitcons_float, lof_z:precessive, phylop_100way, grantham, cadd_phred, fathmm_mkl_coding_score, linsight, metalr_score, lrt_converted_rankscore, metasvm_rankscore, eigen_phred, mutationtaster_score, provean_converted_rankscore, genesplicer, spliceregion) %>% 
  filter(max_aaf_all < 0.01) # remove any common variants

# fill missing with -1
gnomad_processed[is.na(gnomad_processed)] <- -1

# add 250X the number of clinvar pathogenic variants 
set.seed(13457)
gnomad_processed_sub <- gnomad_processed %>% sample_n((ML_set__clinvar %>% filter(Status=='Pathogenic') %>% nrow()) * 250)
gnomad_processed_other <- gnomad_processed %>% filter(!variant_id %in% gnomad_processed_sub$variant_id) # not used for model building, for potential validation purposes
gnomad_processed_sub$Source <- 'gnomAD'
gnomad_processed_other$Source <- 'gnomAD'
################################################
# Combine UK10K and ClinVar and gnomAD data
################################################

ML_set__all <- bind_rows(ML_set__clinvar %>% select_(.dots = colnames(ML_set__UK10K %>% select(-Complicated_Status))), 
                         ML_set__UK10K %>% select(-Complicated_Status),
                         gnomad_processed_sub %>% select_(.dots = colnames(ML_set__UK10K %>% select(-Complicated_Status))))

ML_set__other <- bind_rows(gnomad_processed_other, ML_set__clinvar__otherPath)

################################
# one hot encode
##################################

temp <- ML_set__all %>% dplyr::select(-Status, -Source, -variant_id)
temp <- dummy.data.frame(temp, sep='_')
ML_set_dummy <- temp %>% mutate(variant_id = ML_set__all$variant_id, Status = ML_set__all$Status, Source = ML_set__all$Source)
# secondary set with non eye path variants
temp <- ML_set__other  %>% dplyr::select(-Status, -Source, -variant_id)
temp <- dummy.data.frame(temp, sep='_')
ML_set_dummy__secondary <- temp %>% mutate(variant_id = ML_set__other $variant_id, Status = ML_set__other $Status, Source = ML_set__other $Source)


# center scale 
# ML_set_dummy_CS <- preProcess(ML_set_dummy, method = c('center','scale')) %>% predict(., ML_set_dummy)

##################################
# train, validate, and test sets
# half to train
# 1/4 to validate
# 1/4 to test (only use when I think I'm done building models)
##################################

set.seed(115470)
train_set <- ML_set_dummy %>% 
  group_by(Status, Source) %>% 
  # filter(!Complicated_Status=='Comp_Het') %>% # remove comp hets for now
  sample_frac(0.5) %>% ungroup()

set.seed(115470)
validate_set <- ML_set_dummy %>% 
  filter(!variant_id %in% train_set$variant_id) %>% 
  group_by(Status, Source) %>% 
  sample_frac(0.5) %>% ungroup()

test_set <- ML_set_dummy %>% 
  filter(!variant_id %in% c(train_set$variant_id, validate_set$variant_id))


########################################
# Set up fitcontrols for caret modeling
#########################################
# CV on rf seems to overfit
fitControl_min <- trainControl(
  classProbs=T,
  savePredictions = T,
  allowParallel = T,
  summaryFunction = prSummary,
  returnData = T)

# fitControl_min_CS <- trainControl(
#   classProbs=T,
#   savePredictions = T,
#   allowParallel = T,
#   summaryFunction = prSummary,
#   returnData = T,
#   preProcOptions = c('center','scale'))

# fitControl <- trainControl(## 5-fold CV
#   method = "cv",
#   number = 5,
#   classProbs=T,
#   savePredictions = T,
#   allowParallel = T,
#   summaryFunction = prSummary,
#   returnData = T)

###########################################
# SAVE DATA
##########################################
model_run <- list()
model_run$most_imp_predictors <- most_imp_predictors
model_run$most_imp_predictors_no_disease_class <- most_imp_predictors_no_disease_class
model_run$train_set <- train_set 
model_run$validate_set <- validate_set
model_run$test_set <- test_set
model_run$ML_set_dummy <- ML_set_dummy__secondary
model_run$sessionInfo <- sessionInfo()
save(model_run, file='model_run__2018_03_24.Rdata')

###########################################
# multi processing
##########################################
cluster <- makeCluster(30) 
registerDoParallel(cluster)

##############################################
# BUILD MODELS!!!!!!!!!!!
#############################################
rfFit_all <- caret::train(Status ~ ., data=train_set %>% select(-variant_id, -Source), 
                      method = "rf", metric='F',
                      trControl = fitControl_min)
# use the first rf model to pick the useful predictors and limit the models to these
most_imp_predictors <- varImp(rfFit_all)$importance  %>% rownames_to_column('Predictors') %>% arrange(-Overall) %>% filter(Overall > 2) %>% pull(Predictors)
# variant with no disease class predictos
most_imp_predictors_no_disease_class <- most_imp_predictors[!grepl('DiseaseClass', most_imp_predictors)]


rfFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                      method = "rf", metric='F',
                      trControl = fitControl_min)

rfFit_noDC <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors_no_disease_class)), 
                           method = "rf", metric='F',
                           trControl = fitControl_min)

bglmFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                        method = "bayesglm", metric='F',
                        trControl = fitControl_min)

bglmFit_noDC <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors_no_disease_class)), 
                             method = "bayesglm", metric='F',
                             trControl = fitControl_min)

glmboostFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                            method = "glmboost", metric='F',
                            trControl = fitControl_min,
                            preProcess = c('center','scale'))

LogitBoostFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                              method = "LogitBoost", metric='F',
                              trControl = fitControl_min,
                              preProcess = c('center','scale'))

avNNetFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                          method = "avNNet", metric='F',
                          trControl = fitControl_min,
                          preProcess = c('center','scale'))

avNNetFit_noDC <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors_no_disease_class)), 
                               method = "avNNet", metric='F',
                               trControl = fitControl_min,
                               preProcess = c('center','scale'))

# tossed, terrible performance
# xgbTreeFit <- caret::train(Status ~ ., data=train_set %>% select(-variant_id, -Source), 
#                      
#                       method = "xgbTree",  metric='AUC',
#                       trControl = fitControl)

# stepLDAFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
#                            method = "stepLDA",  metric='AUC',
#                            trControl = fitControl_min_CS)

# tossed, terrible performance
# naive_bayesFit <- caret::train(Status ~ ., data=train_set %>% select(-variant_id, -Source), 
#                      
#                       method = "naive_bayes",  metric='AUC',
#                       trControl = fitControl)

# dnnFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
#                        method = "dnn",  metric='AUC',
#                        trControl = fitControl_min,
#                        preProcess = c('center','scale'))
# 
# dnnFit_noDC <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors_no_disease_class)), 
#                        method = "dnn",  metric='AUC',
#                        trControl = fitControl_min,
#                        preProcess = c('center','scale'))


# monmlpFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
#                           method = "monmlp",  metric='AUC',
#                           trControl = fitControl_min)


# monmlpFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
#                           method = "monmlp",  metric='F',
#                           trControl = fitControl_min)

# mlpKerasDropoutFit <- caret::train(x = train_set %>% select(-variant_id, -Source, -Status),
#                                  y= train_set$Status,
#                      
#                       method = "mlpKerasDropout",  metric='F',
#                       trControl = fitControl)

# svmLinearWeightsFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
#                                     method = "svmLinearWeights",  metric='Precision',
#                                     trControl = fitControl_min,
#                                     preProcess = c('center','scale'))

caddFit <- caret::train(Status ~ ., data=train_set %>% select(cadd_phred, Status), 
                        method = "glm",  metric='AUC',
                        trControl = fitControl_min)

cadd_plus_DiseaseClassFit <- caret::train(Status ~ ., data=train_set %>% select(cadd_phred, Status, `DiseaseClass_-1`), 
                                          method = "glm",  metric='AUC',
                                          trControl = fitControl_min)

revelFit <- caret::train(Status ~ ., data=train_set %>% select(revel, Status), 
                         method = "glm", metric='AUC',
                         trControl = fitControl_min)

revel_plus_DiseaseClassFit <- caret::train(Status ~ ., data=train_set %>% select(revel, Status, `DiseaseClass_-1`), 
                                           method = "glm", metric='AUC',
                                           trControl = fitControl_min)

dannFit <- caret::train(Status ~ ., data=train_set %>% select(dann, Status), 
                        method = "glm",
                        trControl = fitControl_min)

dann_plus_DiseaseClassFit <- caret::train(Status ~ ., data=train_set %>% select(dann, Status, `DiseaseClass_-1`), 
                                          method = "glm",
                                          trControl = fitControl_min)

##############################
# SAVE MODELS
###############################

for (i in ls()[grepl('Fit',ls())]) {model_run[[i]] <- get(i)}
save(model_run, file='model_run__2018_03_26.Rdata')
