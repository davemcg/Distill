######################
# biowulf2 note
# run with:
# sinteractive --gres=gpu:k80:1,lscratch:10 --mem=64g -c2
# module load cuDNN/7.0/CUDA-9.0 CUDA R/3.5.0 python/3.5
# to load keras/tensorflow 
# https://hpc.nih.gov/apps/caret.html
######################

# create unified df for model_assess.Rmd

library(tidyverse)
library(caret)
library(ModelMetrics)
library(randomForest)

# Load processed data and models

load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/model_data_2018_06_20.Rdata')

load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_6mtry.Rdata')
load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_9mtry.Rdata')
load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_12mtry.Rdata')
load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_15mtry.Rdata')
# # ogvfb cohort
# load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/ogvfb_exome_cohort.Rdata')
# # grimm
# load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/grimm_ML.Rdata')
# # unifun
# load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/unifun_ML.Rdata')
# # ddl nisc panel
# load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/ddl_nisc_panel_variants.Rdata')
# #VPaC_deeprnn
# load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/DeepRNN_model.Rdata')
# # samocha
# load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/samocha_ML.Rdata')
# # homsy
# load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/homsy_ML.Rdata')

# all raw
######
# only load if you need it, it's HUGE
######
load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/master/raw_data.Rdata')


numeric_predictors<-c('is_exonic','is_coding','is_lof','is_splicing','exon','aa_length','impact_severity','polyphen_score','sift_score','dann','eigen_pc_raw_rankscore','eigen_phred','eigen_raw','eigen_coding_or_noncoding','fathmm_converted_rankscore','fathmm_pred','fathmm_score','gerp','genocanyon_score','genocanyon_score_rankscore','hgmd_overlap','linsight','lrt_omega','lrt_converted_rankscore','lrt_score','m_cap_rankscore','m_cap_score','mpc','metalr_rankscore','metalr_score','metasvm_rankscore','metasvm_score','mutationassessor_score','mutationassessor_score_rankscore','mutationtaster_converted_rankscore','mutationtaster_score','provean_converted_rankscore','provean_score','revel','vest3_rankscore','vest3_score','aaf_1kg_afr','aaf_1kg_all','aaf_1kg_amr','aaf_1kg_eas','aaf_1kg_eur','aaf_1kg_sas','aaf_esp_aa','aaf_esp_all','aaf_esp_ea','ac_exac_afr','ac_exac_all','ac_exac_amr','ac_exac_eas','ac_exac_fin','ac_exac_nfe','ac_exac_oth','ac_exac_sas','adj_exp_lof','adj_exp_mis','adj_exp_syn','af_exac_afr','af_exac_all','af_exac_amr','af_exac_eas','af_exac_nfe','af_exac_oth','af_exac_sas','an_exac_afr','an_exac_all','an_exac_amr','an_exac_eas','an_exac_fin','an_exac_nfe','an_exac_oth','an_exac_sas','ccr_pct_v1','cpg_island','epilogos_bivflnk','epilogos_enh','epilogos_enhbiv','epilogos_enhg','epilogos_het','epilogos_quies','epilogos_reprpc','epilogos_reprpcwk','epilogos_tss','epilogos_tssaflnk','epilogos_tssbiv','epilogos_tx','epilogos_txflnk','epilogos_txwk','epilogos_znf','exac_num_het','exac_num_hom_alt','fathmm_mkl_coding_group','fathmm_mkl_coding_pred','fathmm_mkl_coding_rankscore','fathmm_mkl_coding_score','fitcons','geno2mp','gerp_elements','gno_ac_afr','gno_ac_all','gno_ac_amr','gno_ac_asj','gno_ac_eas','gno_ac_fin','gno_ac_nfe','gno_ac_oth','gno_ac_popmax','gno_af_afr','gno_af_all','gno_af_amr','gno_af_asj','gno_af_eas','gno_af_fin','gno_af_nfe','gno_af_oth','gno_af_popmax','gno_an_afr','gno_an_all','gno_an_amr','gno_an_asj','gno_an_eas','gno_an_fin','gno_an_nfe','gno_an_oth','gno_an_popmax','gno_id','gno_popmax','hapmap1','hapmap2','in_1kg','in_esp','in_exac','lof_z','max_aaf_all','mis_z','n_lof','n_mis','n_syn','pli','pnull','precessive','phylop_100way','segway_mean_score','segway_sum_score','stam_mean','syn_z','grantham','maxentscan','cadd_raw','cadd_phred')
############ 
###  DDL ###
############
panel <- readxl::read_excel('~/git/eye_var_Pathogenicity/data/NISC100_Variant_Interpretation_June01_2018.xlsx')
ddl_path_cdot <- panel %>% filter(grepl('Path', `Interpretation Summary`, ignore.case = T)) %>% select(`#Chr`, End, Ref, Alt, avsnp147) %>% mutate(pos_id=paste0(`#Chr`, ':', End, '_', Ref, '_', Alt))
##############
### UK10K ####
##############
metadata <- readxl::read_excel(path='~/git/EGA_EGAD00001002656_NGS_reanalyze/data/1-s2.0-S0002929716305274-mmc3.xlsx') %>% mutate(Sample=Patient)
sample_gene_comp_het <- metadata %>% filter(Status=='Solved'  & Variant_HGVSc!='NA' & GT=='0/1') %>% group_by(Sample, Gene) %>% summarise(Count=n()) %>% filter(Count>1) 
metadata <- left_join(metadata, sample_gene_comp_het) %>% mutate(Comp_Het_Path = case_when(Count >= 2 ~ 'CH', 
                                                                                           TRUE ~ 'No')) %>% 
  select(-Count)

allX <- raw_data %>% 
  mutate(DataSet_o=DataSet) %>%
  mutate(start = as.numeric(as.character(start))) %>% 
  mutate(Variant_genomic = paste0(chrom, ':', start + 1, ref, '>', alt)) %>% 
  mutate(DataSet = case_when(DataSet_o == 'ddl_nisc_100_panel' ~ 'DDL NISC RD Cohort',
                             DataSet_o == 'clinvar' & status != 'PATHOGENIC_OTHER' ~ 'ClinVar HC',
                             DataSet_o == 'clinvar' & status == 'PATHOGENIC_OTHER' ~ 'ClinVar LC Path',
                             DataSet_o == 'grimm' & source == 'humvar' ~ 'Grimm HumVar',
                             DataSet_o == 'grimm' & source == 'exovar' ~ 'Grimm ExoVar',
                             DataSet_o == 'grimm' & source == 'exovar__humvar' ~ 'Grimm ExoVar/HumVar',
                             DataSet_o == 'grimm' & source == 'predictSNP' ~ 'Grimm PredictSNP',
                             DataSet_o == 'grimm' & source == 'swissvar' ~ 'Grimm SwissVar',
                             DataSet_o == 'grimm' & source == 'varibench' ~ 'Grimm VariBench',
                             grepl('homsy', DataSet_o) ~ 'Homsy',
                             grepl('unifun', DataSet_o) ~ 'UniFun',
                             grepl('samocha', DataSet_o) ~ 'Samocha',
                             DataSet_o == 'gnomad' & (pos_id %in% model_data$ML_set__general_TT$train_set$pos_id ||
                                                        pos_id %in% model_data$ML_set__general_TT$test_set$pos_id ||
                                                        pos_id %in% model_data$ML_set__other_TT$train_set$pos_id ||
                                                        pos_id %in% model_data$ML_set__other_TT$test_set$pos_id) ~ 'gnomAD Benign',
                             DataSet_o == 'UK10K' ~ 'UK10K',
                             TRUE ~ 'Other')) %>%
  mutate(Status = case_when((DataSet == 'DDL NISC RD Cohort' & pos_id %in% ddl_path_cdot$pos_id) | 
                            (DataSet == 'DDL NISC RD Cohort' & end %in% ddl_path_cdot$End)  ~ 'Pathogenic',
                            DataSet == 'UK10K' & (Variant_genomic %in% (metadata %>% filter(Status == 'Solved') %>% pull(Variant_genomic))) ~ 'Pathogenic',
                            DataSet == 'UK10K' & (Variant_genomic %in% (metadata %>% filter(Status == 'Partially solved') %>% pull(Variant_genomic))) ~ 'Maybe Pathogenic',
                            grepl('pathogenic', DataSet_o, ignore.case = T) | grepl('pathogenic', status, ignore.case = T) ~ 'Pathogenic',
                            TRUE ~ 'NotPathogenic')) %>% 
  mutate(Status = case_when(Status = grepl('Grimm', DataSet) ~ status,
                            TRUE ~ Status)) %>% 
  mutate_at(vars(one_of(numeric_predictors)), funs(as.numeric(.))) %>% 
  #select(-status) %>% 
  filter(DataSet != 'Other') %>% 
  filter(Status != 'Maybe Pathogenic')

#rm(raw_data)
allX$fitcons_float <- allX$fitcons
allX[is.na(allX)] <- -1


# calculate VPaC scording for allX
allX$VPaC_m06 <- sqrt(predict(VPaC_6mtry, allX, type='prob')[,1])
allX$VPaC_m15 <- sqrt(predict(VPaC_15mtry, allX, type='prob')[,1])
allX$VPaC_m12 <- sqrt(predict(VPaC_12mtry, allX, type='prob')[,1])
allX$VPaC_m09 <- sqrt(predict(VPaC_9mtry, allX, type='prob')[,1])


###########################
## Deep LSTM ##
###########################

##############################
# keras model predict 
#############################
library(DMwR)
set.seed(89345)
train_sub <- model_data$ML_set__general_TT$train_set %>% dplyr::select(one_of(most_imp_predictors_no_disease_class), 'Status')
test_sub <- model_data$ML_set__general_TT$test_set %>% dplyr::select(one_of(most_imp_predictors_no_disease_class),'Status')

train_sub <- SMOTE(Status ~ ., as.data.frame(train_sub))
status_train <- train_sub$Status
status_train01 <- case_when(status_train == 'Pathogenic' ~ 1,
                            TRUE ~ 0)
status_test <- test_sub$Status
status_test01 <- case_when(status_test == 'Pathogenic' ~ 1,
                           TRUE ~ 0)
train_sub <- train_sub %>% dplyr::select(-Status)
test_sub <- test_sub %>% dplyr::select(-Status)
mean <- apply(train_sub %>% dplyr::select_(.dots = most_imp_predictors_no_disease_class), 2, mean)
std <- apply(train_sub %>% dplyr::select_(.dots = most_imp_predictors_no_disease_class), 2, sd)
train_data <- scale(train_sub, center=mean, scale=std)
test_data <- scale(test_sub, center=mean,scale=std)

# reshape
dim(train_data) <- c(nrow(train_data),1,length(most_imp_predictors_no_disease_class))
dim(test_data) <- c(nrow(test_data),1,length(most_imp_predictors_no_disease_class))

model <- keras_model_sequential() %>% 
  layer_lstm(1.5*length(most_imp_predictors_no_disease_class), recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1.5*length(most_imp_predictors_no_disease_class), recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1.5*length(most_imp_predictors_no_disease_class), recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1.5*length(most_imp_predictors_no_disease_class), recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1.5*length(most_imp_predictors_no_disease_class), recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1.5*length(most_imp_predictors_no_disease_class), recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1.5*length(most_imp_predictors_no_disease_class), recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1.5*length(most_imp_predictors_no_disease_class), recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class))) %>%
  layer_dropout(0.1) %>%
  layer_dense(units=1, activation='sigmoid')

model %>% compile(
  optimizer = optimizer_adam(),
  loss = 'binary_crossentropy',
  metric = c('accuracy')
)

history <- model %>% fit(train_data, status_train01, epochs = 10, batch_size=50)

DeepRNN <- list()
DeepRNN$model <- model
DeepRNN$mean <- mean
DeepRNN$std <- std
DeepRNN$predictors <- most_imp_predictors_no_disease_class
# save_model_hdf5(model, "/Volumes/data/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/DeepRNN_2018_06_22.h5")
# save(DeepRNN, file='/Volumes/data/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/DeepRNN_2018_06_22.Rdata')


scale_predict <- function(df, model, predictors, mean, std){
  # takes data frame, scales, predicts DeepRNN scoring,
  # and outputs DeepRNN scores
  set <- df %>% select(one_of(predictors))
  set_scale <- scale(set, center=mean, scale=std)
  dim(set_scale) <- c(nrow(set_scale),1,length(DeepRNN$predictors))
  preds <- model %>% predict(set_scale)
  preds <- preds[,1]
  preds[is.na(preds)] <- 0
  preds
}

########
# TEST #
########
test_set <- model_data$ML_set__general_TT$test_set %>% select(one_of(DeepRNN$predictors))
test_set$DeepRNN <- scale_predict(test_set, model, DeepRNN$predictors, DeepRNN$mean, DeepRNN$std)
# RF based prediction
test_set$VPaC_m06 <- sqrt(predict(VPaC_6mtry, test_set, type='prob')[,1])
test_set$VPaC_m15 <- sqrt(predict(VPaC_15mtry, test_set, type='prob')[,1])
test_set$VPaC_m12 <- sqrt(predict(VPaC_12mtry, test_set, type='prob')[,1])
test_set$VPaC_m09 <- sqrt(predict(VPaC_9mtry, test_set, type='prob')[,1])
test_set$Status <- model_data$ML_set__general_TT$test_set$Status
#########
# TRAIN #
#########
train_set <- model_data$ML_set__general_TT$train_set %>% select(one_of(DeepRNN$predictors))
train_set$DeepRNN <- scale_predict(train_set, model, DeepRNN$predictors, DeepRNN$mean, DeepRNN$std)
# RF based prediction
train_set$VPaC_m06 <- sqrt(predict(VPaC_6mtry, train_set, type='prob')[,1])
train_set$VPaC_m15 <- sqrt(predict(VPaC_15mtry, train_set, type='prob')[,1])
train_set$VPaC_m12 <- sqrt(predict(VPaC_12mtry, train_set, type='prob')[,1])
train_set$VPaC_m09 <- sqrt(predict(VPaC_9mtry, train_set, type='prob')[,1])
train_set$Status <- model_data$ML_set__general_TT$train_set$Status

test_set$VPaC_m15 <- sqrt(predict(VPaC_15mtry, test_set, type='prob')[,1])
test_set$DeepRNN <- keras_preds
test_set$Status <- model_data$ML_set__general_TT$test_set$Status

#############################
### create DeepVPaC score ###
#############################
fitControl_min <- trainControl(
  classProbs=T,
  savePredictions = T,
  allowParallel = T,
  summaryFunction = prSummary,
  returnData = T)
DeepVPaC <- caret::train(Status ~ ., data=test_set %>% select_(.dots=c('Status','VPaC_m15','DeepRNN')), 
                         method = "glm", metric='F', trControl=fitControl_min)

########################
# predict DeepVPaC on train/test
########################
test_set$DeepVPaC <- predict(DeepVPaC, test_set, type='prob')[,1]
train_set$DeepVPaC <- predict(DeepVPaC, train_set, type='prob')[,1]


# predict DeepVPaC on allX
# but first, scale data for DeepRNN model to work on
all_sub <- allX %>% select_(.dots=most_imp_predictors_no_disease_class)
all_sub$DeepRNN <- scale_predict(all_sub, model, DeepRNN$predictors, DeepRNN$mean, DeepRNN$std)
all_sub$VPaC_m15 <- allX$VPaC_m15
all_sub$DeepVPaC <- predict(DeepVPaC, all_sub, type='prob')[,1]

allX$DeepRNN <- all_sub$DeepRNN
allX$Distill <- all_sub$DeepVPaC

# merge test and train set with allX
allX2 <- bind_rows(allX, 
                   test_set %>% mutate(DataSet = 'Test Set', Distill = DeepVPaC), 
                   train_set %>% mutate(DataSet = 'Train Set', Distill = DeepVPaC))
allX2[is.na(allX2)] <- -1
allX <- allX2

save(allX, file='/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/allX_2018_06_22.Rdata')
