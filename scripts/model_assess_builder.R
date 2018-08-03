######################
# biowulf2 note
# run with:
# sinteractive --gres=gpu:k80:1,lscratch:10 --mem=96g -c8
# module load cuDNN/7.0/CUDA-9.0 CUDA R/3.5.0 python/3.5
# to load keras/tensorflow 
# https://hpc.nih.gov/apps/caret.html
######################

# create unified df for model_assess.Rmd

library(tidyverse)
library(caret)
library(ModelMetrics)
library(randomForest)
library(keras)
library(tensorflow)
library(xgboost)
library(PRROC)

# Load processed data and models

load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/model_data_2018_08_01.Rdata')

# load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_6mtry.Rdata')
# load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_6mtry_v8.Rdata')
# load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_9mtry_v8.Rdata')
load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_12mtry_v11
     .Rdata')
#VPaC <- VPaC
# all raw
load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/master/raw_data_2018_08_01.Rdata')





eyeIntegration <- c('eyeintegration_rnaseq_adipose_subcutaneous','eyeintegration_rnaseq_tpm_adipose_visceral_omentum','eyeintegration_rnaseq_tpm_adrenalgland','eyeintegration_rnaseq_tpm_artery_aorta','eyeintegration_rnaseq_tpm_artery_coronary','eyeintegration_rnaseq_tpm_artery_tibial','eyeintegration_rnaseq_tpm_brain_amygdala','eyeintegration_rnaseq_tpm_brain_anteriorcingulatecortex_ba24','eyeintegration_rnaseq_tpm_brain_caudate_basalganglia','eyeintegration_rnaseq_tpm_brain_cerebellarhemisphere','eyeintegration_rnaseq_tpm_brain_cerebellum','eyeintegration_rnaseq_tpm_brain_cortex','eyeintegration_rnaseq_tpm_brain_frontalcortex_ba9','eyeintegration_rnaseq_tpm_brain_hippocampus','eyeintegration_rnaseq_tpm_brain_hypothalamus','eyeintegration_rnaseq_tpm_brain_nucleusaccumbens_basalganglia','eyeintegration_rnaseq_tpm_brain_putamen_basalganglia','eyeintegration_rnaseq_tpm_brain_spinalcord_cervicalc_1','eyeintegration_rnaseq_tpm_brain_substantianigra','eyeintegration_rnaseq_tpm_breast_mammarytissue','eyeintegration_rnaseq_tpm_cells_ebv_transformedlymphocytes','eyeintegration_rnaseq_tpm_cells_transformedfibroblasts','eyeintegration_rnaseq_tpm_colon_sigmoid','eyeintegration_rnaseq_tpm_colon_transverse','eyeintegration_rnaseq_tpm_esc_stemcellline','eyeintegration_rnaseq_tpm_esophagus_gastroesophagealjunction','eyeintegration_rnaseq_tpm_esophagus_mucosa','eyeintegration_rnaseq_tpm_esophagus_muscularis','eyeintegration_rnaseq_tpm_heart_atrialappendage','eyeintegration_rnaseq_tpm_heart_leftventricle','eyeintegration_rnaseq_tpm_kidney_cortex','eyeintegration_rnaseq_tpm_liver','eyeintegration_rnaseq_tpm_lung','eyeintegration_rnaseq_tpm_minorsalivarygland','eyeintegration_rnaseq_tpm_muscle_skeletal','eyeintegration_rnaseq_tpm_nerve_tibial','eyeintegration_rnaseq_tpm_pancreas','eyeintegration_rnaseq_tpm_pituitary','eyeintegration_rnaseq_tpm_skin_notsunexposed_suprapubic','eyeintegration_rnaseq_tpm_skin_sunexposed_lowerleg','eyeintegration_rnaseq_tpm_smallintestine_terminalileum','eyeintegration_rnaseq_tpm_spleen','eyeintegration_rnaseq_tpm_stomach','eyeintegration_rnaseq_tpm_thyroid','eyeintegration_rnaseq_tpm_wholeblood')

numeric_predictors <- c('ccr_pct_v1','cadd_raw','vest3_rankscore','cadd_phred','mis_z','pli','lof_z','phylop_100way','revel','hapmap2','hapmap1','n_mis','epilogos_quies','n_lof','precessive','pnull','adj_exp_lof','adj_exp_syn','dann','adj_exp_mis','syn_z','n_syn','epilogos_txwk','fitcons','m_cap_score','m_cap_rankscore','eigen_phred','eigen_raw','epilogos_tx','is_lof','eigen_pc_raw_rankscore','epilogos_reprpcwk','fathmm_mkl_coding_rankscore','metalr_score','fathmm_mkl_coding_score','metalr_rankscore','impact_severity','metasvm_rankscore','metasvm_score','epilogos_enh','genocanyon_score','fathmm_converted_rankscore','mpc','epilogos_enhg','af_exac_all','epilogos_reprpc','max_aaf_all','mutationassessor_score','gerp','polyphen_score','gerp_elements','mutationassessor_score_rankscore','stam_mean','an_exac_all','af_exac_nfe','provean_converted_rankscore','an_exac_nfe','lrt_score','lrt_omega','grantham','lrt_converted_rankscore','genocanyon_score_rankscore','an_exac_afr','an_exac_amr','an_exac_sas','epilogos_het','ac_exac_all','linsight','gno_an_popmax','exac_num_het','an_exac_eas','gno_an_all','ac_exac_nfe','mutationtaster_converted_rankscore','an_exac_oth','an_exac_fin','gno_an_nfe','gno_af_all','gno_an_afr','epilogos_tssaflnk','gno_af_popmax','epilogos_znf','segway_sum_score','aaf_esp_ea','epilogos_txflnk','provean_score','segway_mean_score','epilogos_tss','aaf_esp_all','af_exac_amr','gno_af_nfe','epilogos_enhbiv','af_exac_sas','sift_score','fathmm_score','ac_exac_amr','aaf_esp_aa','gno_ac_all','gno_af_afr','ac_exac_sas','af_exac_eas','gno_an_fin','af_exac_afr','gno_an_eas','gno_an_oth','gno_ac_nfe','gno_ac_popmax','ac_exac_eas','ac_exac_afr','epilogos_tssbiv','gno_ac_afr','vest3_score','sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01')


nn_predictors <- c('ccr_pct_v1','cadd_raw','vest3_rankscore','cadd_phred','mis_z','pli','lof_z','phylop_100way','revel','hapmap2','hapmap1','n_mis','epilogos_quies','n_lof','precessive','pnull','adj_exp_lof','adj_exp_syn','dann','adj_exp_mis','syn_z','n_syn','epilogos_txwk','fitcons','m_cap_score','m_cap_rankscore','eigen_phred','eigen_raw','epilogos_tx','is_lof','eigen_pc_raw_rankscore','epilogos_reprpcwk','fathmm_mkl_coding_rankscore','metalr_score','fathmm_mkl_coding_score','metalr_rankscore','impact_severity','metasvm_rankscore','metasvm_score','epilogos_enh','genocanyon_score','fathmm_converted_rankscore','mpc','epilogos_enhg','af_exac_all','epilogos_reprpc','max_aaf_all','mutationassessor_score','gerp','polyphen_score','gerp_elements','mutationassessor_score_rankscore','stam_mean','an_exac_all','af_exac_nfe','provean_converted_rankscore','an_exac_nfe','lrt_score','lrt_omega','grantham','lrt_converted_rankscore','genocanyon_score_rankscore','an_exac_afr','an_exac_amr','an_exac_sas','epilogos_het','ac_exac_all','linsight','gno_an_popmax','exac_num_het','an_exac_eas','gno_an_all','ac_exac_nfe','mutationtaster_converted_rankscore','an_exac_oth','an_exac_fin','gno_an_nfe','gno_af_all','gno_an_afr','epilogos_tssaflnk','gno_af_popmax','epilogos_znf','segway_sum_score','aaf_esp_ea','epilogos_txflnk','provean_score','segway_mean_score','epilogos_tss','aaf_esp_all','af_exac_amr','gno_af_nfe','epilogos_enhbiv','af_exac_sas','sift_score','fathmm_score','ac_exac_amr','aaf_esp_aa','gno_ac_all','gno_af_afr','ac_exac_sas','af_exac_eas','gno_an_fin','af_exac_afr','gno_an_eas','gno_an_oth','gno_ac_nfe','gno_ac_popmax','ac_exac_eas','ac_exac_afr','epilogos_tssbiv','gno_ac_afr','vest3_score','sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01')

cm_maker <- function(predictor = 'cadd_phred', data, cutoff=0.5, mode = 'prec_recall') {
  if (class(predictor)!='character'){
    print("Running in predictor is a model mode")
    new_predictions <- predict(predictor, data, type='prob') %>% data.frame() %>% 
      mutate(Answers = data$Status, Prediction = case_when(Pathogenic > cutoff ~ 'Pathogenic', TRUE ~ 'NotPathogenic'))
    new_predictions <- new_predictions %>% mutate(preds = case_when(Prediction == 'Pathogenic' ~ 1,
                                                                    TRUE ~ 0),
                                                  actuals = case_when(Answers == 'Pathogenic' ~ 1,
                                                                      TRUE ~ 0))
    out <- caret::confusionMatrix(data = as.factor(new_predictions$Prediction), reference = as.factor(new_predictions$Answers), mode= mode)
    out$MCC <- mcc(new_predictions$preds, new_predictions$actuals, cutoff=cutoff)
  } else {
    print("Running in predictor is a precomputed column in data mode")
    new_predictions <- data 
    new_predictions$Prediction <- 'NotPathogenic'
    new_predictions[(new_predictions[,predictor] > cutoff), 'Prediction'] <- "Pathogenic"
    new_predictions <- new_predictions %>% mutate(preds = case_when(Prediction == 'Pathogenic' ~ 1,
                                                                    TRUE ~ 0),
                                                  actuals = case_when(Status == 'Pathogenic' ~ 1,
                                                                      TRUE ~ 0))
    out <- caret::confusionMatrix(data = as.factor(new_predictions$Prediction), reference = as.factor(new_predictions$Status), mode= mode)
    out$MCC <- mcc(new_predictions$preds, new_predictions$actuals, cutoff=cutoff)
  }
  out
}

pr_maker <- function(predictor, data, cutoff=0.5) {
  
  if (class(predictor)=='character'){
    predictor = enquo(predictor)
    pr.curve(scores.class0 = data %>% filter(Status=='Pathogenic') %>% pull(!!predictor),
             scores.class1 = data %>% filter(Status=='NotPathogenic') %>% pull(!!predictor),
             curve = T)
  }
  else {
    new_predictions <- predict(predictor, data, type = 'prob') %>% data.frame() %>% 
      mutate(Answers = data$Status, Prediction = case_when(Pathogenic > cutoff ~ 'Pathogenic', TRUE ~ 'NotPathogenic'))
    pr.curve(scores.class0 = new_predictions %>% filter(Answers=='Pathogenic') %>% pull(Pathogenic),
             scores.class1 = new_predictions %>% filter(Answers=='NotPathogenic') %>% pull(Pathogenic),
             curve = T)
  }
}

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
metadata <- left_join(metadata, sample_gene_comp_het) %>% 
  mutate(Comp_Het_Path = case_when(Count >= 2 ~ 'CH', 
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
#allX$VPaC_m06_v1 <- sqrt(predict(VPaC_6mtry, allX, type='prob')[,1])
allX$VPaC<- sqrt(predict(VPaC_12mtry_v11, allX, type='prob')[,1])
#allX$VPaC_m09_v8 <- sqrt(predict(VPaC_9mtry_v8, allX, type='prob')[,1])
#allX$VPaC_m06_v8 <- sqrt(predict(VPaC_6mtry_v8, allX, type='prob')[,1])

print('Data loaded')

###########################
## Deep LSTM ##
###########################

##############################
# keras model predict 
#############################
library(DMwR)
set.seed(89345)
############# 
# first remove correlated and near zero var predictors
# scratch this, redues F1 score
#############
#corred <- findCorrelation(model_data$ML_set__general_TT$train_set %>% dplyr::select(one_of(numeric_predictors)))
#numeric_predictors <- numeric_predictors[-corred]
#nZV <- nearZeroVar(model_data$ML_set__general_TT$train_set %>% dplyr::select(one_of(numeric_predictors)))
#numeric_predictors <- numeric_predictors[-nZV]

train_sub <- model_data$ML_set__general_TT$train_set %>% dplyr::select(one_of(nn_predictors),'Status') %>% 
  mutate_at(vars(one_of(nn_predictors)), funs(as.numeric(.)))
train_sub[is.na(train_sub)] <- -1
test_sub <- model_data$ML_set__general_TT$test_set %>% dplyr::select(one_of(nn_predictors),'Status') %>% 
  mutate_at(vars(one_of(nn_predictors)), funs(as.numeric(.)))
test_sub[is.na(test_sub)] <- -1

train_sub <- train_sub %>% sample_frac(1) # resample
train_sub <- SMOTE(Status ~ ., as.data.frame(train_sub))

status_train <- train_sub$Status
status_train01 <- case_when(status_train == 'Pathogenic' ~ 1,
                            TRUE ~ 0)
status_test <- test_sub$Status
status_test01 <- case_when(status_test == 'Pathogenic' ~ 1,
                           TRUE ~ 0)
train_sub <- train_sub %>% dplyr::select(-Status)
test_sub <- test_sub %>% dplyr::select(-Status)
mean <- apply(train_sub %>% dplyr::select_(.dots = nn_predictors), 2, mean)
std <- apply(train_sub %>% dplyr::select_(.dots = nn_predictors), 2, sd)
train_data <- scale(train_sub, center=mean, scale=std)
test_data <- scale(test_sub, center=mean,scale=std)

# reshape
dim(train_data) <- c(nrow(train_data),1,length(nn_predictors))
dim(test_data) <- c(nrow(test_data),1,length(nn_predictors))


model <- keras_model_sequential() %>% 
  layer_lstm(1.5*length(nn_predictors), recurrent_dropout=0.2, input_shape=c(1, length(nn_predictors)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1.5*length(nn_predictors), recurrent_dropout=0.2, input_shape=c(1, length(nn_predictors)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1.5*length(nn_predictors), recurrent_dropout=0.2, input_shape=c(1, length(nn_predictors)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1.5*length(nn_predictors), recurrent_dropout=0.2, input_shape=c(1, length(nn_predictors)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1.5*length(nn_predictors), recurrent_dropout=0.2, input_shape=c(1, length(nn_predictors)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1.5*length(nn_predictors), recurrent_dropout=0.2, input_shape=c(1, length(nn_predictors)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1.5*length(nn_predictors), recurrent_dropout=0.2, input_shape=c(1, length(nn_predictors)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1.5*length(nn_predictors), recurrent_dropout=0.2, input_shape=c(1, length(nn_predictors))) %>%
  layer_dropout(0.1) %>%
  layer_dense(units=1, activation='sigmoid')


model <- keras_model_sequential() %>% 
  layer_lstm(1.5*length(nn_predictors), recurrent_dropout=0.2, input_shape=c(1, length(nn_predictors)), return_sequences = T) %>%
  layer_dropout(0.7) %>%
  layer_lstm(1.5*length(nn_predictors), recurrent_dropout=0.2, input_shape=c(1, length(nn_predictors)), return_sequences = T) %>%
  layer_dropout(0.6) %>%
  layer_lstm(1.5*length(nn_predictors), recurrent_dropout=0.2, input_shape=c(1, length(nn_predictors)), return_sequences = T) %>%
  layer_dropout(0.5) %>%
  layer_lstm(1.5*length(nn_predictors), recurrent_dropout=0.2, input_shape=c(1, length(nn_predictors)), return_sequences = T) %>%
  layer_dropout(0.4) %>%
  layer_lstm(1.5*length(nn_predictors), recurrent_dropout=0.2, input_shape=c(1, length(nn_predictors)), return_sequences = T) %>%
  layer_dropout(0.3) %>%
  layer_lstm(1.5*length(nn_predictors), recurrent_dropout=0.2, input_shape=c(1, length(nn_predictors)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1.5*length(nn_predictors), recurrent_dropout=0.2, input_shape=c(1, length(nn_predictors))) %>%
  layer_dropout(0.1) %>%
  layer_dense(units=1, activation='sigmoid')



model %>% compile(
  optimizer = optimizer_adam(),
  loss = 'binary_crossentropy',
  metric = c('accuracy')
)

history <- model %>% fit(train_data, status_train01, epochs = 10, batch_size=50)

print('LSTM trained')

DeepRNN <- list()
DeepRNN$model <- model
DeepRNN$mean <- mean
DeepRNN$std <- std
DeepRNN$predictors <- nn_predictors
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

######################################
# xgboost
######################################
train_data <- model_data$ML_set__general_TT$train_set %>% 
  dplyr::select(one_of(numeric_predictors),'Status') %>% 
  mutate_at(vars(one_of(numeric_predictors)), funs(as.numeric(.)))
train_data[is.na(train_data)] <- -1
y <- recode(train_data$Status,'Pathogenic'=1, 'NotPathogenic'=0)
set.seed(35249)
xgbTree <- xgboost(label = y, 
                   eta = 0.4, 
                   max_depth = 3,
                   gamma = 0, 
                   colsample_bytree = 0.8,
                   min_child_weight = 1, 
                   subsample = 0.75,
                   data = train_data %>% select_if(is.numeric) %>% as.matrix(), 
                   nrounds = 100, 
                   objective = "binary:logistic", 
                   eval_metric = 'aucpr', 
                   nthread = 16)

print('xgboost trained')

########
# TEST #
########
test_set <- model_data$ML_set__general_TT$test_set %>% 
  mutate_at(vars(one_of(c(nn_predictors,numeric_predictors))), funs(as.numeric(.)))
test_set[is.na(test_set)] <- -1

test_set$DeepRNN <- scale_predict(test_set, model, DeepRNN$predictors, DeepRNN$mean, DeepRNN$std)
# RF based prediction
test_set$fitcons_float <- test_set$fitcons
#test_set$VPaC_m06_v1 <- sqrt(predict(VPaC_6mtry, test_set, type='prob')[,1])
test_set$VPaC <- sqrt(predict(VPaC_12mtry_v11, test_set, type='prob')[,1])
#test_set$VPaC_m09_v8 <- sqrt(predict(VPaC_9mtry_v8, test_set, type='prob')[,1])
#test_set$VPaC_m06_v8 <- sqrt(predict(VPaC_6mtry_v8, test_set, type='prob')[,1])
test_set$xgbTree <- sqrt(predict(xgbTree, test_set %>% dplyr::select(one_of(numeric_predictors)) %>% as.matrix()))


#########
# TRAIN #
#########
train_set <- model_data$ML_set__general_TT$train_set %>% 
  mutate_at(vars(one_of(c(nn_predictors,numeric_predictors))), funs(as.numeric(.)))
train_set[is.na(train_set)] <- -1

train_set$DeepRNN <- scale_predict(train_set, model, DeepRNN$predictors, DeepRNN$mean, DeepRNN$std)
# RF based prediction
train_set$fitcons_float <- train_set$fitcons
#train_set$VPaC_m06_v1 <- sqrt(predict(VPaC_6mtry, train_set, type='prob')[,1])
train_set$VPaC <- sqrt(predict(VPaC_12mtry_v11, train_set, type='prob')[,1])
#train_set$VPaC_m09_v8 <- sqrt(predict(VPaC_9mtry_v8, train_set, type='prob')[,1])
#train_set$VPaC_m06_v8 <- sqrt(predict(VPaC_6mtry_v8, train_set, type='prob')[,1])
train_set$xgbTree <- sqrt(predict(xgbTree, train_set %>% dplyr::select(one_of(numeric_predictors)) %>% as.matrix()))


#########
# OTHER #
#########
other_set <- bind_rows(model_data$ML_set__other_TT$train_set %>% 
                         mutate_at(vars(one_of(c(nn_predictors,numeric_predictors))), funs(as.numeric(.))),
                       model_data$ML_set__other_TT$test_set %>% 
                         mutate_at(vars(one_of(c(nn_predictors,numeric_predictors))), funs(as.numeric(.))))

other_set[is.na(other_set)] <- -1

other_set$DeepRNN <- scale_predict(other_set, model, DeepRNN$predictors, DeepRNN$mean, DeepRNN$std)
# RF based prediction
other_set$fitcons_float <- other_set$fitcons
#other_set$VPaC_m06_v1 <- sqrt(predict(VPaC_6mtry, other_set, type='prob')[,1])
other_set$VPaC <- sqrt(predict(VPaC_12mtry_v11, other_set, type='prob')[,1])
#other_set$VPaC_m09_v8 <- sqrt(predict(VPaC_9mtry_v8, other_set, type='prob')[,1])
#other_set$VPaC_m06_v8 <- sqrt(predict(VPaC_6mtry_v8, other_set, type='prob')[,1])
other_set$xgbTree <- sqrt(predict(xgbTree, other_set %>% dplyr::select(one_of(numeric_predictors)) %>% as.matrix()))

other_set$Status <- c(as.character(model_data$ML_set__other_TT$train_set$Status), 
                      as.character(model_data$ML_set__other_TT$test_set$Status))



#########
# UK10K Withheld #
#########
withheld_set <- model_data$Test_set__UK10K %>% 
  mutate_at(vars(one_of(c(nn_predictors,numeric_predictors))), funs(as.numeric(.)))
withheld_set[is.na(withheld_set)] <- -1

withheld_set$DeepRNN <- scale_predict(withheld_set, model, DeepRNN$predictors, DeepRNN$mean, DeepRNN$std)
# RF based prediction
withheld_set$fitcons_float <- withheld_set$fitcons
#withheld_set$VPaC_m06_v1 <- sqrt(predict(VPaC_6mtry, withheld_set, type='prob')[,1])
withheld_set$VPaC <- sqrt(predict(VPaC_12mtry_v11, withheld_set, type='prob')[,1])
#withheld_set$VPaC_m09_v8 <- sqrt(predict(VPaC_9mtry_v8, withheld_set, type='prob')[,1])
#withheld_set$VPaC_m06_v8 <- sqrt(predict(VPaC_6mtry_v8, withheld_set, type='prob')[,1])
withheld_set$xgbTree <- sqrt(predict(xgbTree, withheld_set %>% dplyr::select(one_of(numeric_predictors)) %>% as.matrix()))

withheld_set$Status <- as.character(model_data$Test_set__UK10K$Status)

print('test, train, other, uk10k withheld calculated')


#############################
### create DeepVPaC score ###
#############################

# split test into tune/test
set.seed(783)
tune_set <- test_set %>% sample_frac(0.5)
test_setN <- test_set %>% filter(!pos_id %in% tune_set$pos_id)
test_set <- test_setN

# fitControl_min <- trainControl(
#   classProbs=T,
#   savePredictions = T,
#   allowParallel = T,
#   summaryFunction = prSummary,
#   returnData = T)

# library(parallel)
# library(doParallel)
# cluster <- makeCluster(16) 
# registerDoParallel(cluster)

# use grid search to linearly blend VPaC, DeepRNN, and xgbTree on the tune set
# build naive
params <- expand.grid(seq(0,0.3, by = 0.01), seq(0.1,0.6, by = 0.01)) %>% data.frame()
params$Var3 = 1-(params$Var1 + params$Var2)                
#params

mcc <- c()
for (i in seq(1:nrow(params))){
  tune_set$grid <- (tune_set$DeepRNN * params[i,1]) + 
    (tune_set$VPaC * params[i,2]) +
    (tune_set$xgbTree * params[i,3])
  cm_out <- cm_maker('grid', tune_set)
  mcc <- c(mcc, cm_out$MCC)
  #print(paste(params[i,], cm_out$MCC))
}
params$MCC <- mcc
print(params %>% arrange(-mcc) %>% head(20))
#######
# DeepRNN * (params %>% arrange(-mcc) %>% head(1))[1] + VPaC * (params %>% arrange(-mcc) %>% head(1))[2] + xgbTree * (params %>% arrange(-mcc) %>% head(1))[3]
#######

########################
# Calculate DeepVPaC on train/test/other
########################


# but predict xgbTree, then DeepRNN with scale_predict on AllX
allX <- allX %>% mutate_at(vars(one_of(c(nn_predictors,numeric_predictors))), funs(as.numeric(.)))
allX[is.na(allX)] <- -1
allX$xgbTree <- sqrt(predict(xgbTree, allX %>% dplyr::select(one_of(numeric_predictors)) %>% as.matrix()))

all_sub <- allX
all_sub$DeepRNN <- scale_predict(all_sub, model, DeepRNN$predictors, DeepRNN$mean, DeepRNN$std)

allX$DeepRNN <- all_sub$DeepRNN
# merge test and train set with allX
allX2 <- bind_rows(allX %>% mutate_all(as.character), 
                   withheld_set %>% mutate(DataSet = 'UK10K Withheld') %>% mutate_all(as.character),
                   tune_set %>% mutate(DataSet = 'Tune Set') %>% mutate_all(as.character),
                   test_set %>% mutate(DataSet = 'Test Set') %>% mutate_all(as.character), 
                   train_set %>% mutate(DataSet = 'Train Set') %>% mutate_all(as.character),
                   other_set %>% mutate(DataSet = 'Other Set') %>% mutate_all(as.character)) %>% 
  mutate_at(vars(one_of(c(numeric_predictors, nn_predictors))), funs(as.numeric(.)))
allX2[is.na(allX2)] <- -1
allX <- allX2

allX <- allX %>% 
  mutate_at(vars(contains('VP')), as.numeric) %>% 
  mutate_at(vars(contains('Deep')), as.numeric) %>% 
  mutate_at(vars(contains('xgb')), as.numeric) %>% 
  mutate(Status=factor(Status, levels=c('Pathogenic','NotPathogenic')))

print('allX reformed')
#############
# build assess data 
#############

SuperGrimm <- allX %>% filter(grepl('Grimm', DataSet)) %>% 
  filter(!pos_id %in% (allX %>% filter(!grepl('Grimm', DataSet)) %>% 
                         pull(pos_id)))

HalfOther_1 <- allX %>% 
  filter(DataSet == 'Other Set', Status=='NotPathogenic') %>% 
  filter(!pos_id %in% (allX %>% filter(DataSet != 'Other Set') %>% pull(pos_id))) %>% 
  mutate(DataSet = 'Other Set') %>% 
  sample_frac(0.5)

HalfOther_2 <- allX %>% 
  filter(DataSet == 'Other Set') %>% 
  filter(!pos_id %in% (allX %>% filter(DataSet != 'Other Set') %>% pull(pos_id))) %>% 
  filter(!pos_id %in% HalfOther_1)

assess_set <- bind_rows(SuperGrimm %>% mutate(DataSet = 'SuperGrimm'), 
                        HalfOther_1 %>% mutate(DataSet = 'SuperGrimm'),
                        HalfOther_2 %>% mutate(DataSet = 'Extra NotPathogenic'),
                        allX %>% filter(DataSet == 'Test Set'),
                        allX %>% filter(DataSet == 'DDL NISC RD Cohort'),
                        allX %>% filter(DataSet == 'Unifun'),
                        allX %>% filter(DataSet == 'Homsy'),
                        allX %>% filter(DataSet == 'Samocha'),
                        allX %>% filter(DataSet == 'UK10K Withheld') %>% filter(!pos_id %in% c(train_set$pos_id, tune_set$pos_id)))
#allX %>% filter(DataSet == 'UK10K') %>% 
#  filter(!pos_id %in% (model_data$ML_set__general_TT$train_set$pos_id)) %>% 
#  filter(!pos_id %in% ((model_data$ML_set__general_TT$tune_set$pos_id))))
print('assess made')


### calculate DIstill on AllX and Assess Set
allX$Distill <- (allX$DeepRNN * (params %>% arrange(-mcc) %>% head(1))[1] %>% as.numeric()) + 
  (allX$VPaC * (params %>% arrange(-mcc) %>% head(1))[2] %>% as.numeric()) + 
  (allX$xgbTree * (params %>% arrange(-mcc) %>% head(1))[3] %>% as.numeric())

assess_set$Distill <- (assess_set$DeepRNN * (params %>% arrange(-mcc) %>% head(1))[1] %>% as.numeric()) + 
  (assess_set$VPaC * (params %>% arrange(-mcc) %>% head(1))[2] %>% as.numeric()) + 
  (assess_set$xgbTree * (params %>% arrange(-mcc) %>% head(1))[3] %>% as.numeric())

###
save(allX, file='/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/allX_2018_08_02.Rdata')
save(assess_set, file='/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/assess_2018_08_02.Rdata')
