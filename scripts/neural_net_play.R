load('/Volumes/data/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/model_data_2018_06_19.Rdata')
most_imp_predictors_no_disease_class <- c('is_lof','impact_severity','mis_z','ccr_pct_v1','cadd_phred','phylop_100way','n_mis','revel','fitcons','precessive','n_lof','m_cap_rankscore','dann','vest3_rankscore','n_syn','pnull','pli','lof_z','fathmm_mkl_coding_rankscore','an_exac_all','eigen_pc_raw_rankscore','gerp_elements','mutationassessor_score_rankscore','mpc','metasvm_rankscore','polyphen_score','metalr_rankscore','lrt_converted_rankscore','genocanyon_score_rankscore','mutationtaster_converted_rankscore','gno_an_popmax','grantham','max_aaf_all','ac_exac_all','fathmm_converted_rankscore','aaf_esp_all','sift_score','ac_exac_sas','linsight',                          
                                          # 'exon',
                                           'aa_length',
                                           'cpg_island',
                                          'epilogos_bivflnk',
                                          'epilogos_enh',
                                          'epilogos_enhbiv',
                                          'epilogos_enhg',
                                          'epilogos_het',
                                          'epilogos_quies',
                                          'epilogos_reprpc',
                                          'epilogos_reprpcwk',
                                          'epilogos_tss',
                                          'epilogos_tssaflnk',
                                          'epilogos_tssbiv',
                                          'epilogos_tx',
                                          'epilogos_txflnk',
                                          'epilogos_txwk',
                                          'epilogos_znf'
                                          )

most_imp_predictors_no_disease_class <- c('ccr_pct_v1',
                                          'cadd_raw',
                                          'vest3_rankscore',
                                          'cadd_phred',
                                          'mis_z',
                                          'pli',
                                          'lof_z',
                                          'phylop_100way',
                                          'revel',
                                          'hapmap2',
                                          'hapmap1',
                                          'n_mis',
                                          'epilogos_quies',
                                          'n_lof',
                                          'precessive',
                                          'pnull',
                                          'adj_exp_lof',
                                          'adj_exp_syn',
                                          'dann',
                                          'adj_exp_mis',
                                          'syn_z',
                                          'n_syn',
                                          'epilogos_txwk',
                                          'fitcons',
                                          'm_cap_score',
                                          'm_cap_rankscore',
                                          'eigen_phred',
                                          'eigen_raw',
                                          'epilogos_tx',
                                          'is_lof',
                                          'eigen_pc_raw_rankscore',
                                          'epilogos_reprpcwk',
                                          'fathmm_mkl_coding_rankscore',
                                          'metalr_score',
                                          'fathmm_mkl_coding_score',
                                          'metalr_rankscore',
                                          'impact_severity',
                                          'metasvm_rankscore',
                                          'metasvm_score',
                                          'epilogos_enh',
                                          'genocanyon_score',
                                          'fathmm_converted_rankscore',
                                          'mpc',
                                          'epilogos_enhg',
                                          'af_exac_all',
                                          'epilogos_reprpc',
                                          'max_aaf_all',
                                          'mutationassessor_score',
                                          'gerp',
                                          'polyphen_score',
                                          'gerp_elements',
                                          'mutationassessor_score_rankscore',
                                          'stam_mean',
                                          'an_exac_all',
                                          'af_exac_nfe',
                                          'provean_converted_rankscore',
                                          'an_exac_nfe',
                                          'lrt_score',
                                          'lrt_omega',
                                          'grantham',
                                          'lrt_converted_rankscore',
                                          'genocanyon_score_rankscore',
                                          'an_exac_afr',
                                          'an_exac_amr',
                                          'an_exac_sas',
                                          'epilogos_het',
                                          'ac_exac_all',
                                          'linsight',
                                          'gno_an_popmax',
                                          'exac_num_het',
                                          'an_exac_eas',
                                          'gno_an_all',
                                          'ac_exac_nfe',
                                          'mutationtaster_converted_rankscore',
                                          'an_exac_oth',
                                          'an_exac_fin',
                                          'gno_an_nfe',
                                          'gno_af_all',
                                          'gno_an_afr',
                                          'epilogos_tssaflnk',
                                          'gno_af_popmax',
                                          'epilogos_znf',
                                          'segway_sum_score',
                                          'aaf_esp_ea',
                                          'epilogos_txflnk',
                                          'provean_score',
                                          'segway_mean_score',
                                          'epilogos_tss',
                                          'aaf_esp_all',
                                          'af_exac_amr',
                                          'gno_af_nfe',
                                          'epilogos_enhbiv',
                                          'af_exac_sas',
                                          'sift_score',
                                          'fathmm_score',
                                          'ac_exac_amr',
                                          'aaf_esp_aa',
                                          'gno_ac_all',
                                          'gno_af_afr',
                                          'ac_exac_sas',
                                          'af_exac_eas',
                                          'gno_an_fin',
                                          'af_exac_afr',
                                          'gno_an_eas',
                                          'gno_an_oth',
                                          'gno_ac_nfe',
                                          'gno_ac_popmax',
                                          'ac_exac_eas',
                                          'ac_exac_afr',
                                          'epilogos_tssbiv',
                                          'gno_ac_afr',
                                          'vest3_score')

library(caret)
library(ModelMetrics)
library(randomForest)
library(tidyverse)
library(keras)
library(tensorflow)
# confusion matrix maker for built models
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

##################
#custom  metrics
# which don't seem to be useful
# sensitivity <- function(y_true, y_pred){
#   true_positives = k_sum(k_round(k_clip(y_true * y_pred, 0, 1)))
#   possible_positives = k_sum(k_round(k_clip(y_true, 0, 1)))
#   true_positives / (possible_positives + k_epsilon())
# }
# 
# specificity <- function(y_true, y_pred){
#   true_negatives = k_sum(k_round(k_clip((1-y_true) * (1-y_pred), 0, 1)))
#   possible_negatives = k_sum(k_round(k_clip(1-y_true, 0, 1)))
#   true_negatives / (possible_negatives + k_epsilon())
# }
# 
# precision <- function(y_true, y_pred){
#   true_positives = k_sum(k_round(k_clip(y_true * y_pred, 0, 1)))
#   predicted_positives = k_sum(k_round(k_clip(y_pred, 0, 1)))
#   precision = true_positives / (predicted_positives + k_epsilon())
#   precision
# }
# 
# recall <- function(y_true, y_pred){
#   true_positives = k_sum(k_round(k_clip(y_true * y_pred, 0, 1)))
#   possible_positives = k_sum(k_round(k_clip(y_true, 0, 1)))
#   recall = true_positives / (possible_positives + k_epsilon())
#   recall
# }
# 
# 
# f1 <- function(y_true, y_pred){
#   true_positives = k_sum(k_round(k_clip(y_true * y_pred, 0, 1)))
#   predicted_positives = k_sum(k_round(k_clip(y_pred, 0, 1)))
#   precision = true_positives / (predicted_positives + k_epsilon())
# 
#   true_positives = k_sum(k_round(k_clip(y_true * y_pred, 0, 1)))
#   possible_positives = k_sum(k_round(k_clip(y_true, 0, 1)))
#   recall = true_positives / (possible_positives + k_epsilon())
# 
#   f1 = 2 * ((precision * recall)/(precision + recall))
#   f1
# }
###############
# multiperceptron
# doesn't work too well
###############

library(tidyverse)
# keras tensorflow deep learning
library(keras)


set.seed(89345)
train_sub <- model_data$ML_set__general_dummy_TT$train_set %>% group_by(Status) %>% group_by(Status) %>% sample_n(887) %>% ungroup() %>% select_(.dots=most_imp_predictors_no_disease_class, 'Status')
test_sub <- model_data$ML_set__general_dummy_TT$test_set %>% select_(.dots=most_imp_predictors_no_disease_class,'Status')
status_train <- train_sub$Status
status_train01 <- case_when(status_train == 'Pathogenic' ~ 1,
                            TRUE ~ 0)
status_test <- test_sub$Status
status_test01 <- case_when(status_test == 'Pathogenic' ~ 1,
                           TRUE ~ 0)
train_sub <- train_sub %>% select(-Status)
test_sub <- test_sub %>% select(-Status)
mean <- apply(train_sub %>% select_(.dots = most_imp_predictors_no_disease_class), 2, mean)
std <- apply(train_sub %>% select_(.dots = most_imp_predictors_no_disease_class), 2, sd)
train_data <- scale(train_sub, center=mean, scale=std)
test_data <- scale(test_sub, center=mean,scale=std)


model <- keras_model_sequential() %>% 
  layer_dense(units = 12*length(most_imp_predictors_no_disease_class), activation = 'relu', input_shape=c(length(most_imp_predictors_no_disease_class))) %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 12*length(most_imp_predictors_no_disease_class), activation = 'relu') %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 12*length(most_imp_predictors_no_disease_class), activation = 'relu') %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 12*length(most_imp_predictors_no_disease_class), activation = 'relu') %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 12*length(most_imp_predictors_no_disease_class), activation = 'relu') %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units=1, activation='sigmoid')  


model %>% compile(
  optimizer = 'rmsprop',
  loss = 'binary_crossentropy',
  metric = c('accuracy')
)

history <- model %>% fit(train_data, status_train01, epochs = 25, batch_size=50, validation_data = list(test_data, status_test01))

validation_score <- model %>% evaluate(test_data, status_test01)
test_score <- model %>% predict(test_data)
test_score[is.na(test_score)] <- 0
model_data$ML_set__general_dummy_TT$test_set$keras <- test_score[,1]
cm_maker('keras', model_data$ML_set__general_dummy_TT$test_set, cutoff=0.999)

###############
# SMOTE 
###############
library(DMwR)
set.seed(89345)
train_sub <- model_data$ML_set__general_dummy_TT$train_set %>% select_(.dots=most_imp_predictors_no_disease_class, 'Status')
test_sub <- model_data$ML_set__general_dummy_TT$test_set %>% select_(.dots=most_imp_predictors_no_disease_class,'Status')

train_sub <- SMOTE(Status ~ ., as.data.frame(train_sub))
status_train <- train_sub$Status
status_train01 <- case_when(status_train == 'Pathogenic' ~ 1,
                            TRUE ~ 0)
status_test <- test_sub$Status
status_test01 <- case_when(status_test == 'Pathogenic' ~ 1,
                           TRUE ~ 0)
train_sub <- train_sub %>% select(-Status)
test_sub <- test_sub %>% select(-Status)
mean <- apply(train_sub %>% select_(.dots = most_imp_predictors_no_disease_class), 2, mean)
std <- apply(train_sub %>% select_(.dots = most_imp_predictors_no_disease_class), 2, sd)
train_data <- scale(train_sub, center=mean, scale=std)
test_data <- scale(test_sub, center=mean,scale=std)


model <- keras_model_sequential() %>% 
  layer_dense(units = 400, activation = 'relu', input_shape=c(length(most_imp_predictors_no_disease_class))) %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 400, activation = 'relu') %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units=1, activation='sigmoid')  


model %>% compile(
  optimizer = 'rmsprop',
  loss = 'binary_crossentropy',
  metric = c('accuracy')
)

history <- model %>% fit(train_data, status_train01, epochs = 25, batch_size=50, validation_data = list(test_data, status_test01))

validation_score <- model %>% evaluate(test_data, status_test01)
test_score <- model %>% predict(test_data)
test_score[is.na(test_score)] <- 0
model_data$ML_set__general_dummy_TT$test_set$keras <- test_score[,1]
cm_maker('keras', model_data$ML_set__general_dummy_TT$test_set, cutoff=0.9)


###############
# rnn deep
# upsample path
# gives highest f1 score (~0.41 for test data)
# doesn't merge too well when ensembling with VPaC RF model
###############
set.seed(89345)
train_sub <- model_data$ML_set__general_dummy_TT$train_set  %>% select_(.dots=most_imp_predictors_no_disease_class, 'Status') 
test_sub <- model_data$ML_set__general_dummy_TT$test_set %>% select_(.dots=most_imp_predictors_no_disease_class,'Status')

train_sub <- upSample(as.data.frame(train_sub), train_sub$Status)
status_train <- train_sub$Class
status_train01 <- case_when(status_train == 'Pathogenic' ~ 1,
                            TRUE ~ 0)
status_test <- test_sub$Status
status_test01 <- case_when(status_test == 'Pathogenic' ~ 1,
                           TRUE ~ 0)
train_sub <- train_sub %>% select(-Class, -Status)
test_sub <- test_sub %>% select(-Status)
mean <- apply(train_sub %>% select_(.dots = most_imp_predictors_no_disease_class), 2, mean)
std <- apply(train_sub %>% select_(.dots = most_imp_predictors_no_disease_class), 2, sd)
train_data <- scale(train_sub, center=mean, scale=std)
test_data <- scale(test_sub, center=mean,scale=std)

# reshape
dim(train_data) <- c(nrow(train_data),1,length(most_imp_predictors_no_disease_class))
dim(test_data) <- c(nrow(test_data),1,length(most_imp_predictors_no_disease_class))




model <- keras_model_sequential() %>% 
  layer_lstm(196, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.4) %>%
  layer_lstm(196, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.3) %>%
  layer_lstm(196, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(196, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class))) %>%
  layer_dropout(0.1) %>%
  layer_dense(units=1, activation='sigmoid')



model %>% compile(
  optimizer = 'rmsprop',
  loss = 'binary_crossentropy',
  metric = c('accuracy')
)

history <- model %>% fit(train_data, status_train01, epochs = 10, batch_size=5000)

#validation_score <- model %>% evaluate(test_data, status_test01)
test_score <- model %>% predict(test_data)
test_score[is.na(test_score)] <- 0
model_data$ML_set__general_dummy_TT$test_set$keras <- test_score[,1]
cm_maker('keras', model_data$ML_set__general_dummy_TT$test_set, cutoff=0.99)

####
# save
####
DeepRNN  <- model
save(DeepRNN, file = '/Volumes/data/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/DeepRNN_model.Rdata')


###############
# rnn smote 1
# adds new minor class (path) and trim down major class
# currently this is the model that works the best when ensembl-ing with RF-based VPaC
###############
library(DMwR)
set.seed(89345)
train_sub <- model_data$ML_set__general_dummy_TT$train_set %>% select(one_of(most_imp_predictors_no_disease_class), 'Status')
test_sub <- model_data$ML_set__general_dummy_TT$test_set %>% select(one_of(most_imp_predictors_no_disease_class),'Status')
                      
train_sub <- SMOTE(Status ~ ., as.data.frame(train_sub))
status_train <- train_sub$Status
status_train01 <- case_when(status_train == 'Pathogenic' ~ 1,
                                                  TRUE ~ 0)
status_test <- test_sub$Status
status_test01 <- case_when(status_test == 'Pathogenic' ~ 1,
                                                 TRUE ~ 0)
train_sub <- train_sub %>% select(-Status)
test_sub <- test_sub %>% select(-Status)
mean <- apply(train_sub %>% select_(.dots = most_imp_predictors_no_disease_class), 2, mean)
std <- apply(train_sub %>% select_(.dots = most_imp_predictors_no_disease_class), 2, sd)
train_data <- scale(train_sub, center=mean, scale=std)
test_data <- scale(test_sub, center=mean,scale=std)
                      
# reshape
dim(train_data) <- c(nrow(train_data),1,length(most_imp_predictors_no_disease_class))
dim(test_data) <- c(nrow(test_data),1,length(most_imp_predictors_no_disease_class))
                      
model <- keras_model_sequential() %>% 
  layer_lstm(1*length(most_imp_predictors_no_disease_class), recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1*length(most_imp_predictors_no_disease_class), recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1*length(most_imp_predictors_no_disease_class), recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1*length(most_imp_predictors_no_disease_class), recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1*length(most_imp_predictors_no_disease_class), recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1*length(most_imp_predictors_no_disease_class), recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1*length(most_imp_predictors_no_disease_class), recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(1*length(most_imp_predictors_no_disease_class), recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class))) %>%
  layer_dropout(0.1) %>%
  layer_dense(units=1, activation='sigmoid')



model %>% compile(
  optimizer = optimizer_adam(),
  loss = 'binary_crossentropy',
  metric = c('accuracy')
)

history <- model %>% fit(train_data, status_train01, epochs = 10, batch_size=50)

#validation_score <- model %>% evaluate(test_data, status_test01)
test_score <- model %>% predict(test_data)
test_score[is.na(test_score)] <- 0
model_data$ML_set__general_dummy_TT$test_set$keras <- test_score[,1]
cm_maker('keras', model_data$ML_set__general_dummy_TT$test_set, cutoff=0.95)



###############
# rnn smote 2
# adds new minor class (path) and trim down major class
# screwing with dim and depth
###############


model <- keras_model_sequential() %>% 
  layer_lstm(96, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(96, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(96, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(96, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class))) %>%
  layer_dropout(0.1) %>%
  layer_dense(units=1, activation='sigmoid')



model %>% compile(
  optimizer = optimizer_adam( lr= 0.01),
  loss = 'binary_crossentropy',
  metric = c('accuracy')
)

history <- model %>% fit(train_data, status_train01, epochs = 10, batch_size=50)

#validation_score <- model %>% evaluate(test_data, status_test01)
test_score <- model %>% predict(test_data)
test_score[is.na(test_score)] <- 0
model_data$ML_set__general_dummy_TT$test_set$keras <- test_score[,1]
cm_maker('keras', model_data$ML_set__general_dummy_TT$test_set, cutoff=0.99)





###############################
# pull grimm path
# rnn smote
###############################

more <- all %>% filter(grepl('grimm',DataSet, ignore.case = T)) 
set.seed(1783) 
half_more <- more %>% group_by(Status) %>% sample_frac(0.5)


###############
library(DMwR)
set.seed(89345)
train_sub <- model_data$ML_set__general_dummy_TT$train_set %>% select_(.dots=most_imp_predictors_no_disease_class, 'Status')
train_sub <- bind_rows(train_sub, half_more %>% select_(.dots=most_imp_predictors_no_disease_class, 'Status'))
test_sub <- model_data$ML_set__general_dummy_TT$test_set %>% select_(.dots=most_imp_predictors_no_disease_class,'Status')

train_sub <- SMOTE(Status ~ ., as.data.frame(train_sub))
status_train <- train_sub$Status
status_train01 <- case_when(status_train == 'Pathogenic' ~ 1,
                            TRUE ~ 0)
status_test <- test_sub$Status
status_test01 <- case_when(status_test == 'Pathogenic' ~ 1,
                           TRUE ~ 0)
train_sub <- train_sub %>% select(-Status)
test_sub <- test_sub %>% select(-Status)
mean <- apply(train_sub %>% select_(.dots = most_imp_predictors_no_disease_class), 2, mean)
std <- apply(train_sub %>% select_(.dots = most_imp_predictors_no_disease_class), 2, sd)
train_data <- scale(train_sub, center=mean, scale=std)
test_data <- scale(test_sub, center=mean,scale=std)

# reshape
dim(train_data) <- c(nrow(train_data),1,length(most_imp_predictors_no_disease_class))
dim(test_data) <- c(nrow(test_data),1,length(most_imp_predictors_no_disease_class))

model <- keras_model_sequential() %>% 
  layer_lstm(196, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.5) %>%
  layer_lstm(196, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.4) %>%
  layer_lstm(196, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.3) %>%
  layer_lstm(196, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(196, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class))) %>%
  layer_dropout(0.1) %>%
  layer_dense(units=1, activation='sigmoid')



model %>% compile(
  optimizer = 'rmsprop',
  loss = 'binary_crossentropy',
  metric = c('accuracy')
)

history <- model %>% fit(train_data, status_train01, epochs = 10, batch_size=5000)

#validation_score <- model %>% evaluate(test_data, status_test01)
test_score <- model %>% predict(test_data)
test_score[is.na(test_score)] <- 0
model_data$ML_set__general_dummy_TT$test_set$keras <- test_score[,1]
cm_maker('keras', model_data$ML_set__general_dummy_TT$test_set, cutoff=0.99)



