load('/Volumes/data/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/model_data.Rdata')
most_imp_predictors_no_disease_class <- c('is_lof','impact_severity','mis_z','ccr_pct_v1','cadd_phred','phylop_100way','n_mis','revel','fitcons_float','precessive','n_lof','m_cap_rankscore','dann','vest3_rankscore','n_syn','pnull','pli','lof_z','fathmm_mkl_coding_rankscore','an_exac_all','eigen_pc_raw_rankscore','gerp_elements','mutationassessor_score_rankscore','mpc','metasvm_rankscore','polyphen_score','metalr_rankscore','lrt_converted_rankscore','genocanyon_score_rankscore','mutationtaster_converted_rankscore','gno_an_popmax','grantham','max_aaf_all','ac_exac_all','fathmm_converted_rankscore','aaf_esp_all','sift_score','ac_exac_sas','linsight')


library(caret)
library(ModelMetrics)
library(randomForest)
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
    out <- caret::confusionMatrix(data = new_predictions$Prediction, reference = new_predictions$Answers, mode= mode)
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
    out <- caret::confusionMatrix(data = new_predictions$Prediction, reference = new_predictions$Status, mode= mode)
    out$MCC <- mcc(new_predictions$preds, new_predictions$actuals, cutoff=cutoff)
  }
  out
}

##################
#custom  metric
sensitivity <- function(y_true, y_pred){
  true_positives = k_sum(k_round(k_clip(y_true * y_pred, 0, 1)))
  possible_positives = k_sum(k_round(k_clip(y_true, 0, 1)))
  true_positives / (possible_positives + k_epsilon())
}

specificity <- function(y_true, y_pred){
  true_negatives = k_sum(k_round(k_clip((1-y_true) * (1-y_pred), 0, 1)))
  possible_negatives = k_sum(k_round(k_clip(1-y_true, 0, 1)))
  true_negatives / (possible_negatives + k_epsilon())
}

precision <- function(y_true, y_pred){
  true_positives = k_sum(k_round(k_clip(y_true * y_pred, 0, 1)))
  predicted_positives = k_sum(k_round(k_clip(y_pred, 0, 1)))
  precision = true_positives / (predicted_positives + k_epsilon())
  precision
}

recall <- function(y_true, y_pred){
  true_positives = k_sum(k_round(k_clip(y_true * y_pred, 0, 1)))
  possible_positives = k_sum(k_round(k_clip(y_true, 0, 1)))
  recall = true_positives / (possible_positives + k_epsilon())
  recall
}


f1 <- function(y_true, y_pred){
  true_positives = k_sum(k_round(k_clip(y_true * y_pred, 0, 1)))
  predicted_positives = k_sum(k_round(k_clip(y_pred, 0, 1)))
  precision = true_positives / (predicted_positives + k_epsilon())
  
  true_positives = k_sum(k_round(k_clip(y_true * y_pred, 0, 1)))
  possible_positives = k_sum(k_round(k_clip(y_true, 0, 1)))
  recall = true_positives / (possible_positives + k_epsilon())
  
  f1 = 2 * ((precision * recall)/(precision + recall))
  f1
}
###############




library(tidyverse)
# keras tensorflow deep learning
library(keras)


set.seed(89345)
train_sub <- model_data$ML_set__eye_dummy_TT$train_set %>% group_by(Status) %>% group_by(Status) %>% sample_n(887) %>% ungroup() %>% select_(.dots=most_imp_predictors_no_disease_class, 'Status')
test_sub <- model_data$ML_set__eye_dummy_TT$test_set %>% select_(.dots=most_imp_predictors_no_disease_class,'Status')
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
  layer_dense(units = 164, activation = 'relu', input_shape=c(39)) %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 164, activation = 'relu') %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 164, activation = 'relu') %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 164, activation = 'relu') %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 164, activation = 'relu') %>% 
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
model_data$ML_set__eye_dummy_TT$test_set$keras <- test_score[,1]
cm_maker('keras', model_data$ML_set__eye_dummy_TT$test_set, cutoff=0.999)

###############
# SMOTE 
###############
library(DMwR)
set.seed(89345)
train_sub <- model_data$ML_set__eye_dummy_TT$train_set %>% select_(.dots=most_imp_predictors_no_disease_class, 'Status')
test_sub <- model_data$ML_set__eye_dummy_TT$test_set %>% select_(.dots=most_imp_predictors_no_disease_class,'Status')

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
  layer_dense(units = 164, activation = 'relu', input_shape=c(39)) %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 164, activation = 'relu') %>% 
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
model_data$ML_set__eye_dummy_TT$test_set$keras <- test_score[,1]
cm_maker('keras', model_data$ML_set__eye_dummy_TT$test_set, cutoff=0.9)


###############
# rnn
###############
set.seed(89345)
train_sub <- model_data$ML_set__eye_dummy_TT$train_set  %>% select_(.dots=most_imp_predictors_no_disease_class, 'Status') 
test_sub <- model_data$ML_set__eye_dummy_TT$test_set %>% select_(.dots=most_imp_predictors_no_disease_class,'Status')

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
dim(train_data) <- c(nrow(train_data),1,39)
dim(test_data) <- c(nrow(test_data),1,39)




model <- keras_model_sequential() %>% 
  layer_lstm(196, recurrent_dropout=0.2, input_shape=c(1, 39), return_sequences = T) %>%
  layer_dropout(0.6) %>%
  layer_lstm(196, recurrent_dropout=0.2, input_shape=c(1, 39), return_sequences = T) %>%
  layer_dropout(0.5) %>%
  layer_lstm(196, recurrent_dropout=0.2, input_shape=c(1, 39), return_sequences = T) %>%
  layer_dropout(0.4) %>%
  layer_lstm(196, recurrent_dropout=0.2, input_shape=c(1, 39), return_sequences = T) %>%
  layer_dropout(0.3) %>%
  layer_lstm(196, recurrent_dropout=0.2, input_shape=c(1, 39), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(196, recurrent_dropout=0.2, input_shape=c(1, 39)) %>%
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
model_data$ML_set__eye_dummy_TT$test_set$keras <- test_score[,1]
cm_maker('keras', model_data$ML_set__eye_dummy_TT$test_set, cutoff=0.99)
