load('/Volumes/data/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/model_spread.Rdata')
most_imp_predictors_no_disease_class <- c('is_lof','impact_severity','mis_z','ccr_pct_v1','cadd_phred','phylop_100way','n_mis','revel','fitcons','precessive','n_lof','m_cap_rankscore','dann','vest3_rankscore','n_syn','pnull','pli','lof_z','fathmm_mkl_coding_rankscore','an_exac_all','eigen_pc_raw_rankscore','gerp_elements','mutationassessor_score_rankscore','mpc','metasvm_rankscore','polyphen_score','metalr_rankscore','lrt_converted_rankscore','genocanyon_score_rankscore','mutationtaster_converted_rankscore','gno_an_popmax','grantham','max_aaf_all','ac_exac_all','fathmm_converted_rankscore','aaf_esp_all','sift_score','ac_exac_sas','linsight','segway_mean_score','segway_sum_score','epilogos_bivflnk','epilogos_enh','epilogos_enhbiv','epilogos_enhg','epilogos_het','epilogos_quies','epilogos_reprpc','epilogos_reprpcwk','epilogos_tss','epilogos_tssaflnk','epilogos_tssbiv','epilogos_tx','epilogos_txflnk','epilogos_txwk','epilogos_znf')

most_imp_predictors_no_disease_class <- c('is_lof','impact_severity','mis_z','ccr_pct_v1','cadd_phred','phylop_100way','n_mis','revel','fitcons','precessive','n_lof','m_cap_rankscore','dann','vest3_rankscore','n_syn','pnull','pli','lof_z','fathmm_mkl_coding_rankscore','an_exac_all','eigen_pc_raw_rankscore','gerp_elements','mutationassessor_score_rankscore','mpc','metasvm_rankscore','polyphen_score','metalr_rankscore','lrt_converted_rankscore','genocanyon_score_rankscore','mutationtaster_converted_rankscore','gno_an_popmax','grantham','max_aaf_all','ac_exac_all','fathmm_converted_rankscore','aaf_esp_all','sift_score','ac_exac_sas','linsight')

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
    out <- caret::confusionMatrix(data = as.factor(new_predictions$Prediction), reference = as.factor(new_predictions$Status), mode= mode)
    out$MCC <- mcc(new_predictions$preds, new_predictions$actuals, cutoff=cutoff)
  }
  out
}

################
# build 3d array
# matching variant with surrounding spread data
#################

clinvar_spread$train_set$chr = sapply(clinvar_spread$train_set$pos_id, function(x) str_split(x, ':')[[1]][1])
clinvar_spread$train_set$pos = sapply(clinvar_spread$train_set$pos_id, function(x) (str_split(x, ':')[[1]][2] %>% str_split(. ,'_'))[[1]][1] %>% as.numeric() )
clinvar_spread$train_set$ref = sapply(clinvar_spread$train_set$pos_id, function(x) (str_split(x, ':')[[1]][2] %>% str_split(. ,'_'))[[1]][2] )
clinvar_spread$train_set$alt = sapply(clinvar_spread$train_set$pos_id, function(x) (str_split(x, ':')[[1]][2] %>% str_split(. ,'_'))[[1]][3] )


clinvar_spread$test_set$chr = sapply(clinvar_spread$test_set$pos_id, function(x) str_split(x, ':')[[1]][1])
clinvar_spread$test_set$pos = sapply(clinvar_spread$test_set$pos_id, function(x) (str_split(x, ':')[[1]][2] %>% str_split(. ,'_'))[[1]][1] %>% as.numeric())
clinvar_spread$test_set$ref = sapply(clinvar_spread$test_set$pos_id, function(x) (str_split(x, ':')[[1]][2] %>% str_split(. ,'_'))[[1]][2] )
clinvar_spread$test_set$alt = sapply(clinvar_spread$test_set$pos_id, function(x) (str_split(x, ':')[[1]][2] %>% str_split(. ,'_'))[[1]][3] )

clinvar_spread$spread$chr = sapply(clinvar_spread$spread$pos_id, function(x) str_split(x, ':')[[1]][1])
clinvar_spread$spread$pos = sapply(clinvar_spread$spread$pos_id, function(x) (str_split(x, ':')[[1]][2] %>% str_split(. ,'_'))[[1]][1] %>% as.numeric() )
clinvar_spread$spread$ref = sapply(clinvar_spread$spread$pos_id, function(x) (str_split(x, ':')[[1]][2] %>% str_split(. ,'_'))[[1]][2] )
clinvar_spread$spread$alt = sapply(clinvar_spread$spread$pos_id, function(x) (str_split(x, ':')[[1]][2] %>% str_split(. ,'_'))[[1]][3] )

# calc mean and std for scaling
train_sub <-clinvar_spread$train_set %>% select_(.dots=most_imp_predictors_no_disease_class, 'Status')
mean <- apply(train_sub %>% select_(.dots = most_imp_predictors_no_disease_class), 2, mean)
std <- apply(train_sub %>% select_(.dots = most_imp_predictors_no_disease_class), 2, sd)


test_array <- array(NA, c(nrow(clinvar_spread$test_set), 13, length(most_imp_predictors_no_disease_class)))
for (i in seq(1:nrow(clinvar_spread$test_set))){
  end <- clinvar_spread$test_set[i, 'pos']
  chrom <-  clinvar_spread$test_set[i, 'chr']
  real <-  clinvar_spread$test_set[i, most_imp_predictors_no_disease_class]
  spread_set <- clinvar_spread$spread %>% 
    filter(chr== as.character(chrom) & pos >= as.numeric(end) - 2 & pos<= as.numeric(end) +2) %>% 
    filter(pos != as.numeric(end)) %>% 
    arrange(pos, alt) %>% select_(.dots = most_imp_predictors_no_disease_class)
  array_spread <- as.matrix(rbind(real, spread_set))

  array_spread <- array(array_spread, c(13, length(most_imp_predictors_no_disease_class)))
  
  #print(dim(array_spread))
  #print(i)
  array_spread <- scale(array_spread, center=mean, scale=std)
  
  test_array[i,,]  <- as.matrix(array_spread)
}

train_array <- array(NA, c(nrow(clinvar_spread$train_set), 13, length(most_imp_predictors_no_disease_class)))
for (i in seq(1:nrow(clinvar_spread$train_set))){
  end <- clinvar_spread$train_set[i, 'pos']
  chrom <-  clinvar_spread$train_set[i, 'chr']
  real <-  clinvar_spread$train_set[i, most_imp_predictors_no_disease_class]
  spread_set <- clinvar_spread$spread %>% 
    filter(chr== as.character(chrom) & pos >= as.numeric(end) - 2 & pos<= as.numeric(end) +2) %>% 
    filter(pos != as.numeric(end)) %>% 
    arrange(pos, alt) %>% select_(.dots = most_imp_predictors_no_disease_class)
  array_spread <- as.matrix(rbind(real, spread_set))
  array_spread <- array(array_spread, c(13, length(most_imp_predictors_no_disease_class)))
  
  #print(dim(array_spread))
  #print(i)
  array_spread <- scale(array_spread, center=mean, scale=std)
  train_array[i,,]  <- as.matrix(array_spread)
}



###############
# rnn smote 1
# adds new minor class (path) and trim down major class
# currently this is the model that works the best when ensembl-ing with RF-based VPaC
###############
set.seed(89345)
train_sub <-clinvar_spread$train_set %>% select_(.dots=most_imp_predictors_no_disease_class, 'Status')
test_sub <-clinvar_spread$test_set %>% select_(.dots=most_imp_predictors_no_disease_class,'Status')

#train_sub <- SMOTE(Status ~ ., as.data.frame(train_sub))
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
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(1, length(most_imp_predictors_no_disease_class))) %>%
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
clinvar_spread$test_set$keras <- test_score[,1]
cm_maker('keras',clinvar_spread$test_set, cutoff=0.9)

##########################################
# multidimensional
##########################################

set.seed(89345)
#train_array 
status_train <- clinvar_spread$train_set$Status
status_train01 <- case_when(status_train == 'Pathogenic' ~ 1,
                            TRUE ~ 0)
status_test <- clinvar_spread$test_set$Status
status_test01 <- case_when(status_test == 'Pathogenic' ~ 1,
                           TRUE ~ 0)


model <- keras_model_sequential() %>% 
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(13, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(13, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(13, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(13, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(13, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(13, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(13, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(13, length(most_imp_predictors_no_disease_class))) %>%
  layer_dropout(0.1) %>%
  layer_dense(units=1, activation='sigmoid')



model %>% compile(
  optimizer = optimizer_adam(),
  loss = 'binary_crossentropy',
  metric = c('accuracy')
)

history <- model %>% fit(train_array, status_train01, epochs = 10, batch_size=50)

#validation_score <- model %>% evaluate(test_data, status_test01)
test_score <- model %>% predict(test_array)
test_score[is.na(test_score)] <- 0
clinvar_spread$test_set$kerasSPREAD <- test_score[,1]
cm_maker('kerasSPREAD',clinvar_spread$test_set, cutoff=0.6)


#######################
# multidimensional
# cnn
########################





model <- keras_model_sequential() %>% 
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(13, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(13, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(13, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(13, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(13, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(13, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(13, length(most_imp_predictors_no_disease_class)), return_sequences = T) %>%
  layer_dropout(0.2) %>%
  layer_lstm(48, recurrent_dropout=0.2, input_shape=c(13, length(most_imp_predictors_no_disease_class))) %>%
  layer_dropout(0.1) %>%
  layer_dense(units=1, activation='sigmoid')

  
  # Start with hidden 2D convolutional layer being fed 32x32 pixel images
  layer_conv_2d(
    filter = 32, kernel_size = c(3,3), padding = "same", 
    input_shape = c(32, 32, 3)
  ) %>%
  layer_activation("relu") %>%
  
  # Second hidden layer
  layer_conv_2d(filter = 32, kernel_size = c(3,3)) %>%
  layer_activation("relu") %>%
  
  # Use max pooling
  layer_max_pooling_2d(pool_size = c(2,2)) %>%
  layer_dropout(0.25) %>%
  
  
  
