
library(tidyverse)
library(randomForest)

load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/model_data.Rdata')

########################
# Build VPaC
# generic pathogenicity score
#########################


############################
# Only makes 10 trees
# Run 150 times and combine the 150 10-tree forests with `combine`
#############################

############################
# give 16GB for biowulf2 run
############################

######################################
# Predictors
#######################################
most_imp_predictors <- c('is_lof','impact_severity','`DiseaseClass_-1`','mis_z','ccr_pct_v1','cadd_phred','phylop_100way','n_mis','revel','fitcons_float','precessive','n_lof','`DiseaseClass_Stargardt,RD`','m_cap_rankscore','dann','DiseaseClass_RD','vest3_rankscore','n_syn','pnull','pli','lof_z','fathmm_mkl_coding_rankscore','an_exac_all','eigen_pc_raw_rankscore','gerp_elements','mutationassessor_score_rankscore','mpc','metasvm_rankscore','polyphen_score','metalr_rankscore','lrt_converted_rankscore','genocanyon_score_rankscore','mutationtaster_converted_rankscore','gno_an_popmax','grantham','max_aaf_all','ac_exac_all','fathmm_converted_rankscore','aaf_esp_all','sift_score','ac_exac_sas', 'linsight')
most_imp_predictors_no_disease_class <- c('is_lof','impact_severity','mis_z','ccr_pct_v1','cadd_phred','phylop_100way','n_mis','revel','fitcons_float','precessive','n_lof','m_cap_rankscore','dann','vest3_rankscore','n_syn','pnull','pli','lof_z','fathmm_mkl_coding_rankscore','an_exac_all','eigen_pc_raw_rankscore','gerp_elements','mutationassessor_score_rankscore','mpc','metasvm_rankscore','polyphen_score','metalr_rankscore','lrt_converted_rankscore','genocanyon_score_rankscore','mutationtaster_converted_rankscore','gno_an_popmax','grantham','max_aaf_all','ac_exac_all','fathmm_converted_rankscore','aaf_esp_all','sift_score','ac_exac_sas','linsight')


##############################################
# BUILD MODEL!!!!!!!!!!!
#############################################

rand_num <- as.integer(paste(sample(0:9, 6, replace=F), collapse = ''))
set.seed(rand_num)
model_data$ML_set__general_dummy_TT$train_set$Status <- factor(model_data$ML_set__general_dummy_TT$train_set$Status, levels=c('Pathogenic','NotPathogenic'))
rfFit_VPaC <- randomForest(Status ~ ., data=model_data$ML_set__general_dummy_TT$train_set %>% select_(.dots=c('Status',most_imp_predictors_no_disease_class)), 
                           ntree=33,
                           mtry=6,
                           importance = TRUE,
                           norm.votes = FALSE)

rand_name = paste0('VPaC__6mtry_v2', rand_num)
assign(rand_name, rfFit_VPaC)


##############################
# SAVE MODEL
###############################
save(list=rand_name, file=paste0('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_pieces/', rand_name, '.Rdata'))
