
library(tidyverse)
library(randomForest)

load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/model_data_2018_08_01.Rdata')

########################
# Build VPaC
# generic pathogenicity score
#########################


############################
# Only makes 5 trees
# Run 300 times and combine the 300 5-tree forests with `combine`
#############################

############################
# give 32GB for biowulf2 run
############################

######################################
# Predictors
#######################################
most_imp_predictors <- c('is_lof','impact_severity','`DiseaseClass_-1`','mis_z','ccr_pct_v1','cadd_phred','phylop_100way','n_mis','revel','fitcons_float','precessive','n_lof','`DiseaseClass_Stargardt,RD`','m_cap_rankscore','dann','DiseaseClass_RD','vest3_rankscore','n_syn','pnull','pli','lof_z','fathmm_mkl_coding_rankscore','an_exac_all','eigen_pc_raw_rankscore','gerp_elements','mutationassessor_score_rankscore','mpc','metasvm_rankscore','polyphen_score','metalr_rankscore','lrt_converted_rankscore','genocanyon_score_rankscore','mutationtaster_converted_rankscore','gno_an_popmax','grantham','max_aaf_all','ac_exac_all','fathmm_converted_rankscore','aaf_esp_all','sift_score','ac_exac_sas', 'linsight')

numeric_predictors <- c('ccr_pct_v1','cadd_raw','vest3_rankscore','cadd_phred','mis_z','pli','lof_z','phylop_100way','revel','hapmap2','hapmap1','n_mis','epilogos_quies','n_lof','precessive','pnull','adj_exp_lof','adj_exp_syn','dann','adj_exp_mis','syn_z','n_syn','epilogos_txwk','fitcons','m_cap_score','m_cap_rankscore','eigen_phred','eigen_raw','epilogos_tx','is_lof','eigen_pc_raw_rankscore','epilogos_reprpcwk','fathmm_mkl_coding_rankscore','metalr_score','fathmm_mkl_coding_score','metalr_rankscore','impact_severity','metasvm_rankscore','metasvm_score','epilogos_enh','genocanyon_score','fathmm_converted_rankscore','mpc','epilogos_enhg','af_exac_all','epilogos_reprpc','max_aaf_all','mutationassessor_score','gerp','polyphen_score','gerp_elements','mutationassessor_score_rankscore','stam_mean','an_exac_all','af_exac_nfe','provean_converted_rankscore','an_exac_nfe','lrt_score','lrt_omega','grantham','lrt_converted_rankscore','genocanyon_score_rankscore','an_exac_afr','an_exac_amr','an_exac_sas','epilogos_het','ac_exac_all','linsight','gno_an_popmax','exac_num_het','an_exac_eas','gno_an_all','ac_exac_nfe','mutationtaster_converted_rankscore','an_exac_oth','an_exac_fin','gno_an_nfe','gno_af_all','gno_an_afr','epilogos_tssaflnk','gno_af_popmax','epilogos_znf','segway_sum_score','aaf_esp_ea','epilogos_txflnk','provean_score','segway_mean_score','epilogos_tss','aaf_esp_all','af_exac_amr','gno_af_nfe','epilogos_enhbiv','af_exac_sas','sift_score','fathmm_score','ac_exac_amr','aaf_esp_aa','gno_ac_all','gno_af_afr','ac_exac_sas','af_exac_eas','gno_an_fin','af_exac_afr','gno_an_eas','gno_an_oth','gno_ac_nfe','gno_ac_popmax','ac_exac_eas','ac_exac_afr','epilogos_tssbiv','gno_ac_afr','vest3_score','sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01', 'eyeintegration_rnaseq_adipose_subcutaneous','eyeintegration_rnaseq_tpm_adipose_visceral_omentum','eyeintegration_rnaseq_tpm_adrenalgland','eyeintegration_rnaseq_tpm_artery_aorta','eyeintegration_rnaseq_tpm_artery_coronary','eyeintegration_rnaseq_tpm_artery_tibial','eyeintegration_rnaseq_tpm_brain_amygdala','eyeintegration_rnaseq_tpm_brain_anteriorcingulatecortex_ba24','eyeintegration_rnaseq_tpm_brain_caudate_basalganglia','eyeintegration_rnaseq_tpm_brain_cerebellarhemisphere','eyeintegration_rnaseq_tpm_brain_cerebellum','eyeintegration_rnaseq_tpm_brain_cortex','eyeintegration_rnaseq_tpm_brain_frontalcortex_ba9','eyeintegration_rnaseq_tpm_brain_hippocampus','eyeintegration_rnaseq_tpm_brain_hypothalamus','eyeintegration_rnaseq_tpm_brain_nucleusaccumbens_basalganglia','eyeintegration_rnaseq_tpm_brain_putamen_basalganglia','eyeintegration_rnaseq_tpm_brain_spinalcord_cervicalc_1','eyeintegration_rnaseq_tpm_brain_substantianigra','eyeintegration_rnaseq_tpm_breast_mammarytissue','eyeintegration_rnaseq_tpm_cells_ebv_transformedlymphocytes','eyeintegration_rnaseq_tpm_cells_transformedfibroblasts','eyeintegration_rnaseq_tpm_colon_sigmoid','eyeintegration_rnaseq_tpm_colon_transverse','eyeintegration_rnaseq_tpm_esc_stemcellline','eyeintegration_rnaseq_tpm_esophagus_gastroesophagealjunction','eyeintegration_rnaseq_tpm_esophagus_mucosa','eyeintegration_rnaseq_tpm_esophagus_muscularis','eyeintegration_rnaseq_tpm_heart_atrialappendage','eyeintegration_rnaseq_tpm_heart_leftventricle','eyeintegration_rnaseq_tpm_kidney_cortex','eyeintegration_rnaseq_tpm_liver','eyeintegration_rnaseq_tpm_lung','eyeintegration_rnaseq_tpm_minorsalivarygland','eyeintegration_rnaseq_tpm_muscle_skeletal','eyeintegration_rnaseq_tpm_nerve_tibial','eyeintegration_rnaseq_tpm_pancreas','eyeintegration_rnaseq_tpm_pituitary','eyeintegration_rnaseq_tpm_skin_notsunexposed_suprapubic','eyeintegration_rnaseq_tpm_skin_sunexposed_lowerleg','eyeintegration_rnaseq_tpm_smallintestine_terminalileum','eyeintegration_rnaseq_tpm_spleen','eyeintegration_rnaseq_tpm_stomach','eyeintegration_rnaseq_tpm_thyroid','eyeintegration_rnaseq_tpm_wholeblood')

##############################################
# BUILD MODEL!!!!!!!!!!!
#############################################
rf_data <- model_data$ML_set__general_TT$train_set %>% select_(.dots=c('Status',numeric_predictors)) %>% 
  mutate_at(vars(one_of(numeric_predictors)), funs(as.numeric(.))) 
rf_data[is.na(rf_data)] <- 1
rf_data$Status <- factor(rf_data$Status, levels=c('Pathogenic','NotPathogenic'))

rand_num <- as.integer(paste(sample(0:9, 6, replace=F), collapse = ''))
set.seed(rand_num)

rand_num <- as.integer(paste(sample(0:9, 6, replace=F), collapse = ''))
set.seed(rand_num)
rfFit_VPaC <- randomForest(Status ~ ., data=rf_data, 
                           ntree=5,
                           mtry=18,
                           importance = TRUE,
                           norm.votes = FALSE)

rand_name2 = paste0('VPaC__18mtry_v14_', rand_num)
assign(rand_name2, rfFit_VPaC)

rand_num <- as.integer(paste(sample(0:9, 6, replace=F), collapse = ''))
set.seed(rand_num)

rfFit_VPaC <- randomForest(Status ~ ., data=rf_data, 
                           ntree=5,
                           mtry=15,
                           importance = TRUE,
                           norm.votes = FALSE)

rand_name3 = paste0('VPaC__15mtry_v14_', rand_num)
assign(rand_name3, rfFit_VPaC)

rfFit_VPaC <- randomForest(Status ~ ., data=rf_data, 
                           ntree=5,
                           mtry=12,
                           importance = TRUE,
                           norm.votes = FALSE)

rand_name4 = paste0('VPaC__12mtry_v14_', rand_num)
assign(rand_name4, rfFit_VPaC)



##############################
# SAVE MODEL
###############################
#save(list=rand_name1, file=paste0('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_pieces/', rand_name1, '.Rdata'))
save(list=rand_name2, file=paste0('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_pieces/', rand_name2, '.Rdata'))
save(list=rand_name3, file=paste0('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_pieces/', rand_name3, '.Rdata'))
save(list=rand_name4, file=paste0('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_pieces/', rand_name4, '.Rdata'))
