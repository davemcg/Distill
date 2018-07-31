
library(tidyverse)
library(randomForest)

load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/model_data_2018_07_31.Rdata')

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

numeric_predictors <- c('is_exonic','is_coding','is_lof','is_splicing','impact_severity','polyphen_score','sift_score','dann','eigen_phred','fathmm_converted_rankscore','gerp','genocanyon_score_rankscore','hgmd_overlap','linsight','lrt_omega','lrt_converted_rankscore','m_cap_rankscore','mpc','metalr_rankscore','metasvm_rankscore','mutationassessor_score_rankscore','mutationtaster_converted_rankscore','provean_converted_rankscore','provean_score','revel','vest3_rankscore','aaf_1kg_afr','aaf_1kg_all','aaf_1kg_amr','aaf_1kg_eas','aaf_1kg_eur','aaf_1kg_sas','aaf_esp_aa','aaf_esp_all','aaf_esp_ea','ac_exac_afr','ac_exac_all','ac_exac_amr','ac_exac_eas','ac_exac_fin','ac_exac_nfe','ac_exac_oth','ac_exac_sas','adj_exp_lof','adj_exp_mis','adj_exp_syn','af_exac_afr','af_exac_all','af_exac_amr','af_exac_eas','af_exac_nfe','af_exac_oth','af_exac_sas','an_exac_afr','an_exac_all','an_exac_amr','an_exac_eas','an_exac_fin','an_exac_nfe','an_exac_oth','an_exac_sas','ccr_pct_v1','cpg_island','epilogos_bivflnk','epilogos_enh','epilogos_enhbiv','epilogos_enhg','epilogos_het','epilogos_quies','epilogos_reprpc','epilogos_reprpcwk','epilogos_tss','epilogos_tssaflnk','epilogos_tssbiv','epilogos_tx','epilogos_txflnk','epilogos_txwk','epilogos_znf','exac_num_het','exac_num_hom_alt','fathmm_mkl_coding_rankscore','fitcons','geno2mp','gerp_elements','gno_ac_afr','gno_ac_all','gno_ac_amr','gno_ac_asj','gno_ac_eas','gno_ac_fin','gno_ac_nfe','gno_ac_oth','gno_ac_popmax','gno_af_afr','gno_af_all','gno_af_amr','gno_af_asj','gno_af_eas','gno_af_fin','gno_af_nfe','gno_af_oth','gno_af_popmax','gno_an_afr','gno_an_all','gno_an_amr','gno_an_asj','gno_an_eas','gno_an_fin','gno_an_nfe','gno_an_oth','gno_an_popmax','gno_id','gno_popmax','hapmap1','hapmap2','in_1kg','in_esp','in_exac','lof_z','max_aaf_all','mis_z','n_lof','n_mis','n_syn','pli','pnull','precessive','phylop_100way','segway_mean_score','segway_sum_score','stam_mean','syn_z','grantham','cadd_phred', 'sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01','eyeintegration_rnaseq_adipose_subcutaneous','eyeintegration_rnaseq_tpm_adipose_visceral_omentum','eyeintegration_rnaseq_tpm_adrenalgland','eyeintegration_rnaseq_tpm_artery_aorta','eyeintegration_rnaseq_tpm_artery_coronary','eyeintegration_rnaseq_tpm_artery_tibial','eyeintegration_rnaseq_tpm_brain_amygdala','eyeintegration_rnaseq_tpm_brain_anteriorcingulatecortex_ba24','eyeintegration_rnaseq_tpm_brain_caudate_basalganglia','eyeintegration_rnaseq_tpm_brain_cerebellarhemisphere','eyeintegration_rnaseq_tpm_brain_cerebellum','eyeintegration_rnaseq_tpm_brain_cortex','eyeintegration_rnaseq_tpm_brain_frontalcortex_ba9','eyeintegration_rnaseq_tpm_brain_hippocampus','eyeintegration_rnaseq_tpm_brain_hypothalamus','eyeintegration_rnaseq_tpm_brain_nucleusaccumbens_basalganglia','eyeintegration_rnaseq_tpm_brain_putamen_basalganglia','eyeintegration_rnaseq_tpm_brain_spinalcord_cervicalc_1','eyeintegration_rnaseq_tpm_brain_substantianigra','eyeintegration_rnaseq_tpm_breast_mammarytissue','eyeintegration_rnaseq_tpm_cells_ebv_transformedlymphocytes','eyeintegration_rnaseq_tpm_cells_transformedfibroblasts','eyeintegration_rnaseq_tpm_colon_sigmoid','eyeintegration_rnaseq_tpm_colon_transverse','eyeintegration_rnaseq_tpm_cornea_adulttissue','eyeintegration_rnaseq_tpm_cornea_cellline','eyeintegration_rnaseq_tpm_cornea_fetaltissue','eyeintegration_rnaseq_tpm_esc_stemcellline','eyeintegration_rnaseq_tpm_esophagus_gastroesophagealjunction','eyeintegration_rnaseq_tpm_esophagus_mucosa','eyeintegration_rnaseq_tpm_esophagus_muscularis','eyeintegration_rnaseq_tpm_heart_atrialappendage','eyeintegration_rnaseq_tpm_heart_leftventricle','eyeintegration_rnaseq_tpm_kidney_cortex','eyeintegration_rnaseq_tpm_liver','eyeintegration_rnaseq_tpm_lung','eyeintegration_rnaseq_tpm_minorsalivarygland','eyeintegration_rnaseq_tpm_muscle_skeletal','eyeintegration_rnaseq_tpm_nerve_tibial','eyeintegration_rnaseq_tpm_pancreas','eyeintegration_rnaseq_tpm_pituitary','eyeintegration_rnaseq_tpm_rpe_adulttissue','eyeintegration_rnaseq_tpm_rpe_cellline','eyeintegration_rnaseq_tpm_rpe_fetaltissue','eyeintegration_rnaseq_tpm_rpe_stemcellline','eyeintegration_rnaseq_tpm_retina_adulttissue','eyeintegration_rnaseq_tpm_retina_stemcellline','eyeintegration_rnaseq_tpm_skin_notsunexposed_suprapubic','eyeintegration_rnaseq_tpm_skin_sunexposed_lowerleg','eyeintegration_rnaseq_tpm_smallintestine_terminalileum','eyeintegration_rnaseq_tpm_spleen','eyeintegration_rnaseq_tpm_stomach','eyeintegration_rnaseq_tpm_thyroid','eyeintegration_rnaseq_tpm_wholeblood')

##############################################
# BUILD MODEL!!!!!!!!!!!
#############################################
rf_data <- model_data$ML_set__general_TT$train_set %>% select_(.dots=c('Status',numeric_predictors))
rf_data[is.na(rf_data)] <- 1
rf_data$Status <- factor(rf_data$Status, levels=c('Pathogenic','NotPathogenic'))

rand_num <- as.integer(paste(sample(0:9, 6, replace=F), collapse = ''))
set.seed(rand_num)

rand_num <- as.integer(paste(sample(0:9, 6, replace=F), collapse = ''))
set.seed(rand_num)
rfFit_VPaC <- randomForest(Status ~ ., data=rf_data, 
                           ntree=10,
                           mtry=18,
                           importance = TRUE,
                           norm.votes = FALSE)

rand_name2 = paste0('VPaC__18mtry_v12_', rand_num)
assign(rand_name2, rfFit_VPaC)

rand_num <- as.integer(paste(sample(0:9, 6, replace=F), collapse = ''))
set.seed(rand_num)

rfFit_VPaC <- randomForest(Status ~ ., data=rf_data, 
                           ntree=15,
                           mtry=12,
                           importance = TRUE,
                           norm.votes = FALSE)

rand_name3 = paste0('VPaC__15mtry_v12_', rand_num)
assign(rand_name3, rfFit_VPaC)

rfFit_VPaC <- randomForest(Status ~ ., data=rf_data, 
                           ntree=12,
                           mtry=9,
                           importance = TRUE,
                           norm.votes = FALSE)

rand_name4 = paste0('VPaC__12mtry_v12_', rand_num)
assign(rand_name4, rfFit_VPaC)



##############################
# SAVE MODEL
###############################
#save(list=rand_name1, file=paste0('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_pieces/', rand_name1, '.Rdata'))
save(list=rand_name2, file=paste0('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_pieces/', rand_name2, '.Rdata'))
save(list=rand_name3, file=paste0('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_pieces/', rand_name3, '.Rdata'))
save(list=rand_name4, file=paste0('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_pieces/', rand_name4, '.Rdata'))
