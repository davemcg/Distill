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


allX <- raw_data %>% 
  mutate(DataSet_o=DataSet) %>%
  mutate(DataSet = case_when(pos_id %in% model_data$ML_set__general_TT$train_set$pos_id ~ 'VPaC Train Set',
                             pos_id %in% model_data$ML_set__general_TT$test_set$pos_id ~ 'VPaC Test Set',
                             pos_id %in% model_data$ML_set__other_TT$train_set$pos_id ~ 'VPaC ClinVar LC',
                             pos_id %in% model_data$ML_set__other_TT$test_set$pos_id ~ 'VPaC ClinVar LC',
                             DataSet == 'ddl_nisc_100_panel' ~ 'DDL NISC RD Cohort',
                             pos_id %in% (raw_data %>% filter(DataSet == 'grimm', source == 'humvar') %>% pull(pos_id)) ~ 'Grimm HumVar',
                             pos_id %in% (raw_data %>% filter(DataSet == 'grimm', source == 'exovar') %>% pull(pos_id)) ~ 'Grimm ExoVar',
                             pos_id %in% (raw_data %>% filter(DataSet == 'grimm', source == 'exovar__humvar') %>% pull(pos_id)) ~ 'Grimm ExoVar/HumVar',
                             pos_id %in% (raw_data %>% filter(DataSet == 'grimm', source == 'predictSNP') %>% pull(pos_id)) ~ 'Grimm PredictSNP',
                             pos_id %in% (raw_data %>% filter(DataSet == 'grimm', source == 'swissvar') %>% pull(pos_id)) ~ 'Grimm SwissVar',
                             pos_id %in% (raw_data %>% filter(DataSet == 'grimm', source == 'varibench') %>% pull(pos_id)) ~ 'Grimm VariBench',
                             pos_id %in% (raw_data %>% filter(grepl('homsy', DataSet)) %>% pull(pos_id)) ~ 'Homsy',
                             pos_id %in% (raw_data %>% filter(grepl('unifun', DataSet)) %>% pull(pos_id)) ~ 'UniFun',
                             pos_id %in% (raw_data %>% filter(grepl('samocha', DataSet)) %>% pull(pos_id)) ~ 'Samocha',
                             TRUE ~ DataSet),
         Status = case_when(grepl('pathogenic', DataSet_o, ignore.case = T) | grepl('pathogenic', status, ignore.case = T) ~ 'Pathogenic',
                            (DataSet == 'DDL NISC RD Cohort' & pos_id %in% ddl_path_cdot$pos_id) | 
                              (DataSet == 'DDL NISC RD Cohort' & end %in% ddl_path_cdot$End)  ~ 'Pathogenic',
                            TRUE ~ 'NotPathogenic')) %>% 
  mutate_at(vars(one_of(numeric_predictors)), funs(as.numeric(.))) %>% 
  filter(!grepl('gnomad', DataSet)) 

#rm(raw_data)
allX$fitcons_float <- allX$fitcons
allX[is.na(allX)] <- -1


# calculate VPaC scording for allX
allX$VPaC_m06 <- sqrt(predict(VPaC_6mtry, allX, type='prob')[,1])
allX$VPaC_m15 <- sqrt(predict(VPaC_15mtry, allX, type='prob')[,1])
allX$VPaC_m12 <- sqrt(predict(VPaC_12mtry, allX, type='prob')[,1])
allX$VPaC_m09 <- sqrt(predict(VPaC_9mtry, allX, type='prob')[,1])

save(allX, file='/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/allX_2018_06_20.Rdata')
