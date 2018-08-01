################################################
# Process the gemini output data for the modeling
################################################

# biowulf paths
## run: Rscript ~/git/eye_var_Pathogenicity/scripts/build_UK10K_data.R
uk10k_data <- '/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/uk10k_gemini_rare_variants_2018_07_31.Rdata'
## run: time gemini query --header -q "SELECT * from variants WHERE (aaf_esp_all < 0.01 AND aaf_1kg_all < 0.01 AND af_exac_all < 0.01  AND (is_coding=1 OR is_splicing=1)) OR impact_severity='HIGH' OR clinvar_sig LIKE '%patho%'" clinvar_RD.PED_faux.gemini.db | bgzip > clinvar.gemini.tsv.gz 
clinvar_file <- '/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/clinvar/clinvar.gemini.tsv.gz'
## run: time gemini query --header -q "SELECT * from variants WHERE (aaf_esp_all < 0.01 AND aaf_1kg_all < 0.01 AND af_exac_all < 0.01  AND (is_coding=1 OR is_splicing=1)) OR impact_severity='HIGH' OR clinvar_sig LIKE '%patho%'" gnomad.exomes.r2.0.2.sites.maxAF_0.01_20percent.PED_faux.gemini.db | bgzip > gnomad_rare_benign_ish.gemini.tsv.gz
gnomad_file <- '/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/gnomad_rare_benign_ish/gnomad_rare_benign_ish.gemini.tsv.gz'
# contains variants adjacent to existing variants in clinvar
#spread_file <- '/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/clinvar/spread/clinvar.gemini.tsv.gz'

library(tidyverse)
library(data.table)
library(dummies)
library(caret)

###############################
# Highly correlated predictors removed by inspecting correlation heatmap in remove_correlated_predictors.Rmd
# 
# or rather unique(ish) predictors kept
# reflected in the `select` commands below
###############################
info <- c('pos_id', 
          'Status', 
          'Complicated_Status')

cat_predictors <- c('genesplicer',
                    'spliceregion')

numeric_predictors<-c('is_exonic','is_coding','is_lof','is_splicing','impact_severity','polyphen_score','sift_score','dann','eigen_phred','fathmm_converted_rankscore','gerp','genocanyon_score_rankscore','hgmd_overlap','linsight','lrt_omega','lrt_converted_rankscore','m_cap_rankscore','mpc','metalr_rankscore','metasvm_rankscore','mutationassessor_score_rankscore','mutationtaster_converted_rankscore','provean_converted_rankscore','provean_score','revel','vest3_rankscore','aaf_1kg_afr','aaf_1kg_all','aaf_1kg_amr','aaf_1kg_eas','aaf_1kg_eur','aaf_1kg_sas','aaf_esp_aa','aaf_esp_all','aaf_esp_ea','ac_exac_afr','ac_exac_all','ac_exac_amr','ac_exac_eas','ac_exac_fin','ac_exac_nfe','ac_exac_oth','ac_exac_sas','adj_exp_lof','adj_exp_mis','adj_exp_syn','af_exac_afr','af_exac_all','af_exac_amr','af_exac_eas','af_exac_nfe','af_exac_oth','af_exac_sas','an_exac_afr','an_exac_all','an_exac_amr','an_exac_eas','an_exac_fin','an_exac_nfe','an_exac_oth','an_exac_sas','ccr_pct_v1','cpg_island','epilogos_bivflnk','epilogos_enh','epilogos_enhbiv','epilogos_enhg','epilogos_het','epilogos_quies','epilogos_reprpc','epilogos_reprpcwk','epilogos_tss','epilogos_tssaflnk','epilogos_tssbiv','epilogos_tx','epilogos_txflnk','epilogos_txwk','epilogos_znf','exac_num_het','exac_num_hom_alt','fathmm_mkl_coding_rankscore','fitcons','geno2mp','gerp_elements','gno_ac_afr','gno_ac_all','gno_ac_amr','gno_ac_asj','gno_ac_eas','gno_ac_fin','gno_ac_nfe','gno_ac_oth','gno_ac_popmax','gno_af_afr','gno_af_all','gno_af_amr','gno_af_asj','gno_af_eas','gno_af_fin','gno_af_nfe','gno_af_oth','gno_af_popmax','gno_an_afr','gno_an_all','gno_an_amr','gno_an_asj','gno_an_eas','gno_an_fin','gno_an_nfe','gno_an_oth','gno_an_popmax','gno_id','gno_popmax','hapmap1','hapmap2','in_1kg','in_esp','in_exac','lof_z','max_aaf_all','mis_z','n_lof','n_mis','n_syn','pli','pnull','precessive','phylop_100way','segway_mean_score','segway_sum_score','stam_mean','syn_z','grantham','cadd_phred', 'sigmaaf_lof_0001', 'sigmaaf_lof_01', 'sigmaaf_missense_0001', 'sigmaaf_missense_01', 'primatedl','eyeintegration_rnaseq_adipose_subcutaneous','eyeintegration_rnaseq_tpm_adipose_visceral_omentum','eyeintegration_rnaseq_tpm_adrenalgland','eyeintegration_rnaseq_tpm_artery_aorta','eyeintegration_rnaseq_tpm_artery_coronary','eyeintegration_rnaseq_tpm_artery_tibial','eyeintegration_rnaseq_tpm_brain_amygdala','eyeintegration_rnaseq_tpm_brain_anteriorcingulatecortex_ba24','eyeintegration_rnaseq_tpm_brain_caudate_basalganglia','eyeintegration_rnaseq_tpm_brain_cerebellarhemisphere','eyeintegration_rnaseq_tpm_brain_cerebellum','eyeintegration_rnaseq_tpm_brain_cortex','eyeintegration_rnaseq_tpm_brain_frontalcortex_ba9','eyeintegration_rnaseq_tpm_brain_hippocampus','eyeintegration_rnaseq_tpm_brain_hypothalamus','eyeintegration_rnaseq_tpm_brain_nucleusaccumbens_basalganglia','eyeintegration_rnaseq_tpm_brain_putamen_basalganglia','eyeintegration_rnaseq_tpm_brain_spinalcord_cervicalc_1','eyeintegration_rnaseq_tpm_brain_substantianigra','eyeintegration_rnaseq_tpm_breast_mammarytissue','eyeintegration_rnaseq_tpm_cells_ebv_transformedlymphocytes','eyeintegration_rnaseq_tpm_cells_transformedfibroblasts','eyeintegration_rnaseq_tpm_colon_sigmoid','eyeintegration_rnaseq_tpm_colon_transverse','eyeintegration_rnaseq_tpm_cornea_adulttissue','eyeintegration_rnaseq_tpm_cornea_cellline','eyeintegration_rnaseq_tpm_cornea_fetaltissue','eyeintegration_rnaseq_tpm_esc_stemcellline','eyeintegration_rnaseq_tpm_esophagus_gastroesophagealjunction','eyeintegration_rnaseq_tpm_esophagus_mucosa','eyeintegration_rnaseq_tpm_esophagus_muscularis','eyeintegration_rnaseq_tpm_heart_atrialappendage','eyeintegration_rnaseq_tpm_heart_leftventricle','eyeintegration_rnaseq_tpm_kidney_cortex','eyeintegration_rnaseq_tpm_liver','eyeintegration_rnaseq_tpm_lung','eyeintegration_rnaseq_tpm_minorsalivarygland','eyeintegration_rnaseq_tpm_muscle_skeletal','eyeintegration_rnaseq_tpm_nerve_tibial','eyeintegration_rnaseq_tpm_pancreas','eyeintegration_rnaseq_tpm_pituitary','eyeintegration_rnaseq_tpm_rpe_adulttissue','eyeintegration_rnaseq_tpm_rpe_cellline','eyeintegration_rnaseq_tpm_rpe_fetaltissue','eyeintegration_rnaseq_tpm_rpe_stemcellline','eyeintegration_rnaseq_tpm_retina_adulttissue','eyeintegration_rnaseq_tpm_retina_stemcellline','eyeintegration_rnaseq_tpm_skin_notsunexposed_suprapubic','eyeintegration_rnaseq_tpm_skin_sunexposed_lowerleg','eyeintegration_rnaseq_tpm_smallintestine_terminalileum','eyeintegration_rnaseq_tpm_spleen','eyeintegration_rnaseq_tpm_stomach','eyeintegration_rnaseq_tpm_thyroid','eyeintegration_rnaseq_tpm_wholeblood')

###############################
# UK10K processing
###############################
load(uk10k_data)
## Set up data for modeling
all_processed <- uk10k_gemini_rare_variants %>% 
  separate(gene_eyediseaseclass, c('RDGene','DiseaseClass'), sep='_') %>%  #split off RD disease type
  select(-RDGene) %>% 
  mutate(pos_id=paste0(chrom, ':', end, '_', ref, '_', alt),
         impact_severity = case_when(impact_severity == 'HIGH' ~ 3, # convert to integer 
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
  mutate_at(vars(one_of(numeric_predictors)), funs(as.numeric(.))) %>%  # convert columns with ac_|whatever to integer (ac is allele count), etc. af is allele frequency
  #select(one_of(info), one_of(numeric_predictors), one_of(cat_predictors)) %>% 
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
#all_set__uk10k$Source <- 'UK10K'
pos_id__source <- all_set__uk10k %>% select(pos_id, Status) %>% mutate(Source='UK10K')

all_PATH <- all_set__uk10k %>% filter(Status=='Pathogenic')
# chance to cut down non pathogenic
# i'm just keeping all right now
set.seed(115470)
all_NOT_PATH__CUT <- all_set__uk10k %>% filter(Status=='NotPathogenic') #%>% sample_n(20000)


#####
# split train and test
#####

ML_set__UK10K <- rbind(all_PATH %>% filter(DataSet == 'UK10 Train'), all_NOT_PATH__CUT %>% filter(DataSet == 'UK10 Train'))
Test_set__UK10K <- rbind(all_PATH %>% filter(DataSet == 'UK10 Test'), all_NOT_PATH__CUT %>% filter(DataSet == 'UK10 Test'))





print('UK10K Loaded')
###########################################
# ClinVar Processing
###########################################
#clinvar <- fread('gzcat ~/git/eye_var_Pathogenicity/processed_data/clinvar.gemini.tsv.gz')
clinvar <- fread(paste0('gzcat ', clinvar_file))

## Prep data for modeling
clinvar_processed <- clinvar %>% 
  #filter(status!='PATHOGENIC_OTHER') %>% # drop non eye pathogenic variants for the model learning 
  separate(gene_eyediseaseclass, c('RDGene','DiseaseClass'), sep='_') %>%  #split off RD disease type
  select(-RDGene) %>% 
  mutate(pos_id=paste0(chrom, ':', end, '_', ref, '_', alt),
         impact_severity = case_when(impact_severity == 'HIGH' ~ 3, # convert to integer 
                                     impact_severity == 'MED' ~ 2, 
                                     TRUE ~ 1),
         Status = case_when(status=='PATHOGENIC_EYE' ~ 'Pathogenic',
                            status=='PATHOGENIC_OTHER' ~ 'Pathogenic_NOTEYE',
                            status=='PATHOGENIC_OTHER_HC' ~ 'Pathogenic_NOTEYE_HC',
                            TRUE ~ 'NotPathogenic'),
         genesplicer = case_when(genesplicer == "" ~ 'No',
                                 grepl('^gain', genesplicer) ~ 'Gain',
                                 grepl('^loss', genesplicer) ~ 'Loss',
                                 grepl('^diff', genesplicer) ~ 'Diff',
                                 TRUE ~ 'Else')) %>% 
  mutate_at(vars(one_of(numeric_predictors)), funs(as.numeric(.))) %>%  # convert columns with ac_|whatever to integer (ac is allele count), etc. af is allele frequency
  #select(one_of(info), one_of(numeric_predictors), one_of(cat_predictors)) %>% 
  filter(max_aaf_all < 0.01) # remove any common variants

# fill missing with -1
clinvar_processed[is.na(clinvar_processed)] <- -1

ML_set__clinvar <- clinvar_processed %>% 
  filter(Status != 'Pathogenic_NOTEYE' & Status != 'Pathogenic_NOTEYE_HC') %>% 
  mutate(Status = factor(Status, levels=c('Pathogenic','NotPathogenic'))) 
#ML_set__clinvar$Source <- 'ClinVar'
pos_id__source <- rbind(pos_id__source, 
                        ML_set__clinvar %>% select(pos_id, Status) %>% mutate(Source='ClinVar'))

ML_set__clinvar__otherPath_HC <- clinvar_processed %>% 
  filter(Status == 'Pathogenic_NOTEYE_HC') %>% 
  mutate(Status = gsub('Pathogenic_NOTEYE_HC','Pathogenic',Status)) %>% 
  mutate(Status = factor(Status, levels=c('Pathogenic','NotPathogenic'))) 

pos_id__source <- rbind(pos_id__source, 
                        ML_set__clinvar__otherPath_HC %>% select(pos_id, Status) %>% mutate(Source='ClinVar_OtherPath_HC'))


ML_set__clinvar__otherPath <- clinvar_processed %>% 
  filter(Status == 'Pathogenic_NOTEYE') %>% 
  mutate(Status = gsub('Pathogenic_NOTEYE','Pathogenic',Status)) %>% 
  mutate(Status = factor(Status, levels=c('Pathogenic','NotPathogenic'))) 
#ML_set__clinvar__otherPath$Source <- 'ClinVar'
pos_id__source <- rbind(pos_id__source, 
                        ML_set__clinvar__otherPath %>% select(pos_id, Status) %>% mutate(Source='ClinVar_OtherPath'))
print('ClinVar Loaded')
###############################################
# gnomAD benign? processing
##############################################
gnomad <- fread(paste0('gzcat ', gnomad_file))
## Prep data for modeling
gnomad_processed <- gnomad %>% 
  #filter(gene %in% (clinvar %>% filter(status=='PATHOGENIC_EYE' | status == 'PATHOGENIC_OTHER_HC') %>% pull(gene) %>% unique())) %>% # only have gnomad matched for path variants in clinvar genes
  separate(gene_eyediseaseclass, c('RDGene','DiseaseClass'), sep='_') %>%  #split off RD disease type
  select(-RDGene) %>% 
  mutate(pos_id=paste0(chrom, ':', end, '_', ref, '_', alt),
         impact_severity = case_when(impact_severity == 'HIGH' ~ 3, # convert to integer 
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
  mutate_at(vars(one_of(numeric_predictors)), funs(as.numeric(.))) %>%  # convert columns with ac_|whatever to integer (ac is allele count), etc. af is allele frequency
  #select(one_of(info), one_of(numeric_predictors), one_of(cat_predictors)) %>% 
  filter(max_aaf_all < 0.01) # remove any common variants 


# fill missing with -1
gnomad_processed[is.na(gnomad_processed)] <- -1

pos_id__source <- rbind(pos_id__source, 
                        gnomad_processed %>% select(pos_id, Status) %>% mutate(Source='gnomAD')) %>% 
  unique()

# no more than 250X gnomad benign relative to the number of clinvar pathogenic variants 
# also remove anything in uk10k or clinvar
set.seed(13457)
gnomad_processed_sub <- gnomad_processed %>% filter(!pos_id %in% (pos_id__source %>% filter(Source != 'gnomAD') %>% pull(pos_id))) %>% 
  sample_n((ML_set__clinvar %>% filter(Status=='Pathogenic') %>% nrow()) * 250)
gnomad_processed_sub_nonEye <- gnomad_processed %>% 
  filter(!pos_id %in% gnomad_processed_sub$pos_id) %>% 
  sample_n((ML_set__clinvar__otherPath_HC %>% filter(Status=='Pathogenic') %>% nrow()) * 250)
# the remainder gnomad variants
gnomad_processed_other <- gnomad_processed %>% filter(!pos_id %in% c(gnomad_processed_sub$pos_id,gnomad_processed_sub_nonEye$pos_id)) # not used for model building, for potential validation purposes

#gnomad_processed_sub$Source <- 'gnomAD'
#gnomad_processed_other$Source <- 'gnomAD'
print('gnomAD Loaded')

################################################
# Combine UK10K and ClinVar and gnomAD data
################################################

# ML_set__eye <- bind_rows(ML_set__clinvar %>% select_(.dots = colnames(ML_set__UK10K %>% select(-Complicated_Status, -Status))), 
#                          ML_set__UK10K %>% select(-Complicated_Status, -Status),
#                          gnomad_processed_sub %>% select_(.dots = colnames(ML_set__UK10K %>% select(-Complicated_Status, -Status)))) %>% 
#   group_by(pos_id) %>% filter(row_number() == 1) %>% ungroup()

ML_set__eye <- bind_rows(ML_set__clinvar %>% mutate_all(as.character), 
                         ML_set__UK10K %>% mutate_all(as.character),
                         gnomad_processed_sub %>% mutate_all(as.character)) %>% 
  mutate_at(vars(one_of(numeric_predictors)), funs(as.numeric(.))) %>% 
  group_by(pos_id) %>% filter(row_number() == 1) %>% ungroup()

# same as eye, but also adding HC clinvar and a bunch more gnomad
ML_set__general <- bind_rows(ML_set__clinvar %>% mutate_all(as.character),
                             ML_set__UK10K %>% mutate_all(as.character),
                             gnomad_processed_sub %>% mutate_all(as.character),
                             ML_set__clinvar__otherPath_HC %>% mutate_all(as.character), 
                             gnomad_processed_sub_nonEye %>% mutate_all(as.character)) %>% 
  mutate_at(vars(one_of(numeric_predictors)), funs(as.numeric(.))) %>%
  group_by(pos_id) %>% filter(row_number() == 1) %>% ungroup()

# the leftovers
ML_set__other <- bind_rows(gnomad_processed_other %>% mutate_all(as.character), 
                           ML_set__clinvar__otherPath %>% mutate_all(as.character)) %>% 
  mutate_at(vars(one_of(numeric_predictors)), funs(as.numeric(.))) %>%
  group_by(pos_id) %>% filter(row_number() == 1) %>% ungroup()

print('Core Sets Created')
################################
# one hot encode
##################################

one_hot_encode <- function(df){
  temp <- df %>% dplyr::select(-pos_id) %>% data.frame()
  temp <- dummy.data.frame(temp, sep='_')
  ML_set_dummy <- temp %>% mutate(pos_id = df$pos_id) %>% unique()
  ML_set_dummy
}

#ML_set__eye_dummy <- one_hot_encode(ML_set__eye)
#ML_set__general_dummy <- one_hot_encode(ML_set__general)
#ML_set__other_dummy <- one_hot_encode(ML_set__other)

#print('One Hot Encoding Done')
##################################
# add back status, remove dups
##################################

# create three columns of different pathogenic types
# pathogenic in eye
# pathogenic in clinvar, not so stringent
# pathogenic in clinvar, 3 or more stars (HC is 'high quality')
pos_id__source <- pos_id__source %>% 
  mutate(Path_Eye = case_when((Source == 'UK10K' | Source == 'ClinVar') & Status == 'Pathogenic' ~ 1, 
                              TRUE ~ 0),
         Path_ClinVar_Other = case_when(Source == 'ClinVar_OtherPath' & Status == 'Pathogenic' ~ 1,
                                        TRUE ~ 0),
         Path_ClinVar_Other_HC = case_when(Source == 'ClinVar_OtherPath_HC' & Status == 'Pathogenic' ~ 1,
                                           TRUE ~ 0)) 

# collapse by position id
pos_id__source <- pos_id__source %>% 
  group_by(pos_id) %>% 
  summarise(Comp_Status=paste(unique(Status), collapse=','), 
            Source=paste(Source,collapse=','),
            Path_Eye = max(Path_Eye),
            Path_ClinVar_Other = max(Path_ClinVar_Other),
            Path_ClinVar_Other_HC = max(Path_ClinVar_Other_HC)) %>% 
  mutate(Status = case_when(grepl('^Pathogenic|,Pathogenic', Comp_Status) ~ 'Pathogenic',
                            TRUE ~ 'NotPathogenic'))

print('Position ID table created')
# eye pathogenicity 
# select variants pathogenic in eye (either from ClinVar eye related or UK10K) or benign (clinvar, uk10k, gnomad)
ML_set__eye <- inner_join(ML_set__eye %>% select(-Status), 
                                pos_id__source %>% 
                                  filter((Path_Eye == 1) | Status=='NotPathogenic') %>% 
                                  dplyr::select(pos_id, Status, Source))
ML_set__eye$Status <- factor(ML_set__eye$Status, levels = c('Pathogenic','NotPathogenic'))
#ML_set__eye$DiseaseClass <- factor(ML_set__eye$DiseaseClass)

# general pathogenicity (can include eye)
ML_set__general <- inner_join(ML_set__general %>% select(-Status), 
                              pos_id__source %>% 
                                filter((Path_Eye == 1 | Path_ClinVar_Other_HC == 1) | Status=='NotPathogenic') %>% 
                                dplyr::select(pos_id, Status, Source))
ML_set__general$Status <- factor(ML_set__general$Status, levels = c('Pathogenic','NotPathogenic'))
#ML_set__general$DiseaseClass <- factor(ML_set__general$DiseaseClass)


# the remainder
ML_set__other <- inner_join(ML_set__other %>% select(-Status), pos_id__source %>% dplyr::select(pos_id, Status, Source))
ML_set__other$Status <- factor(ML_set__other$Status, levels = c('Pathogenic','NotPathogenic'))
#ML_set__other$DiseaseClass <- factor(ML_set__other$DiseaseClass)


# center scale 
# ML_set_dummy_CS <- preProcess(ML_set_dummy, method = c('center','scale')) %>% predict(., ML_set_dummy)

##################################
# train, validate, and test sets
# 70% to train
# 30% to test 
##################################

train_test_maker <- function(df){
  set.seed(115470)
  train_set <- df %>% 
    group_by(Status, Source) %>% 
    sample_frac(0.7) %>% ungroup() 
  
  test_set <- df %>% 
    filter(!pos_id %in% c(train_set$pos_id))
  
  out <- list()
  out$train_set <- train_set
  out$test_set <- test_set
  out
}

ML_set__eye_TT <- train_test_maker(ML_set__eye)
ML_set__general_TT <- train_test_maker(ML_set__general)
ML_set__other_TT <- train_test_maker(ML_set__other)

###########################################
# SAVE DATA
##########################################
model_data <- list()
model_data$ML_set__eye_TT <- ML_set__eye_TT
model_data$ML_set__general_TT <- ML_set__general_TT
model_data$ML_set__other_TT <- ML_set__other_TT
model_data$Test_set__UK10K <- Test_set__UK10K # withheld UK10K data for model perf 
model_data$pos_id__source <- pos_id__source
model_data$predictors <- numeric_predictors
model_data$sessionInfo <- sessionInfo()
save(model_data, file='/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/model_data_2018_08_01.Rdata')

#save(clinvar_spread, file='/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/model_spread.Rdata')
