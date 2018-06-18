################################################
# Process the gemini output data for the modeling
################################################

# biowulf paths
uk10k_data <- '/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/uk10k_gemini_rare_variants.Rdata'
clinvar_file <- '/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/clinvar/clinvar.gemini.tsv.gz'
gnomad_file <- '/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/gnomad_rare_benign_ish/gnomad_rare_benign_ish.gemini.tsv.gz'
# contains variants adjacent to existing variants in clinvar
spread_file <- '/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/clinvar/spread/clinvar.gemini.tsv.gz'

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
predictors <- c('pos_id', 
                'Status', 
                'Complicated_Status', 
                'is_exonic', 
                'is_coding', 
                'is_lof', 
                'is_splicing', 
                'impact_severity', 
                'polyphen_score', 
                'sift_score', 
                'dann', 
                'gerp_elements', 
                'DiseaseClass', 
                'mpc', 
                'revel', 
                'max_aaf_all',
                'gno_ac_afr', 
                'gno_ac_eas', 
                'gno_ac_all', 
                'gno_ac_popmax',
                'ac_exac_sas', 
                'ac_exac_fin',
                'aaf_1kg_all',
                'aaf_esp_all', 
                'ac_exac_all', 
                'ac_exac_amr', 
                'ac_exac_oth', 
                'gno_af_all', 
                'gno_an_popmax', 
                'an_exac_all', 
                'af_exac_all', 
                'fitcons', 
                'linsight', 
                'lof_z:precessive', 
                'phylop_100way',
                'syn_z',
                'grantham', 
                'cadd_phred', 
                'ccr_pct_v1', 
                'genesplicer', 
                'spliceregion',
                'eigen_pc_raw_rankscore',
                'fathmm_converted_rankscore',
                'genocanyon_score_rankscore',
                'lrt_converted_rankscore',
                'm_cap_rankscore',
                'metalr_rankscore',
                'metasvm_rankscore',
                'mutationassessor_score_rankscore',
                'mutationtaster_converted_rankscore',
                'provean_converted_rankscore',
                'vest3_rankscore',
                'fathmm_mkl_coding_rankscore',
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
                'epilogos_znf',
                'segway_mean_score',
                'segway_sum_score')

for_numeric <- 'ac_|an_|^n_|af_|dann|revel|mpc|gerp|polyphen_score|sift_score|fitcons|gerp_elements|^adj|_z$|^pli$|^pnull$|precessive|^phylop|linsight|_rankscore$|ccr_pct_v1|linsight|epilogos|segway'
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
  mutate_at(vars(matches(for_numeric)), funs(as.numeric(.))) %>%  # convert columns with ac_|whatever to integer (ac is allele count), etc. af is allele frequency
  select(one_of(predictors)) %>% 
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

ML_set__UK10K <- rbind(all_PATH, all_NOT_PATH__CUT)
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
  mutate_at(vars(matches('ac_|an_|^n_')), funs(as.integer(.))) %>% # convert columns with ac_|whatever to integer (ac is allele count)
  mutate_at(vars(matches('af_|dann|revel|mpc|gerp|polyphen_score|sift_score|fitcons|gerp_elements|^adj|_z$|^pli$|^pnull$|precessive|^phylop_100|linsight|_rankscore$|ccr_pct_v1')), funs(as.numeric(.))) %>%  # af is allele frequency
  selectmutate_at(vars(matches(for_numeric)), funs(as.numeric(.))) %>%  # convert columns with ac_|whatever to integer (ac is allele count), etc. af is allele frequency
  select(one_of(predictors)) %>%  
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
  filter(gene %in% (clinvar %>% filter(status=='PATHOGENIC_EYE' | status == 'PATHOGENIC_OTHER_HC') %>% pull(gene) %>% unique())) %>% # only have gnomad matched for path variants in clinvar genes
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
  mutate_at(vars(matches(for_numeric)), funs(as.numeric(.))) %>%  # convert columns with ac_|whatever to integer (ac is allele count), etc. af is allele frequency
  select(one_of(predictors)) %>% 
  filter(max_aaf_all < 0.01) # remove any common variants 


# fill missing with -1
gnomad_processed[is.na(gnomad_processed)] <- -1

pos_id__source <- rbind(pos_id__source, 
                        gnomad_processed %>% select(pos_id, Status) %>% mutate(Source='gnomAD')) %>% 
  unique()

# no more than 100X gnomad benign relative to the number of clinvar pathogenic variants 
set.seed(13457)
gnomad_processed_sub <- gnomad_processed %>% sample_n((ML_set__clinvar %>% filter(Status=='Pathogenic') %>% nrow()) * 100)
gnomad_processed_sub_nonEye <- gnomad_processed %>% 
  filter(!pos_id %in% gnomad_processed_sub$pos_id) #%>% 
  #sample_n((ML_set__clinvar__otherPath_HC %>% filter(Status=='Pathogenic') %>% nrow()) * 100)
# the remainder gnomad variants
gnomad_processed_other <- gnomad_processed %>% filter(!pos_id %in% c(gnomad_processed_sub$pos_id,gnomad_processed_sub_nonEye$pos_id)) # not used for model building, for potential validation purposes

#gnomad_processed_sub$Source <- 'gnomAD'
#gnomad_processed_other$Source <- 'gnomAD'
print('gnomAD Loaded')

# ################################################
# # Build ClinVar spread data
# ################################################
# spread <- fread(paste0('gzcat ', spread_file))
# 
# ## Prep data for modeling
# spread_processed <- spread %>% 
#   #filter(status!='PATHOGENIC_OTHER') %>% # drop non eye pathogenic variants for the model learning 
#   #separate(gene_eyediseaseclass, c('RDGene','DiseaseClass'), sep='_') %>%  #split off RD disease type
#   #select(-RDGene) %>% 
#   mutate(pos_id=paste0(chrom, ':', end, '_', ref, '_', alt),
#          impact_severity = case_when(impact_severity == 'HIGH' ~ 3, # convert to integer 
#                                      impact_severity == 'MED' ~ 2, 
#                                      TRUE ~ 1),
#          Status = factor('SPREAD'),
#          genesplicer = case_when(genesplicer == "" ~ 'No',
#                                  grepl('^gain', genesplicer) ~ 'Gain',
#                                  grepl('^loss', genesplicer) ~ 'Loss',
#                                  grepl('^diff', genesplicer) ~ 'Diff',
#                                  TRUE ~ 'Else')) %>% 
#   mutate_at(vars(matches('ac_|an_|^n_')), funs(as.integer(.))) %>% # convert columns with ac_|whatever to integer (ac is allele count)
#   mutate_at(vars(matches('af_|dann|revel|mpc|gerp|polyphen_score|sift_score|fitcons|gerp_elements|^adj|_z$|^pli$|^pnull$|precessive|^phylop_100|linsight|_rankscore$|ccr_pct_v1|epilogos|segway')), funs(as.numeric(.))) %>%  # af is allele frequency
#   select(pos_id,
#          Status,
#          is_exonic,
#          is_coding,
#          is_lof,
#          is_splicing,
#          impact_severity,
#          polyphen_score,
#          sift_score,
#          dann,
#          gerp_elements,
#          mpc,
#          revel,
#          max_aaf_all,
#          gno_ac_afr,
#          gno_ac_eas,
#          gno_ac_all,
#          gno_ac_popmax,
#          ac_exac_sas,
#          ac_exac_fin,
#          aaf_1kg_all,
#          aaf_esp_all,
#          ac_exac_all,
#          ac_exac_amr,
#          ac_exac_oth,
#          gno_af_all,
#          gno_an_popmax,
#          an_exac_all,
#          af_exac_all,
#          fitcons,
#          linsight,
#          lof_z:precessive,
#          phylop_100way,
#          grantham,
#          cadd_phred,
#          contains("_rankscore"),
#          contains("segway"),
#          contains("epilogos"),
#          ccr_pct_v1,
#          genesplicer,
#          spliceregion)
# 
# # fill missing with -1
# spread_processed[is.na(spread_processed)] <- -1
# 
# ML_set__spread <- spread_processed 
# 
# print('ClinVar Spread Loaded')

################################################
# Combine UK10K and ClinVar and gnomAD data
################################################

ML_set__eye <- bind_rows(ML_set__clinvar %>% select_(.dots = colnames(ML_set__UK10K %>% select(-Complicated_Status, -Status))), 
                         ML_set__UK10K %>% select(-Complicated_Status, -Status),
                         gnomad_processed_sub %>% select_(.dots = colnames(ML_set__UK10K %>% select(-Complicated_Status, -Status)))) %>% 
  group_by(pos_id) %>% filter(row_number() == 1) %>% ungroup()

# same as eye, but also adding HC clinvar and a bunch more gnomad
ML_set__general <- bind_rows(ML_set__clinvar %>% select_(.dots = colnames(ML_set__UK10K %>% select(-Complicated_Status, -Status))), 
                             ML_set__UK10K %>% select(-Complicated_Status, -Status),
                             gnomad_processed_sub %>% select_(.dots = colnames(ML_set__UK10K %>% select(-Complicated_Status, -Status))),
                             ML_set__clinvar__otherPath_HC %>% select_(.dots = colnames(ML_set__UK10K %>% select(-Complicated_Status, -Status))), 
                             gnomad_processed_sub_nonEye %>% select_(.dots = colnames(ML_set__UK10K %>% select(-Complicated_Status, -Status)))) %>% 
  group_by(pos_id) %>% filter(row_number() == 1) %>% ungroup()

# the leftovers
ML_set__other <- bind_rows(gnomad_processed_other %>% select(-Status), ML_set__clinvar__otherPath %>% select(-Status)) %>% 
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

ML_set__eye_dummy <- one_hot_encode(ML_set__eye)
ML_set__general_dummy <- one_hot_encode(ML_set__general)
ML_set__other_dummy <- one_hot_encode(ML_set__other)

print('One Hot Encoding Done')
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
ML_set__eye <- inner_join(ML_set__eye, 
                                pos_id__source %>% 
                                  filter((Path_Eye == 1) | Status=='NotPathogenic') %>% 
                                  dplyr::select(pos_id, Status, Source))
ML_set__eye$Status <- factor(ML_set__eye$Status, levels = c('Pathogenic','NotPathogenic'))
ML_set__eye$DiseaseClass <- factor(ML_set__eye$DiseaseClass)

ML_set__eye_dummy <- inner_join(ML_set__eye_dummy, 
                                pos_id__source %>% 
                                  filter((Path_Eye == 1) | Status=='NotPathogenic') %>% 
                                  dplyr::select(pos_id, Status, Source))
ML_set__eye_dummy$Status <- factor(ML_set__eye_dummy$Status, levels = c('Pathogenic','NotPathogenic'))


# general pathogenicity (can include eye)
ML_set__general <- inner_join(ML_set__general, 
                              pos_id__source %>% 
                                filter((Path_Eye == 1 | Path_ClinVar_Other_HC == 1) | Status=='NotPathogenic') %>% 
                                dplyr::select(pos_id, Status, Source))
ML_set__general$Status <- factor(ML_set__general$Status, levels = c('Pathogenic','NotPathogenic'))
ML_set__general$DiseaseClass <- factor(ML_set__general$DiseaseClass)

ML_set__general_dummy <- inner_join(ML_set__general_dummy, 
                                    pos_id__source %>% 
                                      filter((Path_Eye == 1 | Path_ClinVar_Other_HC == 1) | Status=='NotPathogenic') %>% 
                                      dplyr::select(pos_id, Status, Source))
ML_set__general_dummy$Status <- factor(ML_set__general_dummy$Status, levels = c('Pathogenic','NotPathogenic'))

# the remainder
ML_set__other <- inner_join(ML_set__other, pos_id__source %>% dplyr::select(pos_id, Status, Source))
ML_set__other$Status <- factor(ML_set__other$Status, levels = c('Pathogenic','NotPathogenic'))
ML_set__other$DiseaseClass <- factor(ML_set__other$DiseaseClass)

ML_set__other_dummy <- inner_join(ML_set__other_dummy, pos_id__source %>% dplyr::select(pos_id, Status, Source))
ML_set__other_dummy$Status <- factor(ML_set__other_dummy$Status, levels = c('Pathogenic','NotPathogenic'))

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

ML_set__eye_dummy_TT <- train_test_maker(ML_set__eye_dummy)
ML_set__general_dummy_TT <- train_test_maker(ML_set__general_dummy)
ML_set__other_dummy_TT <- train_test_maker(ML_set__other_dummy)

ML_set__eye_TT <- train_test_maker(ML_set__eye)
ML_set__general_TT <- train_test_maker(ML_set__general)
ML_set__other_TT <- train_test_maker(ML_set__other)

##################################
# ClinVar Spread
# train, validate, and test sets
# 70% to train
# 30% to test 
##################################
train_test_maker_SPREAD <- function(df){
  set.seed(115470)
  train_set <- df %>% 
    group_by(Status) %>% 
    sample_frac(0.7) %>% ungroup() 
  
  test_set <- df %>% 
    filter(!pos_id %in% c(train_set$pos_id))
  
  out <- list()
  out$train_set <- train_set
  out$test_set <- test_set
  out
}


clinvar_spread <- train_test_maker_SPREAD(ML_set__clinvar)
clinvar_spread$spread <- ML_set__spread

###########################################
# SAVE DATA
##########################################
model_data <- list()
model_data$ML_set__eye_dummy_TT  <- ML_set__eye_dummy_TT
model_data$ML_set__general_dummy_TT <- ML_set__general_dummy_TT
model_data$ML_set__other_dummy_TT <- ML_set__other_dummy_TT
model_data$ML_set__eye_TT <- ML_set__eye_TT
model_data$ML_set__general_TT <- ML_set__general_TT
model_data$ML_set__other_TT <- ML_set__other_TT
model_data$pos_id__source <- pos_id__source
#model_data$ML_set__spread <- ML_set__spread
model_data$sessionInfo <- sessionInfo()
save(model_data, file='/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/model_data.Rdata')

save(clinvar_spread, file='/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/model_spread.Rdata')
