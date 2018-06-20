# load raw
library(tidyverse)
#gemini_tsvs <- list.files('/Volumes/data/projects/nei/mcgaughey/eye_var_Pathogenicity/data/', recursive = T, pattern='*tsv.gz', full.names = T)
base_path = '/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/'
files = c('ddl_nisc_100_panel/ddl_nisc_100_panel.gemini.tsv.gz',
          'gnomad_rare_benign_ish/gnomad_rare_benign_ish.gemini.tsv.gz',
          'grimm/grimm.gemini.tsv.gz',
          'homsy_benign/homsy_benign.gemini.tsv.gz',
          'homsy_pathogenic/homsy_pathogenic.gemini.tsv.gz',
          'samocha_benign/samocha_benign.gemini.tsv.gz',
          'samocha_pathogenic/samocha_pathogenic.gemini.tsv.gz',
          'unifun_benign/unifun_benign.gemini.tsv.gz',
          'unifun_pathogenic/unifun_pathogenic.gemini.tsv.gz',
          'clinvar/clinvar.gemini.tsv.gz')

data_files = list()
for (i in files){
  name <- (i %>% str_split(., '/'))[[1]] %>% tail(1) %>% str_split(.,'\\.') %>% unlist() %>% head(1)
  data_files[[name]] <- read_tsv(paste0(base_path, i), col_types = cols(.default = "c"))
}
het <- read_tsv('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/UK10K/UK10K_EGAD.hets.gz', col_types = cols(.default = "c"))
hom <- read_tsv('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/UK10K/UK10K_EGAD.homs.gz', col_types = cols(.default = "c"))

raw_data <- bind_rows(data_files, .id = 'DataSet')
raw_data <- bind_rows(raw_data, 
                      het %>% mutate(DataSet = 'UK10K'), 
                      hom %>% mutate(DataSet = 'UK10K'))
#rm(data_files)
raw_data <- raw_data %>% mutate(pos_id=paste0(chrom, ':', end, '_', ref, '_', alt))

save(raw_data, file='/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/master/raw_data.Rdata')
