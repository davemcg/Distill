# create gemini output, then merge together
library(tidyverse)
#gemini_tsvs <- list.files('/Volumes/data/projects/nei/mcgaughey/eye_var_Pathogenicity/data/', recursive = T, pattern='*tsv.gz', full.names = T)
base_path = '/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/'



gemini_dbs = c('ddl_nisc_100_panel/DDL_NISC_targeted_panel.PED_ddl_nisc.gemini.db',
          # 'gnomad_rare_benign_ish/gnomad.exomes.r2.0.2.sites.maxAF_0.01_20percent.PED_faux.gemini.db',
          'grimm/grimm.PED_faux.gemini.db',
          'homsy_benign/homsy.benign.PED_faux.gemini.db',
          'homsy_pathogenic/homsy.pathogenic.PED_faux.gemini.db',
          'ogvfb/VCFs.GATK.PED_master.gemini.db',
          'samocha_benign/samocha.benign.PED_faux.gemini.db',
          'samocha_pathogenic/samocha.pathogenic.PED_faux.gemini.db',
          'unifun_benign/UniFun_benign.fakeS.PED_faux.gemini.db',
          'unifun_pathogenic/UniFun_deleterious.fakeS.PED_faux.gemini.db',
          'clinvar/clinvar_RD.PED_faux.gemini.db')

for (gemini_db in gemini_dbs){
  system(paste0("time gemini query --header -q \"SELECT * from variants\" ", base_path, gemini_db, "| bgzip > ", base_path, gsub('.db','.tsv.gz', gemini_db)))
  print(paste(gemini_db, 'done'))
}

data_files = list()
for (i in gsub('.db','.tsv.gz', gemini_dbs)){
  name <- (i %>% str_split(., '/'))[[1]] %>% tail(2) %>% str_split(.,'\\.') %>% unlist() %>% head(1)
  print(name)
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

save(raw_data, file='/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/master/raw_data_2018_07_13.Rdata')
