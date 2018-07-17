setwd('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/UK10K')
#system('bash ~/git/eye_var_Pathogenicity/scripts/gemini_query_calls_UK10K.sh')


library(tidyverse)
metadata <- readxl::read_excel(path='~/git/EGA_EGAD00001002656_NGS_reanalyze/data/1-s2.0-S0002929716305274-mmc3.xlsx') %>% mutate(Sample=Patient)

library(data.table)
het <- fread('gzcat /data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/UK10K/UK10K_EGAD.hets.gz')
hom <- fread('gzcat /data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/data/UK10K/UK10K_EGAD.homs.gz')

hom_long <- hom %>% 
  mutate(Sample = str_split(hom_alt_samples, ","), Genotype = 'Hom') %>%  # split with str_split to create a list
  unnest(Sample, .drop=F) # unnest will spread the list, creating a row for each entry

het_long <- het %>% 
  mutate(Sample = str_split(het_samples, ","), Genotype = 'Het') %>%  # split with str_split to create a list
  unnest(Sample, .drop=F) # unnest will spread the list, creating a row for each entry

all <- rbind(hom_long, het_long)

sample_gene_comp_het <- metadata %>% filter(Status=='Solved'  & Variant_HGVSc!='NA' & GT=='0/1') %>% group_by(Sample, Gene) %>% summarise(Count=n()) %>% filter(Count>1) 
metadata <- left_join(metadata, sample_gene_comp_het) %>% mutate(Comp_Het_Path = case_when(Count >= 2 ~ 'CH', 
                                                                                           TRUE ~ 'No')) %>% 
  select(-Count)

uk10k_gemini_rare_variants <- all %>% 
  mutate(Variant_genomic = paste0(chrom, ':', start + 1, ref, '>', alt)) %>% 
  mutate(Status = case_when(Variant_genomic %in% (metadata %>% filter(Status == 'Solved') %>% pull(Variant_genomic)) ~ 'Pathogenic',
                            Variant_genomic %in% (metadata %>% filter(Status == 'Partially solved') %>% pull(Variant_genomic)) ~ 'Maybe Pathogenic',
                            TRUE ~ 'NotPathogenic')) %>% 
  filter(Status != 'Maybe Pathogenic') %>% 
  mutate(Complicated_Status = case_when(Variant_genomic %in% (metadata %>% filter(Status != 'Unsolved' & Comp_Het_Path == 'CH') %>% pull(Variant_genomic)) ~ 'Comp_Het',
                                        Variant_genomic %in% (metadata %>% filter(Status != 'Unsolved' & GT=='0/1') %>% pull(Variant_genomic)) ~ 'AD',
                                        Variant_genomic %in% (metadata %>% filter(Status != 'Unsolved' & GT=='1/1') %>% pull(Variant_genomic)) ~ 'AR',
                                        TRUE ~ 'NotPathogenic')) %>% 
  nest(Sample, .key='Samples')


output_file <- '/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/uk10k_gemini_rare_variants.Rdata'
#if(!file.exists(output_file)){
  save(uk10k_gemini_rare_variants, file = output_file)
#}

