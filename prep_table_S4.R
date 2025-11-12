library(tidyverse)
library(jsonlite)
library(ape)

df = read_csv("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\Bella\\Bioinformatics\\spo0B_TM_110525\\TMR_df_all_111025.csv") %>%
  select(-...1)

df_gcf = read_tsv("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\Bella\\Bioinformatics\\spo0B_TM_110525\\spo0B_genome_hits_111025_2.tsv", col_names = c("genome_acc", "prot_acc")) 

#%>%
  #distinct(genome_acc, .keep_all = TRUE)

df_full = full_join(df, df_gcf,by = join_by(accession == prot_acc))  %>%
  filter(aa_length < 300)

df_bacilli = read_tsv("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\Bella\\Bioinformatics\\spo0B_TM_110525\\bacilli.tsv")%>% left_join(df_full, by = join_by(`Assembly Accession` == genome_acc))

file = "C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\Bella\\Bioinformatics\\spo0B_TM_110525\\spo0B_pHMMER_5e2_111025.csv"

df_hmmer = data.frame(readLines(file)) %>%
  slice(4:n()) %>%
  separate(readLines.file., into = c(as.character(1:34)), sep = " +") %>%
  select(-2,-3,-4) %>%
  replace(is.na(.), " ") %>% #Replace NAs with spaces.
  unite("description", 16:31, sep = " ") %>%
  mutate_if(is.character, str_trim) %>% #Trim spaces off description.
  filter(`1` != "#") %>%
  mutate(species = str_extract(description, "(?<=\\[).*(?=\\])")) %>% #Extract species name from between [] in description. 
  `colnames<-`(c("accession", "full_eval", "full_score", "full_bias", "dom_eval", "dom_score", "dom_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description", "species")) %>% #Rename columns.
  mutate_at(vars(full_eval, full_score, full_bias, dom_eval, dom_score, dom_bias), as.numeric) 


df_hmmer_join = right_join(df_bacilli, df_hmmer, by = join_by(accession))


df_hmmer_best_hit = df_hmmer_join %>% 
  group_by(`Assembly Accession`) %>%
  arrange(full_eval) %>%
  filter(row_number()==1) 


#write.csv(df_hmmer_best_hit$accession, "C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\Bella\\Bioinformatics\\spo0B_TM_110525\\spo0B_best_hits_111025.txt", row.names = FALSE, quote = FALSE)

##############

lineages = read_tsv("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\Bella\\Bioinformatics\\spo0B_TM_110525\\bac_taxonomy_all.tsv")

df_lineages = lineages %>%
  select(-Query, -Authority, -Rank, -`Basionym authority`, -`Curator common name`, -`Has type material`, -`Kingdom taxid`, -`Domain taxid`, -`Phylum taxid`, -`Class taxid`, -`Order taxid`, -`Family taxid`, -`Species taxid`)


df_final = left_join(df_bacilli, df_lineages, by = join_by(`Organism Taxonomic ID` == `Taxid`)) %>%
  mutate(spo0B_presence = ifelse(is.na(accession), 0, 1)) %>%
  mutate(TMR_presence = ifelse(TMRs >0, 1, 0))%>%
  mutate(TMRs = replace(TMRs, is.na(TMRs), "no_spo0B"))
  

df_spo0B = df_final %>%
  filter(spo0B_presence == 1) 

table(df_spo0B$`Family name`, df_spo0B$TMR_presence)

table(df_final$`Family name`, df_final$TMR_presence)

summarize(df_spo0B, total_TMR = sum(TMR_presence), total_spo0B = sum(spo0B_presence), perc_TMR = sum(TMR_presence)/sum(spo0B_presence))


df_spo0B %>%
  group_by(`Family name`) %>%
  summarize(total_TMR = sum(TMR_presence), total_spo0B = sum(spo0B_presence), perc_TMR = sum(TMR_presence)/sum(spo0B_presence))


df_spo0B %>%
  filter(`Family name` == "Paenibacillaceae") %>%
  group_by(`Genus name`) %>%
  summarize(total_TMR = sum(TMR_presence), total_spo0B = sum(spo0B_presence), perc_TMR = sum(TMR_presence)/sum(spo0B_presence))

df_spo0B %>%
  filter(`Family name` == "Alicyclobacillaceae") %>%
  group_by(`Genus name`) %>%
  summarize(total_TMR = sum(TMR_presence), total_spo0B = sum(spo0B_presence), perc_TMR = sum(TMR_presence)/sum(spo0B_presence))


table_S3 = df_final %>%
  select(-`Genus taxid`) %>%
  unite("taxonomy", `Order name`:`Species name`, sep = ";") %>%
  select(`Assembly Accession`, `Assembly Name`, taxonomy, spo0B_presence, TMR_presence) %>%
  mutate(TMR_presence = replace(TMR_presence, is.na(TMR_presence), "no_spo0B"))
  

write.csv(table_S3, "C:\\Users\\cassp\\Cornell University\\Heather Feaga - Bella Paenibacillus paper\\Bella Paenibacillus paper\\Table_S3_Bacilli.csv", row.names = FALSE, quote = FALSE) 
