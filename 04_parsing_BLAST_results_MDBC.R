library(tidyverse)

# BLAST ASVs
path <- "path/to/results"
marker <- "28S"

# set blast database path
blast_db_path <- "path/to/NCBI_MZGdb_fish_12S_sequence_database.fa"

# define blastn
blastn <- rBLAST::blast(db = blast_db_path, type = "blastn")

# read in ASV fasta
ASV_seqs <- Biostrings::readDNAStringSet(paste0(path, "/ASVs_MDBC_", marker,".fasta"), format="fasta", use.names = TRUE)

# run blast
blast_results <- predict(blastn, ASV_seqs, custom_format = "qseqid sseqid pident length qcovs mismatch gapopen evalue bitscore")

# Import and clean BLAST results ----

# column names for blast results
# 28S
col_names <- c("ASV" = "qseqid",
               "Accession" = "sseqid1",
               "Domain" = "sseqid2",
               "Kingdom" = "sseqid3",
               "Phylum" = "sseqid4",
               "Class" = "sseqid5",
               "Order" = "sseqid6",
               "Family" = "sseqid7",
               "Genus" = "sseqid8",
               "Species" = "sseqid9",
               "Identity" = "pident",
               "Length" = "length",
               "Coverage" = "qcovs",
               "Mismatch" = "mismatch",
               "Gap" = "gapopen",
               "E-value" = "evalue",
               "Bitscore" = "bitscore")

#18S/CO1
col_names <- c("ASV" = "qseqid",
               "Kingdom" = "sseqid1",
               "Phylum" = "sseqid2",
               "Class" = "sseqid3",
               "Order" = "sseqid4",
               "Family" = "sseqid5",
               "Genus" = "sseqid6",
               "Species" = "sseqid7",
               "Identity" = "pident",
               "Length" = "length",
               "Coverage" = "qcovs",
               "Mismatch" = "mismatch",
               "Gap" = "gapopen",
               "E-value" = "evalue",
               "Bitscore" = "bitscore")
# 12S
col_names <- c("ASV" = "qseqid",
               "Accession" = "sseqid1",
               "Domain" = "sseqid2",
               "Kingdom" = "sseqid3",
               "Phylum" = "sseqid4",
               "Class" = "sseqid5",
               "Order" = "sseqid6",
               "Family" = "sseqid7",
               "Genus" = "sseqid8",
               "Species" = "sseqid9",
               "Identity" = "pident",
               "Length" = "length",
               "Coverage" = "qcovs",
               "Mismatch" = "mismatch",
               "Gap" = "gapopen",
               "E-value" = "evalue",
               "Bitscore" = "bitscore")

# separate taxonomy out into separate columns using delimiter (semi-colon)
# automatically names the column by adding a number to end of "taxonomy"
blast_clean <- blast_results %>% 
  separate_wider_delim(sseqid, delim = ";", names_sep = "", too_few = "align_start") %>%
  dplyr::rename(all_of(col_names)) %>%
  mutate(Species = str_replace_all(Species, "_", " "),
         Species = case_when(stringr::str_detect(Species, "EXT") ~ NA_character_,
                             .default = Species),
         Species = case_when(stringr::str_detect(Species, "NA") ~ NA_character_,
                             .default = Species)) %>%
  # change the columns that are numbers to numeric
  # otherwise the filtering below will not work correctly
  mutate_at(c("Identity", "Length", "Coverage","Mismatch","Gap","E-value","Bitscore"), as.numeric)
write_tsv(blast_clean, paste0(path, "/blast_results_", marker,".tsv"))


## Modified LCA ----
# subset species
# return only hits that are identified to species level with at least 99% identity and 90% coverage
species <- blast_clean %>%
  filter(!is.na(Species) & Identity >= 99 & Coverage >= 90 & Bitscore > 100) %>%
  dplyr::group_by(ASV) %>%
  dplyr::filter(Bitscore > max(Bitscore)-(max(Bitscore)*0.05)) %>%
  ungroup()

# return list of species with multiple good hits
problematic_species <- species %>%
  group_by(ASV) %>%
  # return only hits that have two different species per ASV
  distinct(Species, .keep_all = TRUE) %>%
  # return ONLY ASVs that have multiple "best" blast hits
  filter(n()>1)
write_tsv(problematic_species, paste0(path, "/problematic_matches_", marker,".tsv"))


# Filter by bitscore
# Filter only bitscores within 5% of the maximum bitscore for each ASV
not_species <- blast_clean %>%
  filter(!ASV %in% unique(species$ASV)) %>%
  filter(!is.na(Species) & Identity >= 90 & Coverage >= 90 & Bitscore > 100) %>%
  dplyr::group_by(ASV) %>%
  dplyr::filter(Bitscore > max(Bitscore)-(max(Bitscore)*0.05)) %>%
  select(ASV, Accession, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  ungroup()

# combine species and not_species for LCA
final <- bind_rows(species, not_species) %>%
  select(ASV, Accession, Kingdom, Phylum, Class, Order, Family, Genus, Species)

# run LCA function
final_lca <- aggregate(
  final[, 3:ncol(final)], 
  by = list(ASV = final$ASV), 
  FUN = function(taxonomy_levels) {
    lca <- Reduce(intersect, taxonomy_levels)
    return(ifelse(length(lca) == 0, NA, paste(lca, collapse = NULL))) 
  }
)

# Add a step to aggregate the second column using a comma as a delimiter
accession_column <- aggregate(
  final[, 2], 
  by = list(ASV = final$ASV), 
  FUN = function(info_values) {
    paste(unique(info_values), collapse = ",")
  }
)

# get names and sequences for all ASVs
all_ASVs <- as.data.frame(ASV_seqs) %>%
  rownames_to_column(var = "ASV") %>%
  dplyr::rename("sequence" = x)
  
# rename the first column back to ASV and save
final_lca_clean <- merge(final_lca, accession_column, by = "ASV") %>%
  full_join(all_ASVs, by = "ASV") %>%
  # replace stray NAs with NA character
  mutate(Species = str_replace_all(Species, "NA", NA_character_)) %>%
  arrange(Kingdom, Phylum, Class, Order, Family, Genus, Species)
write_tsv(final_lca_clean, paste0(path, "/lca_", marker,".tsv")) 

# output taxa list
final_taxa_list <- final_lca_clean %>%
  select(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  dplyr::distinct(.keep_all = FALSE) %>%
  arrange(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  drop_na(Kingdom)
write_tsv(final_taxa_list, paste0(path, "/taxa_list_", marker,".tsv"))


