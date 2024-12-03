# This script was used to process additional 12S sequences from GenBank 

library(tidyverse)
library(taxonomizr)
library(Biostrings) # DNA sequence manipulation tools
library(rBLAST) # BLAST searches without leaving R
library(refdb) # reference database cleaning
library(here)

## Get NCBI sequences ----------------------------------------------------------
# include: txid7742 (vertebrata)
# exlude: txid32523 (tetrapoda)

# COMMAND LINE:
# retrieve mitochondrial sequences from GenBank
# all vertebrates except tetrapods
# esearch -db nucleotide -query '("txid7742"[ORGN] AND "mitochondrion"[filter] AND ("150"[SLEN] : "50000"[SLEN])) NOT "txid32523"[ORGN] NOT "PREDICTED"[Title] NOT "UNVERIFIED"[Title]' | \
# efetch -format fasta > fish_seqs_NCBI_2024-08-22.fa

# Convert to BLAST db ----
rBLAST::makeblastdb(file = "NCBI_fish_mitogenomes/fish_seqs_NCBI_2024-08-22.fa", dbtype = "nucl")

blast_db_path <- here::here("NCBI_fish_mitogenomes/fish_seqs_NCBI_2024-08-22.fa")

# Look for primers in sequences ------------------------------------------------
## define primer sequences (MiFish 12S) ----
forward_primer <- "GTGTCGGTAAAACTCGTGCCAGC"
reverse_primer <- "CATAGTGGGGTATCTAATCCCAGTTTG"

# get reverse complement of reverse primer (3' --> 5')
reverse_primer_reversecomp <- bioseq::as_dna(reverse_primer) %>%
  bioseq::seq_complement() %>%
  bioseq::seq_reverse()

## import sequences ----
# as_dna
seqs <- bioseq::read_fasta(file = blast_db_path, type = "DNA")

## extract amplicons using primer sequences ---- 
# allow 15 % error in primers
# return without primer portion of sequence
# this can take a long time
seed_seqs <- 
  bioseq::seq_crop_pattern(
  seqs,
  forward_primer,
  reverse_primer_reversecomp,
  max_error_in = 0.15,
  max_error_out = 0.15,
  include_patterns = FALSE
  # remove sequences that don't have matches 
) %>% na.omit()

## filter seed sequences by length -----
# get lengths of sequences
seq_lengths <- as.data.frame(bioseq::seq_nchar(seed_seqs))
# bind lengths to sequences
seed_seqs_df <- as.data.frame(seed_seqs)
# add accession/description to df
seed_seqs_df$description <- names(seed_seqs)
# clean up the results
seed_seqs_final <- cbind(seq_lengths, seed_seqs_df) %>%
  # move rownames to ID column
  tibble::rownames_to_column() %>%
  tidyr::separate(rowname, into = c("accession", "description"), sep = "\\s", extra = "merge") %>%
  # prettify names
  dplyr::rename("length" = "bioseq::seq_nchar(seed_seqs)", "sequence" = seed_seqs) %>%
  # filter sequences by approximate amplicon length
  dplyr::filter(length < 190, length > 165)

# blast seeds against database ----
# prep for blastn
blastn <- rBLAST::blast(db = blast_db_path, type = "blastn")

# run BLAST
blast_results <- predict(blastn, Biostrings::DNAStringSet(seed_seqs_final$sequence), custom_format = "qseqid sseqid sstart send length qcovs bitscore")

# remove duplicates and filter results
blast_results_filtered <- blast_results %>%
  # remove qseqid
  dplyr::select(-qseqid) %>%
  # filter sequences by length
  dplyr::filter(length < 190, length > 165) %>%
  # sort by longest length
  # there may be duplicate entries and this will ensure only the longest are kept
  dplyr::arrange(desc(bitscore)) %>%
  # remove duplicate entries
  dplyr::distinct(sseqid, .keep_all = TRUE)

# Add taxonomy to NCBI accessions ----------------------------------------------

# grab accessions
accession <- blast_results_filtered$sseqid
taxaId <- taxonomizr::accessionToTaxa(accession,"/path/to/accessionTaxa.sql")

# get taxonomy
taxa <- taxonomizr::getTaxonomy(taxaId,"/path/to/accessionTaxa.sql")

# convert sequences to data frame
seqs_df <- as.data.frame(seqs)
seqs_df$accession <- names(seqs) 
seqs_df <- seqs_df %>%
  # split accession (first value before space) into one column and description into another column
  tidyr::separate(accession, into = c("accession", "description"), sep = "\\s", extra = "merge") %>%
  dplyr::rename("sequence" = "seqs")

# combine accessions, taxa, and blast results
taxa_accessions <- row.names(taxa) %>%
  cbind(accession,taxa) %>%
  tibble::as_tibble() %>%
  dplyr::rename("taxaID" = ".") %>%
  dplyr::mutate(taxaID = as.numeric(taxaID)) %>%
  # add blast results
  dplyr::left_join(blast_results_filtered, by = c("accession" = "sseqid")) %>%
  # add sequences
  dplyr::left_join(seqs_df, by = "accession") %>%
  # check for reversed sequences
  dplyr::mutate(start = dplyr::case_when(sstart > send ~ send, .default = sstart),
                end = dplyr::case_when(sstart > send ~ sstart, .default = send),
                reversed = dplyr::case_when(sstart > send ~ TRUE, .default = FALSE))

# Crop sequences by position ---------------------------------------------------

# initiate tibble for cropped sequences
cropped_sequences <- tibble::tibble(this_accession = character(), sequence = character())

# extract sequence based on BLAST position
for (i in 1:nrow(taxa_accessions)) {
  
  this_accession <- taxa_accessions[i,]$accession
  this_sstart <- taxa_accessions[i,]$start
  this_ssend <- taxa_accessions[i,]$end
  this_seq <- taxa_accessions[i,]$sequence
  
  sequence <- stringr::str_sub(this_seq, start = this_sstart, end = this_ssend)
  new_seq <- cbind(this_accession, sequence)
  cropped_sequences <- rbind(cropped_sequences, new_seq)
  
}

# merge cropped sequences with the rest of the data
sequence_database <- taxa_accessions %>%
  # remove uncropped sequence
  dplyr::select(-c(taxaID, sequence, sstart, send, bitscore, qcovs)) %>%
  # add croped sequence
  dplyr::left_join(cropped_sequences, by = c("accession" = "this_accession"))

# pull out sequences that need to be reversed
sequence_database_reversed <- sequence_database %>%
  dplyr::filter(reversed == TRUE)

# reverse complement sequences
reverse_seq <- bioseq::as_dna(sequence_database_reversed$sequence) %>%
  bioseq::seq_complement() %>%
  bioseq::seq_reverse() %>%
  tibble::as_tibble()

sequence_database_revcomp <- cbind(sequence_database_reversed, reverse_seq) %>%
  # remove old sequence
  dplyr::select(-sequence) %>%
  # rename new sequence column
  dplyr::rename("sequence" = value)

# merge sequences back with data
sequence_database_fixed <- sequence_database %>%
  dplyr::filter(reversed == FALSE) %>%
  dplyr::bind_rows(sequence_database_revcomp)

# Clean reference database -----------------------------------------------------
ref_db <- refdb::refdb_set_fields(
  sequence_database_fixed,
  id = "accession",
  taxonomy = c(superkingdom = "superkingdom",
               phylum = "phylum",
               class = "class",
               order = "order",
               family = "family",
               genus = "genus",
               species = "species"),
  sequence = "sequence"
)

# remove extra info in species name
ref_db <- refdb::refdb_clean_tax_remove_extra(ref_db)
ref_db <- refdb::refdb_clean_tax_remove_uncertainty(ref_db)
ref_db <- refdb::refdb_clean_tax_remove_subsp(ref_db)

# sp. --> NA
ref_db <- refdb::refdb_clean_tax_NA(ref_db)

# remove any extra words
ref_db <- ref_db %>%
  dplyr::mutate(species = stringr::word(species, start = 1, end = 2))

sequence_conflicts <- refdb::refdb_check_seq_conflict(ref_db)

## Save clean database as fasta ----

# format for fasta
# remove rows with no sequence data
ref_db_fasta <- ref_db %>%
  dplyr::filter(!is.na(sequence)) %>%
  dplyr::mutate(kingdom = "Animalia") %>%
  dplyr::rename("domain" = superkingdom) %>%
  dplyr::select(domain, kingdom, phylum, class, order, family, genus, species, sequence) %>%
  # merge taxonomy into singel column
  dplyr::mutate(species = stringr::str_replace_all(species, " ", "_")) %>%
  tidyr::unite(domain:species, col = "taxonomy", sep = ";", remove = FALSE) %>%
  # remove duplicates
  dplyr::distinct()

# convert sequences to string set
new_seqs <- Biostrings::DNAStringSet(ref_db_fasta$sequence)
# associate sequences with other info
names(new_seqs) <- ref_db_fasta$taxonomy
# save sequences to fasta file
Biostrings::writeXStringSet(new_seqs, here::here("new_12S_Mifish_sequences.fasta"), append=FALSE, compress=FALSE, format="fasta")
