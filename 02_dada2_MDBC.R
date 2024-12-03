### NOTES ###

# Approximate expected sizes for on-target ASVs
# 18S: 125-129 bp
# CO1: 313 bp
# 28S: 333-371 bp
# 12S: ~170 bp


### END NOTES ###

### END COMMAND LINE ###

# load libraries
library(dada2)
library(Biostrings)
library(decontam)

# Modified from: https://benjjneb.github.io/dada2/tutorial.html

# Importing reads ----
### CHANGE WORKING PATHS
# path to folder where results will be saved
setwd("/path/to/folder")
# path to raw demultiplexed reads
path <- "path/to/demultiplexed/reads"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnFs <- sort(list.files(path, pattern="_R1.clean.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.clean.fq.gz", full.names = TRUE))

# return sample names without _R1/_R2 or fastq.gz
# this assumes sample names do NOT have an underscore in them!
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Filter and trim ----

# Create subdirectory named filtered and place files in that directory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

### CHANGE TRUNCATION LENGTH
# first value is R1, second is R2
# these might need to be changed depending on sequencing
# 12S/Mifish: truncLen=100
# 28S/corals: truncLen=225
# CO1: truncLen=c(200,150)
# 18S: truncLen=100
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, truncLen = c(100,100), maxEE=c(2,2), rm.phix=TRUE, compress=TRUE, multithread=TRUE)

# remove files with zero reads (otherwise they will cause an error during the inference step)
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]

# Learn error rates ----
# machine learning of error rates in the data set (this may take a little while...)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Infer number of true variants using pooled samples (see https://benjjneb.github.io/dada2/pseudo.html)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = TRUE)

# Merge paired reads ----
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, trimOverhang = TRUE, verbose=TRUE)

# make a sequence table from merged reads
seqtab <- makeSequenceTable(mergers)

# Remove contaminants ----
# See package `decontam` package here: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

### import metadata
# the following code expects the metadata to contain a column named Blank
# which has TRUE or FALSE for each sample:
# blanks = TRUE
# real samples = FALSE
meta <- read.table("SALT3_metadata.tsv", header = TRUE, sep="\t")

# Remove samples from metadata that were removed from analyses
samples <- as.data.frame(names(mergers))
colnames(samples) <- "sample_ID"
meta <- merge(samples, meta, by="sample_ID", all=FALSE)

# extract TRUE/FALSE for blanks
vector_for_decontam <- meta$Blank

# table identifying whether sequences are contaminants or not
contam_df <- isContaminant(seqtab, neg=vector_for_decontam, threshold=0.1)

# stats: sum of total ASVs that are contaminants
sum(contam_df$contaminant == TRUE)

# stats: percentage of ASVs that are contaminants
plyr::count(contam_df, vars = "contaminant")

# pull out ASV sequences that are potential contaminants
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

# remove contaminants from data
seqtab <- seqtab[,!colnames(seqtab) %in% contam_asvs]

# export contaminants to fasta
contam_seq <- Biostrings::DNAStringSet(contam_asvs)
names(contam_seq) <- paste0("contam", 1:nrow(as.matrix(contam_asvs)))
Biostrings::writeXStringSet(contam_seq, "ASV_contaminants.fasta", append=FALSE, compress=FALSE, format="fasta")

# END CONTAMINANTS

# Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# stats: frequency of chimeric sequences
sum(seqtab.nochim)/sum(seqtab)


# Make stats table ----
library(plyr) #needed for join_all and count

getN <- function(x) sum(getUniques(x))

# prepare data frames with statistics for merging
filtered <- as.data.frame(out)
denoisedF <- as.data.frame(sapply(dadaFs, getN))
denoisedR <- as.data.frame(sapply(dadaRs, getN))
merged <- as.data.frame(sapply(mergers, getN))
contams <- as.data.frame(rowSums(seqtab))
nonchim <- as.data.frame(rowSums(seqtab.nochim))

# make row names into a column
filtered$rn <- rownames(filtered)
denoisedF$rn <- rownames(denoisedF)
denoisedR$rn <- rownames(denoisedR)
merged$rn <- rownames(merged)
contams$rn <- rownames(contams)
nonchim$rn <- rownames(nonchim)

# clean up row names in the `filtered` data frame so they match the sample names
# in the rest of the data frames
filtered <- transform(filtered, rn = sub("_R1.clean.fq.gz", "", rn))

# merge the data frames and clean up
track <- join_all(list(filtered,denoisedF,denoisedR,merged,contams,nonchim), by = 'rn', type = 'full')
rownames(track) <- track$rn
track <- subset(track, select = -rn)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "decontam", "nonchim")
write.csv(track, "dada2_denoising_stats.csv", row.names = TRUE)

# extract and rename ASVs
ASVs.nochim <- Biostrings::DNAStringSet(colnames(seqtab.nochim))
names(ASVs.nochim) <- paste0("ASV", 1:ncol(seqtab.nochim))

# export ASVs to fasta
Biostrings::writeXStringSet(ASVs.nochim, "ASVs.fasta", append=FALSE, compress=FALSE, format="fasta")

# Create a table of counts for each ASV sequence and export
seq <- seqtab.nochim
colnames(seq) <- names(ASVs.nochim)
asv.tab <- t(seq)
write.table(asv.tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

