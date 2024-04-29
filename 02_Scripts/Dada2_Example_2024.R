library(dada2)
path<-"01_Raw_Data/"
list.files(path)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names


## Verifcar calidad PHRED

pdf("03_Results/Quality_Forward.pdf",width=13,height = 8)
plotQualityProfile(fnFs[1:10])
dev.off()


pdf("03_Results/Quality_Reverse.pdf",width=13,height = 8)
plotQualityProfile(fnRs[1:10])
dev.off()


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

# Si tu máquina es muy lenta ya generé este objeto y solo cárgalo con 
# readRDS(file="03_Results/errF.RDS)
errF <- learnErrors(filtFs, multithread=TRUE)
#saveRDS(errF,file="03_Results/errF.RDS")

# Si tu máquina es muy lenta ya generé este objeto y solo cárgalo con 
# readRDS(file="03_Results/errR.RDS)
errR <- learnErrors(filtRs, multithread=TRUE)
#saveRDS(errR,file="03_Results/errR.RDS")

png("03_Results/Errores_Forward.png")
plotErrors(errF, nominalQ=TRUE)
dev.off()

png("03_Results/Errores_Reverse.png")
plotErrors(errR, nominalQ=TRUE)
dev.off()


dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
saveRDS(mergers, file="03_Results/mergers.RDS")
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
saveRDS(seqtab, file="03_Results/seqtab.RDS")


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


taxa <- assignTaxonomy(seqtab.nochim, "01_Raw_Data/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
saveRDS(taxa, file="03_Results/taxa.RDS")
#readRDS("03_Results/taxa.RDS")


taxa <- addSpecies(taxa, "01_Raw_Data/silva_species_assignment_v138.1.fa.gz")
saveRDS(taxa, file="03_Results/taxa.RDS")
#readRDS("03_Results/taxa.RDS")



taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
saveRDS(taxa.print, file="03_Results/taxa.print.RDS")
#readRDS("03_Results/taxa.print.RDS")


library(DECIPHER)

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("01_Raw_Data/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; 
rownames(taxid) <- getSequences(seqtab.nochim)
saveRDS(taxid,file="03_Results/taxid.RDS")
readRDS(file="03_Results/taxid.RDS")


unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

