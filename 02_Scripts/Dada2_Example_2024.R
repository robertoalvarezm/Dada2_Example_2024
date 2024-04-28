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



