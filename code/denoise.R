#denoise samples using dada2 
#Modified from DADA2 tutorial (https://benjjneb.github.io/dada2/tutorial.html)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#Steps marked with alternating hash/dash may require adjustment for each data set

#load packages
library(tidyverse)
library(magrittr)
library(dada2)
library(ShortRead)
library(phyloseq)
library(Biostrings)

# Get input/output directories from script arguments (or use defaults)
if(length(commandArgs(trailingOnly=T))==0){
  inDir <- "output/processing/trim"
  outDir <- "output/processing/dada"
}else{
  inDir <- commandArgs(trailingOnly=T)[1]
  outDir <- commandArgs(trailingOnly=T)[2]
}
# Try to create output directory in case it does not exist
dir.create(outDir, showWarnings = F) 

#Identify trimmed fwd and rev reads
path.fwd <- inDir %>% list.files(pattern="R1.fq.gz", full.names = T) %>% sort 
path.rev <- inDir %>% list.files(pattern="R2.fq.gz", full.names = T) %>% sort 
  
#preview read quality of 12 randomly selected samples before trimming
qual.samps <- sample(1:length(path.fwd),9)
(qual.fwd <- plotQualityProfile(path.fwd[qual.samps]) + ggtitle("Qual profiles fwd"))
file.path(outDir, "qual.fwd.pdf") %>% ggsave(width=7,height=5)
(qual.rev <- plotQualityProfile(path.rev[qual.samps]) + ggtitle("Qual profiles rev"))
file.path(outDir, "qual.rev.pdf") %>% ggsave(width=7,height=5)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Truncate reads to fixed length 
# Use quality profiles (above) to set the --length parameter for truncation
# Retained length must be long enough to allow overlap of longer amplicons
(truncDir <- file.path(outDir,"trunc")) %>%
  dir.create(showWarnings = F) # Create output directory for truncated reads 
#Loop over forward reads
for(fwd in path.fwd){
  trunc.fwd <- file.path(truncDir, basename(fwd))
  cutadapt.args <- paste(
    "--quiet -g XX", # dummy adapter 
    "--length 225", # Fwd read truncation
    "-o", trunc.fwd, fwd)
  system2("cutadapt", args = cutadapt.args)
}
#Loop over reverse reads
for(rev in path.rev){
  trunc.rev <- file.path(truncDir, basename(rev))
  cutadapt.args <- paste(
    "--quiet -g XX", # dummy adapter 
    "--length 200", # Rev read truncation
    "-o", trunc.rev, rev)
  system2("cutadapt", args = cutadapt.args)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#  
# Quality filter
# Make folder for filtered reads
(filtDir <- file.path(outDir,"filt")) %>%
  dir.create(showWarnings = F) # Create output directory for truncated reads 
# Organize inputs and outputs
trunc.fwd <- file.path(truncDir, basename(path.fwd)) %>% sort
trunc.rev <- file.path(truncDir, basename(path.rev)) %>% sort
filt.fwd <- file.path(filtDir, basename(path.fwd)) %>% sort
filt.rev <- file.path(filtDir, basename(path.rev)) %>% sort
# Run dada2 quality filtering
trim <- filterAndTrim(trunc.fwd, filt.fwd,
                      trunc.rev, filt.rev,
                      maxN=0, maxEE=c(2,2), truncQ=2,
                      multithread=TRUE)
#Update list of filteres file paths to exclude samples with no reads passing filters
filt.fwd <- list.files(filtDir, pattern = "R1.fq.gz",full.names = T)
filt.rev <- list.files(filtDir, pattern = "R2.fq.gz",full.names = T)
  
#Check quality of trimmed and filtered reads
qual.samps <- sample(1:length(filt.fwd),9)
(qual.filt.fwd <- plotQualityProfile(filt.fwd[qual.samps], aggregate=T) + 
    ggtitle("Qual profiles fwd filtered"))
file.path(outDir, "qual.fwd.filtered.pdf") %>% ggsave(width=7,height=5)
(qual.filt.rev <- plotQualityProfile(filt.rev[qual.samps]) + 
    ggtitle("Qual profiles rev filtered"))
file.path(outDir, "qual.rev.filtered.pdf") %>% ggsave(width=7,height=5)
  
#dereplicate
derep.fwd <- derepFastq(filt.fwd, verbose=F)
derep.rev <- derepFastq(filt.rev, verbose=F) 

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#Trim names of derep objects to just sample names 
names(derep.fwd) %<>% gsub(".R1.fq.gz","",.)
names(derep.rev) %<>% gsub(".R2.fq.gz","",.)
  
#Learn errors 
err.fwd <- learnErrors(filt.fwd, multithread=T) 
err.plot.fwd <- plotErrors(err.fwd, nominalQ=T) + ggtitle("Forward reads error model")
file.path(outDir, "errMod.fwd.pdf") %>% ggsave(width=5,height=5)
err.rev <- learnErrors(filt.rev, multithread=T)
err.plot.rev <- plotErrors(err.rev, nominalQ=T) + ggtitle("Reverse reads error model")
file.path(outDir, "errMod.rev.pdf") %>% ggsave(width=5,height=5)
  
#Denoise
dada.fwd <- dada(derep.fwd, err=err.fwd, multithread=TRUE)
dada.rev <- dada(derep.rev, err=err.rev, multithread=TRUE)

#Merge reads
merged <- mergePairs(dada.fwd, derep.fwd, dada.rev, derep.rev, trimOverhang = T)
  
#Make sequence table
seqtab <- merged %>% makeSequenceTable
  
#Remove chimeras
seqtab.nonChimeras <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
  
#Make summary report of what happened
getN <- function(x) sum(getUniques(x))
trim.summary <- trim %>% data.frame %>% rownames_to_column("Sample") 
trim.summary$Sample %<>% strsplit(., ".",fixed=T) %>% sapply(., `[`, 1)
track <- cbind(sapply(dada.fwd, getN), 
               sapply(dada.rev, getN), 
               sapply(merged, getN),
               rowSums(seqtab.nonChimeras)) %>% 
  data.frame %>% rownames_to_column("Sample") 
track$Sample  %<>% strsplit(., ".",fixed=T) %>% sapply(., `[`, 1)
track %<>% left_join(trim.summary,.)
colnames(track) <- c("sample", "input", "filtered", "denoised.fwd", "denoised.rev", "merged", "ChimeraFiltered")
file.path(outDir,"dadaSummary.csv") %>% write.csv(track,.,row.names = F)
  
# Save output
file.path(outDir,"seqTab.rds") %>% saveRDS(seqtab.nonChimeras,.)
