#Compile MiSeq library and sample data into single phyloseq object 
#Also add taxonomy and remove non-target sequences
 
library(tidyverse)
library(magrittr)
library(foreach)
library(dada2)
library(Biostrings)
library(phyloseq)

# Create scratch directory
dir.create("output/scratch", showWarnings = F)

#read in sample meta data
meta <- read.csv("data/SampleData.csv",as.is=T,row.names=1) #row.names = column with sample names 

#read in denoised sequence table
seqTab <- readRDS("output/processing/dada/seqTab.rds")

#extract sequences from OTU table 
seqs <- getSequences(seqTab) %>% DNAStringSet

#rename OTUs with MD5 hash
colnames(seqTab) %<>% as.character %>% openssl:::md5()
names(seqs) <- colnames(seqTab)

# Write sequences to fasta
writeXStringSet(seqs,file="output/scratch/seqs.fasta",width=500)

#### Filter non-fungal sequences #### 
# by matching against the UNITE all eukaryotes database 

# Download and extract UNITE database files (skip if already done)
if(!file.exists("output/UNITE/UNITE.8.0.Alleuk.fasta.gz")){
  download.file("https://files.plutof.ut.ee/public/orig/1D/B9/1DB95C8AC0A80108BECAF1162D761A8D379AF43E2A4295A3EF353DD1632B645B.gz",
                "output/scratch/tmpDb.gz",
                method="wget", quiet=T)
  # Extract
  untar("output/scratch/tmpDb.gz", exdir="output/scratch")
  # Remove unidentified taxa
  eukDB <- readDNAStringSet("output/scratch/sh_general_release_s_all_04.02.2020/sh_general_release_dynamic_s_all_04.02.2020.fasta")
  eukDB %<>% .[!grepl("^unidentified",names(.))]
  # Remove problematic character from one species id (causes problems for vseach)
  names(eukDB) %<>% gsub("×","",.)
  # Write database to file
  dir.create("output/UNITE/", showWarnings = F)
  writeXStringSet(eukDB, "output/UNITE/UNITE.8.0.Alleuk.fasta.gz", compress=F)
  # Delete temporary files
  unlink("output/scratch/tmpDb.gz")
  unlink("output/scratch/sh_general_release_s_all_04.02.2020/", recursive = T)
  rm(eukDB)
}

# Run search with vsearch returning up to 3 matches
vsearch.flags <- paste("--usearch_global output/scratch/seqs.fasta",
                       "--db output/UNITE/UNITE.8.0.Alleuk.fasta.gz",
                       "--id 0.65", "--strand both",
                       "--userout output/scratch/UNITEmatches.txt",
                       "--userfields query+target+id",
                       "--notmatched output/scratch/noMatch.fasta",
                       "--maxhits 3",
                       "--maxaccepts 500 --maxrejects 0"
)
system2("vsearch", args = vsearch.flags)
# Read in search results and calculate the proportion of top matches that are fungal
UNITEmatches <- read.table("output/scratch/UNITEmatches.txt",as.is=T) %>%
  mutate(isFungi=grepl("k__Fungi",V2)) %>%
  group_by(V1) %>%
  summarise(isFungi=sum(isFungi)/n(),
            tmp=V2[1])
# Make list of OTUs to remove from the data (more than half of the top matches are non-fungal)
fungal <- UNITEmatches$V1[UNITEmatches$isFungi>0.5]

#### Extract ITS region with ITSx ####
# Not possible with all primer combinations or with other (non-ITS) markers

# Write updated (filtered) fasta
writeXStringSet(seqs[fungal],file="output/scratch/seqs.fasta",width=500)

# Run ITSx
ITSx.args <- paste("-i output/scratch/seqs.fasta",
                   "-t 'fungi'",
                   "--preserve T",
                   "--complement T", 
                   "--summary T",
                   "-o output/scratch/ITSx",
                   "--only_full T",
                   "--cpu 2",
                   "-E 1e-2")
system2("ITSx", args = ITSx.args)
  
# Read in new sequences which are trimmed to only ITS2 and are now in the correct orientation
seqs <- readDNAStringSet("output/scratch/ITSx.ITS2.fasta")

#### Create phyloseq object ####
phy <- phyloseq(otu_table(seqTab,taxa_are_rows = F), 
                sample_data(meta),
                refseq(readDNAStringSet("output/scratch/ITSx.ITS2.fasta")))

#### Cluster ####
# Merge OTU with identical (or in this case 99% similar) sequences
cluster <- function(phy.in,method="single",dissimilarity=0){
  require(DECIPHER)
  clusters <- DistanceMatrix(refseq(phy.in), includeTerminalGaps = T, processors=NULL) %>%
    IdClusters(method=method, cutoff=dissimilarity, processors=NULL) 
  clusters[,1] %<>% as.character()
  tax_table(phy.in) <- clusters %>% as.matrix %>% tax_table
  phy.in %<>% speedyseq::tax_glom("cluster")
  tax_table(phy.in) <- NULL
  return(phy.in)
}  
phy %<>% cluster(dissimilarity = 0.01) 

#### Assign taxonomy #### 

# Next get curated database of only fungi for taxonomic prediction
if(!file.exists("output/UNITE/UNITE.8.2.fungi.fasta.gz")){
  download.file("https://files.plutof.ut.ee/public/orig/E7/28/E728E2CAB797C90A01CD271118F574B8B7D0DAEAB7E81193EB89A2AC769A0896.gz",
                "output/scratch/tmpDb.gz",
                method="wget", quiet=T)
  # Extract
  untar("output/scratch/tmpDb.gz", exdir="output/scratch/")
  # Remove taxa only assigned to kingdom Fungi
  fungiDB <- readDNAStringSet("output/scratch/sh_general_release_04.02.2020/sh_general_release_dynamic_04.02.2020.fasta")
  fungiDB %<>% .[!grepl("^Fungi",names(.))]
  # Write database to file
  dir.create("output/UNITE/", showWarnings = F)
  writeXStringSet(fungiDB, "output/UNITE/UNITE.8.2.fungi.fasta.gz", compress=T)
  # Delete temporary files
  unlink("output/scratch/tmpDb.gz")
  unlink("output/scratch/sh_general_release_04.02.2020", recursive=T)
  rm(fungiDB)
}

# Assign taxonomy with DADA2 implementation of RDP classifier
tax <- assignTaxonomy(refseq(phy), "output/UNITE/UNITE.8.2.fungi.fasta.gz", multithread=T)
rownames(tax) <- taxa_names(phy)

# Add taxonomy to phyloseq object
phy %<>% merge_phyloseq(tax_table(tax))

###############
### OUTPUTs ###
###############

# Set final OTU names (in order of relative abundance)
OTU.names <- paste0("OTU.",1:ntaxa(phy)) 
names(OTU.names) <- taxa_sums(phy) %>% sort(decreasing = T) %>% names
taxa_names(phy) %<>% OTU.names[.] %>% unname

# Save phyloseq object
saveRDS(phy,"output/processing/compiled/phy.rds")

# Save components for possible manual inspection
otu_table(phy) %>% write.csv("output/processing/compiled/OTU.table.csv")
tax_table(phy) %>% write.csv("output/processing/compiled/taxonomy.table.csv")

# Remove scratch folder
unlink("output/scratch", T)