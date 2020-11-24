#' ---
#' title: Example analysis script
#' author: Devin R Leopold
#' date: Nov. 23, 2020
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---

#+ warning=F, message=F
#Load packages
library(tidyverse)
library(magrittr)
library(phyloseq)
library(vegan)

sessionInfo()

#' ## Prepare data for analysis

#' ### Read in phyloseq object
(phy <- readRDS("output/compiled/phy.rds"))

#' ### Subset samples to remove controls
phy %<>% subset_samples(Samp_site %in% c("Corvallis_blk2","Westport")) %>% 
  prune_taxa(taxa_sums(.) > 0,.)

#' ### Look at sequencing depth
minDepth <- 250
data.frame(SeqDepth=sort(sample_sums(phy)), Site=sample_data(phy)$Samp_site) %>%
  mutate(cutoff=SeqDepth>minDepth) %>%
  ggplot(aes(x=Site, y=SeqDepth)) +
    geom_violin() +
    geom_point(aes(color=cutoff),position=position_jitter(width=0.1)) +
    theme_classic()
  
#' ### Remove samples below sequencing depth cutoff
(phy %<>% prune_samples(sample_sums(.)>minDepth,.))

#' ### Cluster similar sequence variants
cluster <- function(phy.in,method="single",dissimilarity=0.01){
  require(DECIPHER)
  clusters <- DistanceMatrix(refseq(phy.in), includeTerminalGaps = T, processors=NULL) %>%
    IdClusters(method=method, cutoff=dissimilarity, processors=NULL) %>%
    rownames_to_column("OTUid")
  for(i in unique(clusters$cluster)){
    foo <- clusters$OTUid[which(clusters$cluster==i)] 
    if(length(foo)>1){phy.in %<>% merge_taxa(foo)}
  }
  return(phy.in)
}
(phy %<>% cluster("single"))

#' ### Remove low abundance OTUs by prevelance
# Only keep OTUs present in at least 1% of samples
(phy %>% filter_taxa(., function(x) {sum(x>0) > 0.01*nsamples(.)}, TRUE))

#' ### Convert to proportional abundance
phy %<>% transform_sample_counts(function(x){x*min(sample_sums(.)/sum(x))})

#' ## Community analyses

#' ### Plot an ordination (PCoA)
phy.ord <- phy %>% ordinate("MDS","bray")
plot_ordination(phy,phy.ord,color="Samp_site") + theme_classic()
#Save image to file
ggsave("output/figs/PCoA.pdf")

#' ### Run a perMANOVA
phy.dist <- phy %>% phyloseq::distance("bray")  
adonis(phy.dist~sample_data(phy)$Samp_site)

  
