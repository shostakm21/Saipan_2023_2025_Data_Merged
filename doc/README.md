```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("microbiome")
install.packages("vegan")
install.packages("devtools")
BiocManager::install("microbiome")
```

```{r}
#All the packages you will need for entire code, some of these could be unnecessary pending what you want to graph
library(dada2)
library(BiocManager)
library(ggplot2); packageVersion("ggplot2")
library(devtools)
library(vegan)
library(dbplyr)
library(microbiome)
library(tidyverse)
library(conflicted)
library(dplyr)
library(cowplot)

formatPvalues <- function(pvalue) {
  ra<- ""
  if(pvalue <= 0.1) ra<- "."
  if(pvalue <= 0.05) ra<- "*"
  if(pvalue <= 0.01) ra<- "**"
  if(pvalue <= 0.001) ra<- "***"
  return(ra)
}

formatPvalues <- function(pvalue) {
  ra<- ""
  if(pvalue <= 0.1) ra<- "."
  if(pvalue <= 0.05) ra<- "*"
  if(pvalue <= 0.01) ra<- "**"
  if(pvalue <= 0.001) ra<- "***"
  return(ra)
}
```

#Initial Setup: 
```{r setup, include=FALSE}
se <- function(x, ...){sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))}
ci <- function(x, ...){1.96 * sd(x, na.rm = TRUE)}
#Code Dependencies
library(dada2); packageVersion("dada2")
```

Run through the Filter and Trim step as one large group. 
Sort out filtered files into Forward and Reverse reads and subgroups. I did Sample year, since that is how I collected my data and had already pre-sorted data files in that format. You could also do it by year or by site. The sorting helps break it down so that it can run the files without computationally shutting down.
Run the Errors and Mergers on each individual subset and save each subset as an RDS file. In this code, I changed the file locations for each subset. You could also just set it up, so each subset has its own chunk of code and run it that way.
Merge all subset RDS files into a single Sequence Table
Remove chimeras and assign Taxonomy using the single Sequence Table.

```{r}
#Sequence Data Processing - DADA2 Pipeline
#Identify File Path
pathF <- "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/fastq/FWD"
pathR <- "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/fastq/REV"
filtpathF <- file.path(pathF, "filtered")
filtpathR <- file.path(pathR, "filtered")
fastqFs <- sort(list.files(pathF, pattern = "fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern = "fastq.gz"))
if(length(fastqFs) !=length(fastqRs)) stop("Forward and reverse files do not match.")
```

```{r}
FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "CCGYCAATTYMTTTRAGTTT"
trimLeft = c(FWD, REV)
out <- filterAndTrim(fwd = file.path(pathF, fastqFs), filt = file.path(filtpathF, fastqFs), rev = file.path(pathR, fastqRs), filt.rev = file.path(filtpathR, fastqRs), truncLen = c(280,180), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = FALSE, trimLeft = c(19,20))
head(out)
```

This is where the parallelization started. I sorted filtered files into Forward and Reversed and then again by sample cycle. This was done in the Filtered folder on my drive. The Filtered Forward folder had separate folders for SC0 through SC8, the same with the Filtered Reverse folder. In the highlighted steps, I then changed the path locations at the start and end of this section, so where the files were being pulled from and where the final RDS file for that Sample Cycle was stored.

# 2023 Data Set
```{r}
filtpathF <-"/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/fastq/FWD/filtered/FWD_2023_data"
filtpathR <- "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/fastq/REV/filtered/REV_2023_data"
filtFs <- list.files(filtpathF, pattern = "fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern = "fastq.gz", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), '[', 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), '[', 1)

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match,")

names(filtFs) <- sample.names
names(filtRs) <- sample.namesR

set.seed(100)
errF <- learnErrors(filtFs, nbases = 1e8, multithread = TRUE)
errR <- learnErrors(filtRs, nbases = 1e8, multithread = TRUE)

mergers <-  vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err = errF, multithread = TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err = errR, multithread = TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
seqtab23 <- makeSequenceTable(mergers)
saveRDS(seqtab23, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/seqtab23.rds")
```

# 2025 Data Set
```{r}
filtpathF <-"/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/fastq/FWD/filtered/FWD_2025_data"
filtpathR <- "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/fastq/REV/filtered/REV_2025_data"
filtFs <- list.files(filtpathF, pattern = "fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern = "fastq.gz", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), '[', 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), '[', 1)

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match,")

names(filtFs) <- sample.names
names(filtRs) <- sample.namesR

set.seed(100)
errF <- learnErrors(filtFs, nbases = 1e8, multithread = TRUE)
errR <- learnErrors(filtRs, nbases = 1e8, multithread = TRUE)

mergers <-  vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err = errF, multithread = TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err = errR, multithread = TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
seqtab25 <- makeSequenceTable(mergers)
saveRDS(seqtab25, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/seqtab25.rds")
```

# 2023 & 2025 Combined
```{r}
st.all <- mergeSequenceTables(seqtab23, seqtab25)

seqtab <- removeBimeraDenovo(st.all, method = "consensus", multithread = TRUE)

tax <- assignTaxonomy(seqtab, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/fastq/silva_nr99_v138.2_toGenus_trainset.gz", multithread = TRUE)

taxa.print <- tax
rownames(taxa.print) <- NULL
head(taxa.print)

#write.csv(tax, file = "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/ASVs_taxonomy") 
#rite.csv(tax, file = "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/ASVa_taxonomy.csv") 
saveRDS(tax, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/ASVa_taxonomy.rds") 

asv_headers <- vector(dim(seqtab)[2], mode = "character") 
count.asv.tab <- t(seqtab) 
row.names(count.asv.tab) <- sub(">", "", asv_headers) 

#write.csv(count.asv.tab, file = "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/ASVs_counts.csv") 
saveRDS(count.asv.tab, file = "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/ASVs_counts.rds") 

#write.csv(tax, file = "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/ASVs_taxonomy.csv") 
saveRDS(tax, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/ASVs_taxonomy.rds") 

#write.csv(tax, file = "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/asvs_taxonomy.csv") 
saveRDS(tax, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/asvs_taxonomy.rds") 

asv_headers <- vector(dim(seqtab)[2], mode = "character") 
count.asv.tab <- t(seqtab) 
row.names(count.asv.tab) <- sub(">", "", asv_headers) 

#write.csv(count.asv.tab, file = "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/asvs_count.csv") 
saveRDS(count.asv.tab, file = "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/asvs_counts.rds") 

asv_seqs <- colnames(seqtab) 
asv_headers <- vector(dim(seqtab)[2], mode="character") 
for (i in 1:dim(seqtab)[2]) {asv_headers[i] <- paste(">ASV", i, sep="_")} 

asv_fasta <- c(rbind(asv_headers, asv_seqs)) 
asv_otu <- t(seqtab) 
row.names(asv_otu) <- sub(">", "", asv_headers) 
asv_tax <- tax
row.names(asv_tax) <- sub(">", "", asv_headers) 
OTU_TAX_table <- merge(asv_otu, asv_tax, by=0) 

#write(asv_fasta, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/asv_fasta.fa") 
#write.table(asv_otu, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/asv_otu.csv", sep=",", quote=F, col.names=NA) 
#write.table(asv_tax, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/asv_tax.csv", sep=",", quote=F, col.names=NA) 
#write.table(OTU_TAX_table, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/OTU_TAX_table.csv", sep=",", quote=F, col.names=NA)
```

# Merged Data Analysis
```{r}
metadata <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/metadata_2023_2025.csv")
metadata

otu_counts <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/asv_otu.csv") %>%
  pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_counts

taxonomy <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/asv_tax.csv")
taxonomy
```

```{r}
otu_rel_abund <- inner_join(metadata, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund

write.table(otu_rel_abund, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/otu_rel_abund.csv", sep=",", quote=F, col.names=NA)

otu_rel_abund <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/otu_rel_abund.csv")
otu_rel_abund
```

# Stacked Barcharts
```{r}
## Phylum
otu_rel_abund %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, y="Mean Relative Abundance (%)") + theme_classic()
ggsave("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/results/phylum_stacked_barchart.tiff", width=20, height=10)

## Class
otu_rel_abund %>%
  filter(level=="Class") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/results/class_stacked_barchart.tiff", width=25, height=10)

## Order
otu_rel_abund %>%
  filter(level=="Order") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/results/order_stacked_barchart.tiff", width=55, height=10, limitsize = FALSE)

## Family
otu_rel_abund %>%
  filter(level=="Family") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/results/family_stacked_barchart.tiff", width=50, height=10, limitsize = FALSE)
```
