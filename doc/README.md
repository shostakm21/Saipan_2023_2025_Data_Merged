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

otu_count <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/asv_otu.csv") %>%
  pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_count

write.table(otu_count, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/otu_count.csv", sep=",", quote=F, col.names=NA)

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

# ASV Table Manipulation for OTU Count Tables
```{r}
# Have it so ASV is columns and Sample_ID are rows
library(tidyverse)

#df1 = ASV_1 through ASV_23820
df1 <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/df_otu_1.csv")
df1

#df1 = ASV_23822 through ASV_37917
df2 <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/df_otu_2.csv")
df2

#df3 = ASV_37918 through ASV_99999
df3 <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/df_otu_3.csv")
df3

df_otu <- list(df1, df2, df3)
df_otu

write.table(df_otu,"/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/df_otu.csv", sep=",", col.names=NA)
```

```{r}
df_otu <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/df_otu.csv")
df_otu

otu_count <- df_otu %>% select (-c(X)) %>%
  pivot_longer(-sample_id, names_to = "ASV", values_to = "count")

write.table(otu_count, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/otu_count.csv", sep=",", quote=F, col.names=NA)

otu_count <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/otu_count.csv")
otu_count

otu_count %>%
  group_by(sample_id) %>%
  mutate(total = sum(count)) %>%
  filter(total > 5000) %>%
  group_by(ASV) %>%
  mutate(total=sum(count)) %>% 
  filter(total != 0) %>%
  as.data.frame()
#Going to set threshold at 5000
```

# NMDS Plots
```{r}
meta <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/metadata_2023_2025.csv")
meta

df_otu <-read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/df_otu.csv")
df_otu

nmds_asv_otu <- inner_join(meta, df_otu, by="sample_id")
nmds_asv_otu

write.table(nmds_asv_otu, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/nmds_asv_otu.csv", sep=",", quote=F, col.names=NA)
```

```{r}
pc <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/nmds_asv_otu.csv")
pc

#make community matrix: extract columns with ASV information
com <- pc[,9:ncol(pc)]
com

#turn ASV information into a matrix
m_com <- as.matrix(com)
m_com
```

```{r}
#Run NMDS using Bray-Curtis distance
set.seed(123)
nmds <- metaMDS(m_com, distance="bray") #stress = 0.1244532 
nmds
plot(nmds)

#access the specific points data of the NMDS plot & scores
str(nmds)
nmds$points
scores(nmds)

#extract NMDS scores
data.scores = as.data.frame(scores(nmds)$sites)

#add columns to data frame
data.scores$sample_id = pc$sample_id
data.scores$location = pc$location
data.scores$depth = pc$depth
data.scores$sample_type = pc$sample_type
data.scores$sample_year = pc$sample_year

head(data.scores)
```

```{r}
# Samples by Wreck
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xx
ggsave("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/results/NMDS_sample_location.tiff", width = 10, height = 10)
```

```{r}
# Samples by Type
xx1 = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = sample_type))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Sample Type", y = "NMDS2")
xx1
ggsave("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/results/NMDS_sample_type.tiff", width = 10, height = 10)
```

```{r}
# Samples by Year
xx2 = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = sample_year))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "Sample Year", y = "NMDS2")
xx2
ggsave("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/results/NMDS_sample_year.tiff", width = 10, height = 10)
```

# Diversity Testing
```{r}
richness <- function(x){
  sum(x > 0)}

shannon <- function(x){
  relabund <- x[x>0]/sum(x)
  -sum(relabund * log(relabund))
}

simpson <- function(x){
  n <- sum(x)
  sum(x * (x-1) /(n*n-1))
}


richness_biof <- function(x){
  sum(x > 0)}

shannon_biof <- function(x){
  relabund <- x[x>0]/sum(x)
  -sum(relabund * log(relabund))
}

simpson_biof <- function(x){
  n <- sum(x)
  sum(x * (x-1) /(n*n-1))
}
```

```{r}
otu_count <- otu_count %>%
  group_by(sample_id) %>%
  summarize(richness = richness(count),
            shannon = shannon(count), 
            evenness = shannon/log(richness),
            n=sum(count))

write.table(otu_count, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/otu_count_all_diversity_metrics.csv", sep=",", quote=F, col.names=NA)

diversity_metrics <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/otu_count_all_diversity_metrics.csv")
diversity_metrics

write.table(diversity_metrics, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/saipan_all_diversity_metrics.csv", sep=",", quote=F, col.names=NA)
```

```{r}
# Boxplots
diversity_metrics <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/saipan_all_diversity_metrics.csv")
diversity_metrics

diversity_metrics %>%
  group_by(location) %>%
  pivot_longer(cols=c(richness, shannon, evenness), 
               names_to="metric") %>%
ggplot(aes(x=n, y=value, fill= sample_type)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 8, outlier.size = 2) +
  facet_wrap(~metric, nrow=4, scales="free_y")

ggsave("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/results/alpha_diversity_metrics_all_sample_type.tiff", width = 10, height = 20)

diversity_metrics %>%
  group_by(location) %>%
  pivot_longer(cols=c(richness, shannon, evenness), 
               names_to="metric") %>%
ggplot(aes(x=n, y=value, fill= location)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 8, outlier.size = 2) +
  facet_wrap(~metric, nrow=4, scales="free_y")

ggsave("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/results/alpha_diversity_metrics_all_sample_location.tiff", width = 10, height = 20)
#Each point represents a sample, (n) Sum of Count, (X) Total number of sequences for each sample & (Y) Value of diversity metric
```

# Alternative Beta Diversity Method
```{r}
otu_table <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/asv_otu.csv", header=T, row.names=1, check.names=FALSE)
head(otu_table)
```


```{r}
## Transpose the data to have sample names on rows
otu.table.diver <- t(otu_table)
otu.table.diver <- as.data.frame(otu.table.diver)
head(otu.table.diver)

write.table(otu.table.diver,"/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/otu.table.diver.csv", sep=",", quote=F, col.names=NA)
```

```{r}
otu.table.diver <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/otu.table.diver.csv")
otu.table.diver
```

```{r}
data(otu.table.diver)
H <- diversity(otu.table.diver)
H

richness <- specnumber(otu.table.diver)
richness

evenness <- H/log(richness)
evenness
```

```{r}
metadata <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/metadata_2023_2025.csv")
metadata

alpha <- cbind(shannon = H, richness = richness, pielou = evenness, metadata)
write.csv(alpha, "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/diversity_indices_bio_sed_water.csv")
head(alpha)

## Boxplot by Sample Location
plot.shan <- ggplot(alpha, aes(x = location, y = shannon, fill = location)) +
geom_boxplot(size = 0.5, outlier.color = "black", outlier.shape = 8, outlier.size = 2) +
ylab("Shannon's H'") + 
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.shan

plot.rich <-ggplot(alpha, aes(x = location, y = richness, fill = location)) +
geom_boxplot(size = 0.5, outlier.color = "black", outlier.shape = 8, outlier.size = 2) +
ylab("Species Richness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.rich

plot.even <- ggplot(alpha, aes(x = location, y = pielou, fill = location)) +
geom_boxplot(size = 0.5, outlier.color = "black", outlier.shape = 8, outlier.size = 2) +
ylab("Pielou's Evenness") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
plot.even

legend <- get_legend(plot.even)

plot_grid(plot.shan + theme(legend.position = "none"), plot.rich + theme(legend.position = "none"), plot.even + theme(legend.position = "none"),ncol = 3)

ggsave("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/Shannon_Richness_Eveness_all.tiff")
```

## Phyloseq
```{r}
library(phyloseq)
library(Biostrings)

# Read data into R
otu_tab <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/asv_otu.csv")
otu_tab

tax_tab <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/asv_tax.csv")
tax_tab

samples_df <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/metadata_2023_2025.csv")
samples_df
```

```{r}
# Phyloseq objects need to have row.names
otu_tab <- otu_tab %>%
  tibble::column_to_rownames("ASV")

tax_tab <- tax_tab %>%
  tibble::column_to_rownames("ASV")

samples_df <- samples_df %>%
  tibble::column_to_rownames("sample_id")

# Transform OTU & Tax table into matrices
otu_tab <- as.matrix(otu_tab)
tax_tab <- as.matrix(tax_tab)

#Transform into Phyloseq Objects
ASV = otu_table(otu_tab, taxa_are_rows = TRUE)
TAX = tax_table(tax_tab)
sample_id = sample_data(samples_df)
  
ps <- phyloseq(ASV, TAX, sample_id)
ps

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 43720 taxa and 153 samples ]
#sample_data() Sample Data:       [ 153 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 43720 taxa by 6 taxonomic ranks ]

# Visualize Data
sample_names(ps)
rank_names(ps)
sample_variables(ps)

# Normalize number of reads in each sample using median sequencing depth
total = median(sample_sums(ps))
standf = function(x, t=total) round(t * (x / sum(x)))
ps = transform_sample_counts(ps, standf)

# Bar graphs based on division
plot_bar(ps, fill = "Phylum")

plot_bar(ps, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

# Heatmaps
plot_heatmap(ps, method = "NMDS", distance = "bray")

ps_abund <- filter_taxa(ps, function(x) sum(x > total*0.20) > 0, TRUE)
ps_abund

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 17 taxa and 153 samples ]
#sample_data() Sample Data:       [ 153 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 17 taxa by 6 taxonomic ranks ]
#OTU Table:          [8 taxa and 5 samples]

otu_table(ps)[1:8, 1:5]

#plot_heatmap(ps_abund, method = "NMDS", distance = "bray")
```

# Simper Analysis
```{r}
simper <- simper(otu.table.diver, metadata_5000$location, permutations=999)
options(max.print=500)
summary(simper)
dput(simper, file = "/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/simp_location.txt")
sim <- dget("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/simp_location.txt")
summary(sim)
```

# Rank Abundance Curves
```{r}
#
library(BiodiversityR)
BiodiversityRGUI()

library(vegan)
library(ggplot2)
library(ggrepel)

saipan <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/saipan_abund_count.csv")
saipan

str(saipan)

saipan.env <- read.csv("/Users/maggieshostak/Desktop/Saipan_R_Studio/post_rarefaction/saipan.env.csv")
saipan.env

saipan.env$sample_type<-as.factor(saipan.env$sample_type)
saipan.env$location<-as.factor(saipan.env$location)
str(saipan.env)

saipan.env

RankAbun.1 <- rankabundance(saipan)
RankAbun.1

rankabunplot(RankAbun.1, scale='abundance', addit=FALSE, specnames=c(1,2,3,4,5))
rankabunplot(RankAbun.1, scale='proportion', addit=FALSE, specnames=c(1,2,3,4,5))
rankabunplot(RankAbun.1, scale='logabun', addit=FALSE, specnames=c(1:30),srt=45, ylim=c(1,100))

rankabuncomp(saipan, y=saipan.env, factor='location', scale='proportion', legend=FALSE)
```

```{r}
otu_counts_5 <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025/asv_otu_top_5_2023.csv") %>%
  pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_counts_5

Tax_5 <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025/asv_tax_top_5_2023.csv")
Tax_5

data5 <- otu_counts_5 %>%
  left_join(Tax_5, by="ASV")
data5

meta <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025/metadata_saipan.csv")
meta

data5 <- data5 %>%
  left_join(meta, by="sample_id")
data5

write.table(data5, "/Users/maggieshostak/Desktop/Saipan_2023_2025/asv_count_tax_top_5_metadata_2023.csv", sep=",", quote=F, col.names=NA)

data5 %>%
  ggplot(aes(x = sample_id, y = count)) +
  facet_grid(~ location, scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")

ggsave("/Users/maggieshostak/Desktop/Saipan_2023_2025/master_data_table_top_5_plot.tiff", width = 40, height = 20)

data5 %>%
  ggplot(aes(x = sample_id, y = count)) +
  facet_grid(~ location + metal_type, scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")

ggsave("/Users/maggieshostak/Desktop/Saipan_2023_2025/master_data_table_top_5_plot2.tiff", width = 40, height = 20)

data5 %>%
  ggplot(aes(x = sample_id, y = count)) +
  facet_grid(~ location, scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill", width = 1) +
  scale_y_continuous(name = "Relative Abundance",
                     labels = scales::percent) +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle= 30))#,
       # axis.text.y = element_text(colour = "black"),
       # strip.text = element_text(face = "bold"),
       # strip.background = element_blank ())

ggsave("/Users/maggieshostak/Desktop/Saipan_2023_2025/master_data_table_top_5_plot3.tiff", width = 30, height = 20)

data5 %>%
  ggplot(aes(x = sample_id, y = count)) +
  facet_grid(~ location + metal_type, scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill", width = 1) +
  scale_y_continuous(name = "Relative Abundance",
                     labels = scales::percent) +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle= 30))#,
        #axis.text.y = element_text(colour = "black"),
        #strip.text = element_text(face = "bold"),
        #strip.background = element_blank ())

ggsave("/Users/maggieshostak/Desktop/Saipan_2023_2025/master_data_table_top_5_plot4.tiff", width = 45, height = 20)

data5 %>%
  ggplot(aes(x = sample_id, y = count)) +
  facet_grid(~ location + depth, scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill", width = 1) +
  scale_y_continuous(name = "Relative Abundance",
                     labels = scales::percent) +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle= 30))#,
        #axis.text.y = element_text(colour = "black"),
        #strip.text = element_text(face = "bold"),
        #strip.background = element_blank ())

ggsave("/Users/maggieshostak/Desktop/Saipan_2023_2025/master_data_table_top_5_plot5.tiff", width = 35, height = 20)
```

# ANOSIM
```{r}
pc_ano <- read.csv("/Users/maggieshostak/Desktop/Saipan_2023_2025_R_Studio/Merged_2023_2025_Data/data/nmds_asv_otu.csv")
pc_ano
```

```{r}
# All Samples: Biof vs Sed vs Water by Location
com_ano = pc_ano[,9:ncol(pc_ano)]
m_com_ano = as.matrix(com_ano)
ano_all = anosim(m_com_ano, pc_ano$location, distance = "bray", permutations = 9999)
ano_all

# ANOSIM statistic R: 0.1357
      # Significance: 1e-04
```

```{r}
# All Samples: Biof vs Sed vs Water by Sample Type
com_ano1 = pc_ano[,9:ncol(pc_ano)]
m_com_ano1 = as.matrix(com_ano1)
ano_all1 = anosim(m_com_ano1, pc_ano$sample_type, distance = "bray", permutations = 9999)
ano_all1

# ANOSIM statistic R: 0.5486
      # Significance: 1e-04
```

```{r}
# All Samples: Biof vs Sed vs Water by Depth
com_ano2 = pc_ano[,9:ncol(pc_ano)]
m_com_ano2 = as.matrix(com_ano2)
ano_all2 = anosim(m_com_ano2, pc_ano$depth, distance = "bray", permutations = 9999)
ano_all2

# ANOSIM statistic R: 0.1488
      # Significance: 1e-04
```

```{r}
# All Samples: Biof vs Sed vs Water by Year
com_ano2 = pc_ano[,9:ncol(pc_ano)]
m_com_ano2 = as.matrix(com_ano2)
ano_all2 = anosim(m_com_ano2, pc_ano$sample_year, distance = "bray", permutations = 9999)
ano_all2

# ANOSIM statistic R: 0.744
      # Significance: 1e-04

```
