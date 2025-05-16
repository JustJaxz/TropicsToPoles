#------------------------------------------------------------------------------------------------------------------------------#

# Paper title: Tropics to the poles: Diversity and composition of coastal eukaryotic marine microalgae communities across five ecoregions.
# Paper authors: Jacqui Stuart (1,2), Ken Ryan (1), Natalie Robinson (4), Svenja Halfter (4), Craig Stewart (4), John K. Pearman (1), 
# Jacob Thomson-Laing (2), Kirsty F. Smith (2,3)

# Affiliations:
# 1.	School of Biological Sciences, Victoria University of Wellington, PO Box 600, Wellington 6140, New Zealand. 
# 2.	Cawthron Institute, Private Bag 2, Nelson 7042, New Zealand. 
# 3.	School of Biological Sciences, University of Auckland, Private Bag 92019, Auckland 1142, New Zealand.
# 4.	National Institute of Water and Atmospheric Research (NIWA), Private Bag 14901, Kilbirnie, Wellington 6241

## Assembly processes analysis
# This R code conducts an analysis on eukaryotic microalgae communities in different ecoregions. 
# The code performs several tasks, including data preparation, clustering of sequences at 97% and 99% similarity, 
# creating and analyzing phylogenetic trees, calculating beta diversity metrics (weighted beta MNTD and Raup-Crick indices), 
# and running statistical analyses to assess dispersal limitation, homogenizing dispersal, variable selection, 
# and homogeneous selection among the sampled sites. The results are then summarized within and between ecoregions, 
# providing insights into the ecological processes shaping microalgae communities in the studied areas.

# The final output is a map of the study area displaying assembly processes between sites and ecoregions. 

#------------------------------------------------------------------------------------------------------------------------------#

## Load Libraries
#install.packages("marmap")
library(dplyr)
library(Biostrings)
library(msa)
library(bios2mds)
library(phangorn)
library(picante)
library(iCAMP)
library(tidyr)
library(ggplot2)
library(marmap)

#------------------------------------------------------------------------------------------------------------------------------#

## Data prep
# NOTE: FIRST RUN DATA SETUP SCRIPT IN RMD FILE 

# Set the working directory to the specified path
setwd("C:/Documents/R")

# Load rarefied alpha diversity data into 'ps' object
ps <- alpha.waterSamp_pru.rare

# Extract ASV (amplicon sequence variant) names from the dataset
asv_seqs <- taxa_names(ps)

# Create headers for ASV sequences by appending ">ASV" and a numerical index
asv_headers <- vector(dim(data.frame(taxa_names(ps)))[1], mode="character")
for (i in 1:dim(data.frame(taxa_names(ps)))[1]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# Combine headers and sequences into a FASTA format and write it to a file
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "pole.tropics.ASVs.fa")

# Create a new set of headers without the ">" symbol for taxonomy table
asv_headers1 <- vector(dim(data.frame(taxa_names(ps)))[1], mode="character")
for (i in 1:dim(data.frame(taxa_names(ps)))[1]) {
  asv_headers1[i] <- paste("ASV", i, sep="_")
}

# Replace original ASV headers with the new ones in the 'ps' object
taxa_names(ps) <- asv_headers1

# Create a data frame from the taxonomy table and add ASV as a column
taxa.df <- data.frame(tax_table(ps))
taxa.df$ASV <- rownames(taxa.df)

#------------------------------------------------------------------------------------------------------------------------------#

## Clustering at 99%

# Run on grunty computer!
#~/vsearch-2.25.0-macos-aarch64/bin/vsearch --cluster_fast pole.tropics.ASVs.fa --threads 8 --id 0.99 --strand plus --sizeout --centroids centroids.99.fasta --uc pole.tropics.uc.99.txt

uc.99 <- read.table("pole.tropics.uc.99.txt")
uc.99 <- uc.99%>% filter(V1 != "C") %>% select(V9, V10)
uc.99 <- uc.99 %>% mutate(cluster.99 = factor(ifelse(V10 == "*", V9, V10))) %>% select(-V10)


## Clustering at 97%
# Run on grunty computer!
#~/vsearch-2.25.0-macos-aarch64/bin/vsearch --cluster_fast pole.tropics.ASVs.fa --threads 8 --id 0.97 --strand plus --sizeout --centroids centroids.97.fasta --uc pole.tropics.uc.97.txt

uc.97 <- read.table("pole.tropics.uc.97.txt")
uc.97 <- uc.97 %>% filter(V1 != "C") %>% select(V9, V10)
uc.97 <- uc.97 %>% mutate(cluster.97 = factor(ifelse(V10 == "*", V9, V10))) %>% select(-V10)


# This combines 97 and 99 and taxonomy table into one dataframe
taxa.df <- left_join(taxa.df, uc.97, by = c("ASV" = "V9"))
taxa.df <- left_join(taxa.df, uc.99, by = c("ASV" = "V9"))


rownames(taxa.df) <- taxa.df$ASV
taxa.df <- as.matrix(taxa.df)

#reads it back into phyloseq object
tax_table(ps) <- taxa.df

# merging ASVs that belong to a single OTU
ps.99 <- tax_glom(ps, "cluster.99")
ps.97  <- tax_glom(ps, "cluster.97")

sample_data(ps)


#------------------------------------------------------------------------------------------------------------------------------#
## Calculating pairwise genetic distances and making Tree

file.pole.tropics <- readDNAStringSet("pole.tropics.ASVs.fa") # Read the DNA sequences from the FASTA file into a DNAStringSet object
align.pole.tropics <- msa(file.pole.tropics, method="ClustalOmega", verbose=T) # Perform multiple sequence alignment using ClustalOmega method
cv.align.pole.tropics <- msaConvert(align.pole.tropics, type=c("bios2mds::align")) # Convert the alignment to bios2mds format
export.fasta(cv.align.pole.tropics, outfile = "pole.tropics_aligned.fas", ncol(align.pole.tropics), open = "w") # Export the aligned sequences to a new FASTA file

pole.tropics.aligned <- read.dna("pole.tropics_aligned.fas", format="fasta") # Read the aligned sequences back into a DNA object
dna_dist <- dist.ml(pole.tropics.aligned, model="JC69") # Calculate the pairwise genetic distances using the Jukes-Cantor model
pole.tropics.aligned_upgma  <- upgma(dna_dist) # Create an UPGMA (Unweighted Pair Group Method with Arithmetic Mean) tree based on the distances
plot(pole.tropics.aligned_upgma) # Plot the UPGMA tree

write.tree(pole.tropics.aligned_upgma, file="pole.tropics_upgma.tre") # Write the UPGMA tree to a file in Newick format

# NOTE: Convert this file to '.nexus' file before continuing


#------------------------------------------------------------------------------------------------------------------------------#
## ASVs
ps <- merge_samples(ps, "siteID")
sample_data(ps)

write.csv(otu_table(ps), "pole.tropics_comm.csv")
write.csv(data.frame(sample_data(ps)), "pole.tropics_comm.sd.csv")

tree.pole.tropics <- read.nexus(file="poles.tropics.nexus")
pole.tropics.comm <-read.csv("pole.tropics_comm.csv", h=T, row.names=1)

tree.pole.tropics <- read.nexus(file="poles.tropics.nexus")
pole.tropics.comm <-read.csv("pole.tropics_comm.csv", h=T, row.names=1)
comm.pole.tropics_log <- log(pole.tropics.comm + 1)

comm.pole.tropics_log.p <- match.phylo.comm(tree.pole.tropics, comm.pole.tropics_log)$comm
tree.pole.tropics_log.p <- match.phylo.comm(tree.pole.tropics, comm.pole.tropics_log)$phy

# write.csv(comm.ARMS_log.p, "comm.ARMS_log.p.csv")

otu <- comm.pole.tropics_log.p; ##set your OTU table to the object 'otu' if you have it imported already##
dim(otu); # this gives the dimensions
otu[1:5,1:5]; # this gives a look at the first 5 rows and columns


## read in the phylogeny
phylo <- tree.pole.tropics_log.p
phylo; # a summary of the phylogeny
#plot.phylo(phylo,typ="fan"); # a quick plot


## make sure the names on the phylogeny are ordered the same as the names in otu table
match.phylo.otu = match.phylo.comm(phy =phylo, comm = otu);
str(match.phylo.otu);

## calculate empirical betaMNTD
comm=match.phylo.otu$comm

nworker=8 # parallel computing thread number
rand.time=999 # usually use 1000 for real data.
save.wd="big/pd.big.pole.tropics2" # please change to the folder you want to save the pd.big output.
pd.big=pdist.big(tree = match.phylo.otu$phy, wd=save.wd, nworker = nworker)
bNTI=bNTI.big(comm=comm, pd.desc=pd.big$pd.file,
              pd.spname=pd.big$tip.label,pd.wd=pd.big$pd.wd,
              spname.check=TRUE, nworker=nworker, memo.size.GB=120,
              weighted=TRUE, exclude.consp=FALSE,rand=rand.time,
              output.dtail=FALSE, RC=FALSE, trace=TRUE)

bNTI.mat <- as.matrix(bNTI)
bNTI.mat[!lower.tri(bNTI.mat)] <- NA

weighted.bNTI.melt <- reshape2::melt(bNTI.mat)

write.csv(weighted.bNTI.melt,"weighted_bNTI.melt_pole.tropics.csv",quote=F, row.names = FALSE)


tree.pole.tropics <- read.nexus(file="poles.tropics.nexus")
pole.tropics.comm <-read.csv("pole.tropics_comm.csv", h=T, row.names=1)
comm.pole.tropics_log <- log(pole.tropics.comm + 1)
comm.pole.tropics_log.p <- match.phylo.comm(tree.pole.tropics, comm.pole.tropics_log)$comm
tree.pole.tropics_log.p <- match.phylo.comm(tree.pole.tropics, comm.pole.tropics_log)$phy

rand.time=999 # usually use 1000 for real data.
nworker=7 # parallel computing thread number

set.seed(20)
RaupCrick_results=RC.pc(comm.pole.tropics_log.p, rand = rand.time,
                        nworker = nworker, weighted = TRUE,
                        sig.index="RC")

RaupCrick_results.mat <- as.matrix(RaupCrick_results$index)
RaupCrick_results.mat[!lower.tri(RaupCrick_results.mat)] <- NA

RaupCrick_results.melt <- reshape2::melt(RaupCrick_results.mat)

#RaupCrick_results.melt <- RaupCrick_results.melt %>% select(-X)

write.csv(RaupCrick_results.mat, "RaupCrick_results_pole.tropics.csv")
write.csv(RaupCrick_results.melt, "RaupCrick_results_pole.tropics_melt.csv", row.names = FALSE)


#------------------------------------------------------------------------------------------------------------------------------#
## Run comparison analysis for 97%

ps.97 <- merge_samples(ps.97, "siteID")
sampleData(ps.97)
write.csv(otu_table(ps.97), "pole.tropics.97_comm.csv")
#write.csv(sample_data(ps.97), "pole.tropics.97_comm.sd.csv")

tree.pole.tropics.97 <- read.nexus(file="poles.tropics.nexus")
pole.tropics.97.comm <-read.csv("pole.tropics.97_comm.csv", h=T, row.names=1)

tree.pole.tropics.97 <- read.nexus(file="poles.tropics.nexus")
pole.tropics.97.comm <-read.csv("pole.tropics.97_comm.csv", h=T, row.names=1)
comm.pole.tropics.97_log <- log(pole.tropics.97.comm + 1)

comm.pole.tropics.97_log.p <- match.phylo.comm(tree.pole.tropics.97, comm.pole.tropics.97_log)$comm
tree.pole.tropics.97_log.p <- match.phylo.comm(tree.pole.tropics.97, comm.pole.tropics.97_log)$phy

#write.csv(comm.ARMS_log.p, "comm.ARMS_log.p.csv")
otu <- comm.pole.tropics.97_log.p; ##set your OTU table to the object 'otu' if you have it imported already##
dim(otu); # this gives the dimensions
otu[1:5,1:5]; # this gives a look at the first 5 rows and columns

## read in the phylogeny
phylo <- tree.pole.tropics.97_log.p
phylo; # a summary of the phylogeny

## make sure the names on the phylogeny are ordered the same as the names in otu table
match.phylo.otu = match.phylo.comm(phy =phylo, comm = otu);
str(match.phylo.otu);

## calculate empirical betaMNTD
comm=match.phylo.otu$comm

nworker=8 # parallel computing thread number
rand.time=999 # usually use 1000 for real data.

save.wd="big/pd.big.pole.tropics.97" # please change to the folder you want to save the pd.big output.
pd.big=pdist.big(tree = match.phylo.otu$phy, wd=save.wd, nworker = nworker)
bNTI=bNTI.big(comm=comm, pd.desc=pd.big$pd.file,
              pd.spname=pd.big$tip.label,pd.wd=pd.big$pd.wd,
              spname.check=TRUE, nworker=nworker, memo.size.GB=120,
              weighted=TRUE, exclude.consp=FALSE,rand=rand.time,
              output.dtail=FALSE, RC=FALSE, trace=TRUE)

bNTI.mat <- as.matrix(bNTI)
bNTI.mat[!lower.tri(bNTI.mat)] <- NA

weighted.bNTI.melt <- reshape2::melt(bNTI.mat)

write.csv(weighted.bNTI.melt,"weighted_bNTI.melt_pole.tropics-diatom.97.csv",quote=F, row.names = FALSE)


tree.pole.tropics.97 <- read.nexus(file="poles.tropics.nexus")
pole.tropics.97.comm <-read.csv("pole.tropics.97_comm.csv", h=T, row.names=1)
comm.pole.tropics.97_log <- log(pole.tropics.97.comm + 1)
comm.pole.tropics.97_log.p <- match.phylo.comm(tree.pole.tropics.97, comm.pole.tropics.97_log)$comm
tree.pole.tropics.97_log.p <- match.phylo.comm(tree.pole.tropics.97, comm.pole.tropics.97_log)$phy

rand.time=999 # usually use 1000 for real data.
nworker=7 # parallel computing thread number

set.seed(20)
RaupCrick_results=RC.pc(comm.pole.tropics.97_log.p, rand = rand.time,
                        nworker = nworker, weighted = TRUE,
                        sig.index="RC")

RaupCrick_results.mat <- as.matrix(RaupCrick_results$index)
RaupCrick_results.mat[!lower.tri(RaupCrick_results.mat)] <- NA

RaupCrick_results.melt <- reshape2::melt(RaupCrick_results.mat)

#RaupCrick_results.melt <- RaupCrick_results.melt %>% select(-X)
write.csv(RaupCrick_results.mat, "RaupCrick_results_pole.tropics.dino.97.csv")
write.csv(RaupCrick_results.melt, "RaupCrick_results_pole.tropics.dino.97_melt.csv", row.names = FALSE)


#------------------------------------------------------------------------------------------------------------------------------#

## Run comparison analysis for 99% cluster
ps.99 <- merge_samples(ps.99, "siteID")
write.csv(otu_table(ps.99), "pole.tropics.99_comm.csv")
#write.csv(sample_data(ps.99), "pole.tropics.99_comm.sd.csv")

tree.pole.tropics.99 <- read.nexus(file="poles.tropics.nexus")
pole.tropics.99.comm <-read.csv("pole.tropics.99_comm.csv", h=T, row.names=1)
tree.pole.tropics.99 <- read.nexus(file="poles.tropics.nexus")
pole.tropics.99.comm <-read.csv("pole.tropics.99_comm.csv", h=T, row.names=1)
comm.pole.tropics.99_log <- log(pole.tropics.99.comm + 1)


comm.pole.tropics.99_log.p <- match.phylo.comm(tree.pole.tropics.99, comm.pole.tropics.99_log)$comm
tree.pole.tropics.99_log.p <- match.phylo.comm(tree.pole.tropics.99, comm.pole.tropics.99_log)$phy

#write.csv(comm.ARMS_log.p, "comm.ARMS_log.p.csv")
otu <- comm.pole.tropics.99_log.p; ##set your OTU table to the object 'otu' if you have it imported already##
dim(otu); # this gives the dimensions
otu[1:5,1:5]; # this gives a look at the first 5 rows and columns

## read in the phylogeny

phylo <- tree.pole.tropics.99_log.p
phylo; # a summary of the phylogeny
#plot.phylo(phylo,typ="fan"); # a quick plot

## make sure the names on the phylogeny are ordered the same as the names in otu table

match.phylo.otu = match.phylo.comm(phy =phylo, comm = otu);
str(match.phylo.otu);

## calculate empirical betaMNTD
comm=match.phylo.otu$comm

nworker=8 # parallel computing thread number
rand.time=999 # usually use 1000 for real data.

save.wd="big/pd.big.pole.tropics.99" # please change to the folder you want to save the pd.big output.
pd.big=pdist.big(tree = match.phylo.otu$phy, wd=save.wd, nworker = nworker)
bNTI=bNTI.big(comm=comm, pd.desc=pd.big$pd.file,
              pd.spname=pd.big$tip.label,pd.wd=pd.big$pd.wd,
              spname.check=TRUE, nworker=nworker, memo.size.GB=120,
              weighted=TRUE, exclude.consp=FALSE,rand=rand.time,
              output.dtail=FALSE, RC=FALSE, trace=TRUE)

bNTI.mat <- as.matrix(bNTI)
bNTI.mat[!lower.tri(bNTI.mat)] <- NA

weighted.bNTI.melt <- reshape2::melt(bNTI.mat)

write.csv(weighted.bNTI.melt,"weighted_bNTI.melt_pole.tropics.99.csv",quote=F, row.names = FALSE)



tree.pole.tropics.99 <- read.nexus(file="poles.tropics.nexus")

pole.tropics.99.comm <-read.csv("pole.tropics.99_comm.csv", h=T, row.names=1)

comm.pole.tropics.99_log <- log(pole.tropics.99.comm + 1)

comm.pole.tropics.99_log.p <- match.phylo.comm(tree.pole.tropics.99, comm.pole.tropics.99_log)$comm
tree.pole.tropics.99_log.p <- match.phylo.comm(tree.pole.tropics.99, comm.pole.tropics.99_log)$phy

rand.time=999 # usually use 1000 for real data.
nworker=7 # parallel computing thread number

set.seed(20)
RaupCrick_results=RC.pc(comm.pole.tropics.99_log.p, rand = rand.time,
                        nworker = nworker, weighted = TRUE,
                        sig.index="RC")

RaupCrick_results.mat <- as.matrix(RaupCrick_results$index)
RaupCrick_results.mat[!lower.tri(RaupCrick_results.mat)] <- NA

RaupCrick_results.melt <- reshape2::melt(RaupCrick_results.mat)

#RaupCrick_results.melt <- RaupCrick_results.melt %>% select(-X)

write.csv(RaupCrick_results.mat, "RaupCrick_results_pole.tropics.99.csv")
write.csv(RaupCrick_results.melt, "RaupCrick_results_pole.tropics.99_melt.csv", row.names = FALSE)

#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#

## Analysis
### Stegen Analysis

RaupCrick_results.melt.transect <- read.csv("RaupCrick_results_pole.tropics.dino.97_melt.csv")
BNTI_results.melt.transect <- read.csv("weighted_bNTI.melt_pole.tropics-dino.97.csv")

RaupCrick_results.melt.transect1 <- RaupCrick_results.melt.transect %>% dplyr::filter(!is.na(value)) %>% unite("Comparison", Var1:Var2, remove = TRUE)
colnames(RaupCrick_results.melt.transect1)[2] <- "RC"

BNTI_results.melt.transect1 <- BNTI_results.melt.transect %>% dplyr::filter(!is.na(value)) %>% unite("Comparison", Var1:Var2, remove = TRUE)
colnames(BNTI_results.melt.transect1)[2] <- "BNTI"

Combined.results.transect <- left_join(BNTI_results.melt.transect1, RaupCrick_results.melt.transect1, by = "Comparison")
Combined.results.transect <- Combined.results.transect %>% filter(BNTI != 0)

Combined.results.transect$Type <- "NDP"

Combined.results.transect$Type[Combined.results.transect$RC > 0.95] <- "Dispersal Limitation"
Combined.results.transect$Type[Combined.results.transect$RC < -0.95] <- "Homogenizing Dispersal"

Combined.results.transect$Type[Combined.results.transect$BNTI > 2] <- "Variable Selection"
Combined.results.transect$Type[Combined.results.transect$BNTI < -2] <- "Homogeneous Selection"
#write.csv(Combined.results.transect, "CHP3Combined.results.transect.csv")

siteData <- data.frame(sample_data(alpha.waterSamp_pru.rare))
coordData <- read.csv("CHP3sites_coordinates.csv")
Combined.results.transect <- Combined.results.transect %>% separate(Comparison, c("site.x", "site.y"), remove = FALSE)
combined.site.data <- left_join(Combined.results.transect, coordData, by = c("site.x" = "SiteID"))
combined.site.data <- left_join(combined.site.data, coordData, by = c("site.y" = "SiteID"))
withinregion.summary <- combined.site.data %>% filter(ecoregion.x == ecoregion.y) %>% group_by(ecoregion.x, Type) %>% tally()
betweenregion.summary <- combined.site.data %>% filter(ecoregion.x != ecoregion.y) %>% group_by(ecoregion.x, ecoregion.y, Type) %>% tally()

write.csv(combined.site.data, "CHP3_CombinedResults-dino-97.csv")


#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#

#### Assembly Processes Map

getNOAA.bathy(lon1 = 145,lon2 = -145, lat1 = -10,lat2 = -90, resolution = 1, antimeridian = T) -> Pacific
coordData$long1 <- ifelse(as.numeric(coordData$long) < 0, (180 + (180 + as.numeric(coordData$long))), as.numeric(coordData$long))

# CHANGE ME: this is where you define the extent of your map, 
# make sure that all sites you are plotting sit within the limits of the image created

results <- combined.site.data
results$long.x1 <- ifelse(as.numeric(results$long.x) < 0, (180 + (180 + as.numeric(results$long.x))), as.numeric(results$long.x))
results$long.y1 <- ifelse(as.numeric(results$long.y) < 0, (180 + (180 + as.numeric(results$long.y))), as.numeric(results$long.y))

p <- autoplot.bathy(Pacific, geom=c("contour", "raster"), colour="white", size=0.1, show.legend = FALSE) +
  scale_fill_gradient2(low="#C2E7FF", mid="#EFEEEC", high="#F0FFFF")


palette <- c("Variable Selection" = "#7652A3", "Homogeneous Selection" = "#FF6F5C", 
             "Homogenizing Dispersal" = "#FFB6AD", "Dispersal Limitation" = "#C8B8DB", 
             "NDP" = "#0D2F6D"
             )

results$Type <- factor(results$Type, levels = c("Variable Selection", "Dispersal Limitation", "Homogeneous Selection", "Homogenizing Dispersal", "NDP")) 

## Plotting assembly processes onto map
map1 <- p + geom_point(data=coordData, aes(y = lat, x = long1), colour="#000000", size = 1) +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_curve(aes(x = long.x1, y = lat.x, xend = long.y1, yend = lat.y, colour = Type), linewidth=1.3, data = results) +
  
  scale_color_manual(values = palette, name="Process") +
  theme(legend.position = c(0.8, 0.2),
        legend.background = element_rect(fill="transparent", colour=NA),
        legend.key = element_rect(colour = NA, fill=NA),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))

map1

# CONGRATS All done! You should now have a map with assembly processes plotted onto it, what a treat :D

#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#
