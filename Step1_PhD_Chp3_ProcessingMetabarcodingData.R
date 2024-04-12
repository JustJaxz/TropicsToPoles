#--------------------------------------------------------#
#--- Cutadapt & Dada2 Pipeline for metabarcoding data ---#
#--------------------------------------------------------#

# Code from Pearman J (john.pearman@cawthron.org.nz)
# Altered by Stuart J for Small sub-unit (SSU) v9 data-set (jacqui.stuart@cawthron.org.nz)

# Date: 2023.05.29
# Primer pair: ill-1380F and ill-1510R
# Sequencing completed by Sequench (https://www.sequench.co.nz/)
# Run numbers CAW-22-11 & CAW-23-21

#------------------------------------------------------------------------------------------------------------------------------#

# Paper title: Tropics to the poles: Diversity and composition of coastal eukaryotic marine microalgae communities across five ecoregions.
# Paper authors: Jacqui Stuart (1,2), Ken Ryan (1), Natalie Robinson (3), Svenja Halfter (3), John K. Pearman (1), 
# Jacob Thomson-Laing (2), Kirsty F. Smith (2)

# Affiliations:
# 1.	School of Biological Sciences, Victoria University of Wellington, PO Box 600, Wellington 6140, New Zealand. 
# 2.	Cawthron Institute, Private Bag 2, Nelson 7042, New Zealand. 
# 3.	National Institute of Water and Atmospheric Research (NIWA), Private Bag 14901, Kilbirnie, Wellington 6241

## Sequence Pipeline

# This R script delineates a comprehensive metabarcoding data processing workflow, encompassing library installation, 
# data organization, primer analysis, Cutadapt-based trimming, Dada2 quality checks, filtering, and chimera removal.
# Utilizing checkpoints ensures adaptability in handling interruptions or errors.

# Subsequently, the script navigates through pivotal steps: "Merging seqtab files," which consolidates tables from
# various sequencing runs based on gene region (v9 or v4), and saves them as .rds files. The "Taxonomy" section processes merged
# tables, metadata, and taxonomy files, crafting phyloseq objects.

# Addressing contamination, the "Subtracting negative controls" segment deducts sequences from negative controls, adjusting
# the main data accordingly. The script then eliminates zero samples, creating a new phyloseq object without negative controls.

# Concluding with a comparison of count sums, the script saves the cleaned dataset as .rds and tracks removed sequences
# in a CSV. It culminates by emphasizing the preparedness of the cleaned data for community and diversity analysis.

#------------------------------------------------------------------------------------------------------------------------------#

#--- Library installation and setup ---#
# Installs BiocManager if don't already have it
# DO NOT INSTALL AGAIN IF YOU HAVE IT...STUFF BREAKS

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

# Install phyloseq & dada2
# DO NOT INSTALL AGAIN IF YOU HAVE IT...STUFF BREAKS
BiocManager::install("phyloseq")
BiocManager::install("dada2", force = TRUE)
BiocManager::install("ShortRead")

# Loading libraries
library("phyloseq")
library("ggplot2")
library("ShortRead")
library("Biostrings")
library("dada2")
library("dplyr")

# Set working directory
### CHANGE ME ### This will be unique to your project, make sure it is correct
setwd("C:/Users/user/Documents/CAW-22-11") # Sequence files 1
#setwd("C:/Users/user/Documents/CAW-23-21") # Sequence files 2


#--------------------#
##### DATA SETUP #####
#--------------------#

# If you have received your sequence data nested in multiple levels of folders
# the next step saves you a lot of copy and pasting and puts all the '.fastq' files
# into a single folder. 


#--- Get files out of sub folders ---#

###CHANGE ME###
# STEP 1: Set source destination - this needs to be the folder all the files and folders are nested in
source_folder <- "C:/Users/user/Documents/CAW-22-11/v9" # Sequence files 1
#source_folder <- "C:/Users/user/Documents/CAW-23-21/v9" # Sequence files 2

# STEP 2: Make and Set destination folder - This is the folder you want all the '.fastq' files to end up in
seq_folder_name <- "sequence-files"
dir.create(seq_folder_name)
destination_folder <- "C:/Users/user/Documents/CAW-22-11/sequence-files"
#destination_folder <- "C:/Users/user/Documents/CAW-23-21/sequence-files"

# STEP 3: Get list of files in source folder and subfolders:
file_list <- list.files(source_folder, recursive = TRUE, full.names = TRUE)

# STEP 4: Loop through each file - This finds and moves all the nested files to their new home
for (file in file_list) {
  # Construct the destination path by removing the source folder path
  destination_path <- file.path(destination_folder, basename(file))
  # Move the file to the destination folder
  file.rename(file, destination_path)
}
# File sorting COMPLETE


#--- Delete now empty sub folders ---#
# You don't have to do this step, but as a result of the last one you will now have a whole series of empty
# folders and sub-folder. For the humans out there that like a wee tidy, this goes out to you.

# STEP 1: Get the list of sub-directories in the source folder
# After the file moving loop, the code retrieves a list of sub-directories in the source folder
# Using list.dirs() with recursive = FALSE.
subdirs <- list.dirs(source_folder, recursive = TRUE, full.names = TRUE)

# STEP 2: Sort the sub-directories 
# This sorts remaining sub-directories in reverse order using 'sort()' 
# This is to ensure that the deepest ones are deleted first
subdirs <- sort(subdirs, decreasing = TRUE)

# STEP 3: This is a loop through each subdirectory
for (subdir in subdirs) {
  # Check if the subdirectory is empty
  if (length(list.files(subdir, recursive = TRUE)) == 0) {
    # Delete the empty subdirectory
    unlink(subdir, recursive = TRUE, force = TRUE)
  }
}

# Folder deleting COMPLETE, all neat and tidy



#--- Check files ---#
### CHANGE ME ### This will direct R to your sequence data in your new 'destination_folder', make sure it is correct
path <- "C:/Users/user/Documents/CAW-22-11/sequence-files"
list.files(path)

# List names of forward and reverse reads
fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))

# View the list of names
fnFs
fnRs


#--- Primer pair setup and checking ---#

# This code specifies the primer pair used for your sequencing run. 
### CHANGE ME ### Make sure these are the primers you used.

FWD <- "CCCTGCCHTTTGTACACAC" #FORWARD: 1380F v9 
REV <- "CCTTCYGCAGGTTCACCTAC" #REVERSE: 1510R v9


#Here we find the forward and reverse compliments of primers in all sequences
allOrients <- function(primer) {
  # Then create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients


# This says how many times the primers are found, and in what orientation
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))

# At this point you are hoping for all reads, or the vast majority, to be in the forward direction.
# If this is not the case, find Johnny P and grill him about why...it's most likely a you problem
# But he will know where to direct you to solve it.


#--------------------------------------#
##### PROCESSING STEP 1 - Cutadapt #####
#--------------------------------------#


#--- finds and sets up cutadapt --- #

### CHANGE ME ### 
# Please change the following file path to the cutadapt path on your machine
cutadapt <- "C:/Users/user/Miniconda3/Scripts/cutadapt" 

# Run shell commands from R
system2(cutadapt, args = "--version") 


# This code determines the new file path for the cutadapt folder/files you are creating
path.cut <- file.path(path, "cutadapt18S")
if(!dir.exists(path.cut)) dir.create(path.cut) #does it exist, if it doesn't then make it
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))



#--- Trimming primer sequences off all reads ---#

# Trim FWD off of R1 (forward reads) - 
R1.flags <- paste0("-g", " ^", FWD) 

# Trim REV off of R2 (reverse reads)
R2.flags <- paste0("-G", " ^", REV) 


#--- Run Cutadapt ---#
# This step can take some time

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-e 0.08 --discard-untrimmed", R1.flags, R2.flags,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}

# You will now have a new 'cutadapt18s' folder within the folder your sequences are in.
# These have had the primers trimmed off each sequence.
# Next you are making sure the primer sequences have been removed and that you still have the same amount of files.

# Re-Check how many files still contain the primers after running cutadapt: 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


# Identify the forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

# Check that the number of forward and reverse files are still the same:
if(length(cutRs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if (length(cutRs) != length(cutRs)) stop("Forward and reverse files do not match. Better go back and have a check")

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#Check sample names
sample.names



#----------------------#
##### STOP POINT 1 #####
#----------------------#
# This stop point is designed to save the work-space environment so you cab pick up from here later if needed.
# It means you do not have to run the above code again and also just worth saving sometimes in case shit goes awry


# generate folder in working directory for saved environments
folder_name <- "workspace_environments"
dir.create(folder_name)

# generate dynamic file name to track iterations
current_datetime <- format(Sys.time(), "%Y-%m-%d_%H-%M")
file_name_stop1 <- paste0("C:/Users/user/Documents/CAW-22-11.1/workspace_environments/",current_datetime,"_CAW-22-11.1_stop1", ".RData")

# save workspace environment
save.image(file = file_name_stop1)



#-------------------------------------------------#
##### PROCESSING STEP 2 - Dada2 quality check #####
#-------------------------------------------------#


#--- LOAD WORKSPACE ENVIRONMENT ---#
### CHANGE ME ###
# Make sure you update this file link as you save out your stop point files.
#load("C:/Users/user/Documents/CAW-22-11.1/workspace_environments/2023-05-29_14-38_CAW-22-11.1_stop1.RData")



#--- Dada2 - Quality Checking ---#

require(dada2)
require(ggplot2)

if(length(cutFs) <= 20) {
  fwd_qual_plots <- plotQualityProfile(cutFs) +
    scale_x_continuous(breaks=seq(0,250,10)) +
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  rev_qual_plots <- plotQualityProfile(cutRs) +
    scale_x_continuous(breaks=seq(0,250,10)) +
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
} else {
  rand_samples <- sample(size = 20, 1:length(cutFs)) # grab 20 random
  fwd_qual_plots <- plotQualityProfile(cutFs[rand_samples]) +
    scale_x_continuous(breaks=seq(0,250,10)) +
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  rev_qual_plots <- plotQualityProfile(cutRs[rand_samples]) +
    scale_x_continuous(breaks=seq(0,250,10)) +
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


# Shows the plotted quality profiles
fwd_qual_plots
rev_qual_plots  

# Good idea to export plots as images for later, you never know when you might need to whip them out.


#----------------------#
##### STOP POINT 2 #####
#----------------------#

# This stop point is designed to save the work-space environment so you cab pick up from here later if needed.
# It means you do not have to run the above code again and also just worth saving sometimes in case shit goes awry

#generate dynamic file name to track iterations
current_datetime <- format(Sys.time(), "%Y-%m-%d_%H-%M") # This tags the date and time into the generated file
#This creates a file name with the date and time tag for this stop point
file_name_stop2 <- paste0("C:/Users/user/Documents/CAW-22-11.1/workspace_environments/",current_datetime,"_CAW-22-11.1_stop2", ".RData")
#saves work space environment
save.image(file = file_name_stop2)


#--------------------------------------------------#
##### PROCESSING STEP 3 - Filtering & trimming #####
#--------------------------------------------------#

#--- LOAD WORKSPACE ENVIRONMENT ---#
### CHANGE ME ###
# Make sure you update this file link as you save out your stop point files.

#load("C:/Users/user/Documents/CAW-22-11.1/workspace_environments/2023-05-29_15-28_CAW-22-11.1_stop2")
#--- Filtering and trimming ---#

filtpathF <- file.path(path.cut, "filtered", basename(cutFs))
filtpathR <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtpathF, cutRs, filtpathR,
                     ### CHANGE ME ### 'truncLen = ' This is the expected length of your sequences/what you are trimming them to.
                     # v9 = 150 max length
                     truncLen=c(130,130), maxEE=c(2,4), truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

out

sample.names <- sapply(strsplit(basename(filtpathF), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtpathR), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(identical(sample.names, sample.namesR)) {print("Files are still matching.....congratulations")
} else {stop("Forward and reverse files do not match.")}
names(filtpathF) <- sample.names
names(filtpathR) <- sample.namesR


set.seed(100) # set seed to ensure that randomized steps are replicable


#--- Learn forward error rates ---#
# This can take time, go get a coffee, get out for a run or something
errF <- learnErrors(filtpathF, nbases=1e8, multithread=TRUE, verbose = TRUE)
## 100285500 total bases in 668570 reads from 24 samples will be used for learning the error rates.


#--- Learn reverse error rates ---#
# This will also take time, so maybe crack into some writing since you're all caffeinated or exercised.
errR <- learnErrors(filtpathR, nbases=1e8, multithread=TRUE, verbose = TRUE)

#--- Plotting forward error rates ---#
#Forward
errF_plot <- plotErrors(errF, nominalQ=TRUE)
errF_plot

#--- Plotting reverse error rates ---#
#Reverse
errR_plot <- plotErrors(errR, nominalQ=TRUE)
errR_plot


##### STOP POINT 3 #####
#--- Here you can save the work-space environment and pick up from here later if needed ---#
#--- This means you do not have to run the above code again
#--- Also just worth saving sometimes incase computer crashes


#generate dynamic file name to track iterations
current_datetime <- format(Sys.time(), "%Y-%m-%d_%H-%M")
file_name_stop3 <- paste0("C:/Users/user/Documents/CAW-22-11.1/workspace_environments/",current_datetime,"_CAW-22-11.1_stop3", ".RData")

#save workspace environment
save.image(file = file_name_stop3)



#----------------------------------------------#
##### PROCESSING STEP 4 - Chimera removal ######
#----------------------------------------------#

#--- LOAD WORKSPACE ENVIRONMENT ---#
### CHANGE ME ###
# Make sure you update this file link as you save out your stop point files.

#load("C:/Users/user/Documents/CAW-22-11.1/workspace_environments/2023-05-29_17-51_CAW-22-11.1_stop3.RData")


#----- Chimera removal -----#
# This code performs the de-replication of sequences in a FASTQ file specified by 'filtpathF' using the 'derepFastq()' function. 
# The resulting de-replicated data is stored in the variable 'derepF' or 'derepR'. 
# The 'verbose=TRUE' argument is used to display detailed progress information during the de-replication process.
derepF <- derepFastq(filtpathF, verbose=TRUE)
derepR <- derepFastq(filtpathR, verbose=TRUE)


#--- Checkpoint 1 ---#
#--- This is just to ensure saving of those longer processes, as sometimes the computer crashes after this point.
#generate dynamic file name to track iterations
current_datetime <- format(Sys.time(), "%Y-%m-%d_%H-%M")
file_name_stop3_checkpoint <- paste0("C:/Users/user/Documents/CAW-22-11.1/workspace_environments/",current_datetime,"_CAW-22-11.1_stop3-checkpoint", ".RData")

#save workspace environment
save.image(file = file_name_stop3_checkpoint)

#--------------------#


# This code does the de-noising and chimera removal of amplicon sequence data using the 'dada()' function from the dada2 package. 
# The de-replicated data ('derepF' or 'derepR'), an error model ('errF' or 'errR'), and other parameters are provided as input to the function. 
# The de-noised sequences and related information are stored in the 'dadaF.pseudo' or 'dadaR.pseudo' object.
dadaF.pseudo <- dada(derepF, err=errF, multithread=TRUE, pool="pseudo")
# THIS BIT TAKES A LOOOOOONG TIME ^


dadaR.pseudo <- dada(derepR, err=errR, multithread=TRUE, pool="pseudo")
# THIS BIT TAKES A LOOOOOONG TIME ^


# Merging of paired-end amplicon sequence reads using the 'mergePairs()' function from the dada2 package. 
# The de-noised and de-replicated forward and reverse reads, along with specified parameters, are provided as input. 
# The merged sequences and related information are stored in the mergers object.
mergers <- mergePairs(dadaF.pseudo, derepF, dadaR.pseudo, derepR, maxMismatch = 1, minOverlap = 10, verbose=TRUE)
# THIS TAKES A WHILE TOO  ^


# line of code creates a sequence table (seqtab) from the merged amplicon sequences stored in the mergers object. 
# The sequence table provides information about the abundance of each unique sequence in the merged data, 
# allowing further analysis and downstream processing.
seqtab <- makeSequenceTable(mergers)

# Splits the name of a file or directory specified by path using the "-" character as the delimiter. 
# It then extracts the first two elements from each split part and stores them in the split.dir variable. 
split.dir <- sapply(strsplit(basename(path), "-"), `[`, 1:2)

# Combines the first and second elements from the `split.dir` variable into a single string, 
# separated by a period, and stores it in the `split.dir.name` variable.
split.dir.name  <- paste(split.dir[1], split.dir[2], sep=".")

# saves the R object `seqtab` as an RDS file, using the constructed file path and name that
# includes the directory path, base name, and a combination of split parts from `split.dir.name`.
saveRDS(seqtab, paste0(path, "seqtab", split.dir.name,".rds"))

# generates a frequency table that shows the distribution of sequence lengths in the `seqtab` object. 
# It counts the number of sequences with each length and provides a summary of their occurrence
table(nchar(getSequences(seqtab)))

# generates a histogram that visualizes the distribution of sequence lengths in the `seqtab` object. 
# It plots the frequency of different sequence lengths on the x-axis and the count or density on the y-axis, 
# providing an overview of the distribution pattern
hist(nchar(getSequences(seqtab)), main="Distribution of sequence lengths")

# filters the columns of the `seqtab` object based on the character length of their column names.
# It selects and retains only the columns with column names whose length is within your designated range.
### CHANGE ME ### the numbers in `seq(000,000)` to represent the expected character range for your sequences. 
# The resulting subset is assigned back to the `seqtab` object, effectively updating it with the filtered columns.
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% seq(120,142)]

#  applies the `removeBimeraDenovo()` function to the `seqtab` object to remove chimeric sequences. 
# It utilizes multiple threads for faster computation and provides detailed output during the process. 
# The resulting sequence table without chimeras is assigned to the `seqtab.nochim` variable.
seqtab.nochim <- removeBimeraDenovo(seqtab, multithread=TRUE, verbose=TRUE)

### CHANGE ME ###
# make sure I match your file path to your working folder for this project
# You are saving a copy of the chimera free sequence table to your working directory
saveRDS(seqtab.nochim, "C:/Users/user/Documents/CAW-22-11.1/CAW-22-11.1_jacqui_v9_seqtab.nochim.rds")


### At this point, if you are not doing your own taxonomy files, send the .rds file to John.

##### STOP POINT 4 #####
#--- Here you can save the work-space environment and pick up from here later if needed ---#
#--- This means you do not have to run the above code again
#--- Also just worth saving sometimes incase computer crashes


#generate dynamic file name to track iterations
current_datetime <- format(Sys.time(), "%Y-%m-%d_%H-%M")
file_name_stop4 <- paste0("C:/Users/user/Documents/CAW-22-11.1/workspace_environments/",current_datetime,"_CAW-22-11.1_stop4", ".RData")

#save workspace environment
save.image(file = file_name_stop4)



#------------------------------------------------------------------------------------------------------------------------------#
### Mergin seqtab files ###

# This is to be run when you have multiple sequencing runs that were split to run through the your pipeline
# Only merge seqtab files of the same gene region
# Next step after this is taxonomy

#--- Set working directory ---#
setwd("C:/Users/user/Documents/R")

#--- Read in seqtab files ---#
# Here I am defining the name of each data table `seqtab.nochim.file.name` 
# The reading in the existing sequence table for each sequence run using `readRDS()`
# make sure to change the file path to match where your files are!

# v9 region sequence tables
seqtab.nochim.CAW.22.11 <- readRDS("C:/Users/user/Documents/R/CAW-22-11_jacqui_v9_seqtab.nochim.rds")
seqtab.nochim.CAW.23.21 <- readRDS("C:/Users/user/Documents/R/CAW-23-21_jacqui_v9_seqtab.nochim.rds")

# v4 region sequence tables
seqtab.nochim.CAW.22.12 <- readRDS("C:/Users/user/Documents/R/CAW-22-12_jacqui_v4_seqtab.nochim.rds")
seqtab.nochim.CAW.23.22 <- readRDS("C:/Users/user/Documents/R/CAW-23-22_jacqui_v4_seqtab.nochim.rds")

#----------------------------#
##### Merge seqtab files #####
#----------------------------#

# define input table groups
v9_tables <- list(seqtab.nochim.CAW.22.11, seqtab.nochim.CAW.23.21)
v4_tables <- list(seqtab.nochim.CAW.22.12, seqtab.nochim.CAW.23.22)

# Merge all v9 region sequence tables
v9.seqtab.merge <- mergeSequenceTables(tables = v9_tables)
#check new file rownames order
rownames(v9.seqtab.merge)

# Merge all v4 region sequence tables
v4.seqtab.merge <- mergeSequenceTables(tables = v4_tables)
#check new file rownames order
rownames(v4.seqtab.merge)

# Save `.rds` file of each merged table to working directory
saveRDS(v9.seqtab.merge, "C:/Users/user/Documents/R/2023-06-01_jacqui_v9_merged-seqtab.nochim.rds")
saveRDS(v4.seqtab.merge, "C:/Users/user/Documents/R/2023-06-01_jacqui_v4_merged-seqtab.nochim.rds")

# Now you should have a new merged sequence table in your working directory, or the filepath you defined above
# Woohoo, all done.


#------------------------------------------------------------------------------------------------------------------------------#
### Taxonomy ###
# Here we combine the taxonomy and sample metadata files and
# remove negative control sample data from the data set. Then you will have a delightful clean data-set
# You need to have a metadata file at this point that will have all the base information associated with your sequenced samples. 


#--- Read in data ---#
### CHANGE ME: Make sure file path and name is correct!!!

#--- v9 ---#
# read in `seqtab.nochim` file from your working directory
v9.seqtab.merge <- readRDS("C:/Users/user/Documents/R/2023-06-01_jacqui_v9_merged-seqtab.nochim.rds")
#read in metadata .csv file
v9.metadata <- read.csv("C:/Users/user/Documents/R/JSPhd_sample-metadata.csv", header = TRUE, row.names = 1)
# read in taxonomy .rds file
v9.taxonomy <- readRDS("C:/Users/user/Documents/R/v9.taxonomy.v5.rds")

#--- v4 ---#
# read in `seqtab.nochim` file from your working directory
v4.seqtab.merge <- readRDS("C:/Users/user/Documents/R/2023-06-01_jacqui_v4_merged-seqtab.nochim.rds")
#read in metadata .csv file
v4.metadata <- read.csv("C:/Users/user/Documents/R/JSPhd_sample-metadata.csv", header = TRUE, row.names = 1)
# read in taxonomy .rds file
v4.taxonomy <- readRDS("C:/Users/user/Documents/R/v4.taxonomy.v5.rds")

#--- Rename row names ---#
# rename `seqtab.nochim` rows - NOTE there is a way to just get row-names from the metadata file, commented out below.
# make sure these match the sample order in your sequence table file.

#--- v9 ---#
#re-check sample rowname order
rownames(v9.seqtab.merge)

#define new row names here
v9.new.rownames <- c(
  #CAW.22.11
  "N1a1", "N1a2", "N1a3", "N1b1", "N1b2", "N1b3", "N1c1", "N1c2", "N1c3", "N1d1", "N1d2", "N1d3", "N1e1", "N1e2", "N1e3", "N1f1", "N1f2", "N1f3", "N1g1", "N1g2", "N1g3", "N1h1", "N1h2", "N1h3", "N1i1", "N1i2", "N1i3",
  "N2a1", "N2a2", "N2a3", "N2b1", "N2b2", "N2b3", "N2c1", "N2c2", "N2c3", "N2d1", "N2d2", "N2d3", "N2e1", "N2e2", "N2e3", "N2f1", "N2f2", "N2f3", "N2g1", "N2g2", "N2g3", "N2h1", "N2h2", "N2h3", "N2i1", "N2i2", "N2i3",
  "N3a1", "N3a2", "N3a3", "N3b1", "N3b2", "N3b3", "N3c1", "N3c2", "N3c3", "N3d1", "N3d2", "N3d3", "N3e1", "N3e2", "N3e3", "N3f1", "N3f2", "N3f3", "N3g1", "N3g2", "N3g3", "N3h1", "N3h2", "N3h3", "N3i1", "N3i2", "N3i3",
  "N4a1", "N4a2", "N4a3", "N4b1", "N4b2", "N4b3", "N4c1", "N4c2", "N4c3",
  "N4d1", "N4d2", "N4d3", "ExtBL1", "PCRBL1", "H20BL1" , "N4e1", "N4e2", "N4e3", "N4f1", "N4f2", "N4f3", "N4g1", "N4g2", "N4g3", "N4h1", "N4h2", "N4h3", "N4i1", "N4i2", "N4i3", 
  "N5a1", "N5a2", "N5a3", "N5b1", "N5b2", "N5b3", "N5c1", "N5c2", "N5c3", "N5d1", "N5d2", "N5d3", "N5e1", "N5e2", "N5e3", "N5f1", "N5f2", "N5f3", "N5g1", "N5g2", "N5g3", "N5h1", "N5h2", "N5h3", "N5i1", "N5i2", "N5i3", 
  "N6a1", "N6a2", "N6a3", "N6b1", "N6b2", "N6b3", "N6c1", "N6c2", "N6c3", "N6d1", "N6d2", "N6d3", "N6e1", "N6e2", "N6e3", "N6f1", "N6f2", "N6f3", "N6g1", "N6g2", "N6g3", "N6h1", "N6h2", "N6h3", "N6i1", "N6i2", "N6i3", 
  "H01", "H03", "H25", "H26", "H32", "H37", "H50", "H53", "H84", "H95", "H96", "H101", "H102", "H105", "H108", "H115", "H116", "H119", "H120", "H121", "H124", "H159", "H165", "H171", "ExtBL2", "PCRBL2", "H20BL2",
  "H172", "H178", "H216", "H252", "21_AP1", "21_AP2", "21_AP3", "21_AP8", "21_AP9", "21_AP10", "21_AP11", "H167", "ExtBL3", "PCRBL3", "H20BL3",
  #CAW.23.21
  "P2a-00", "P2a-25", "P2a-50", "P2a-75", "P2a-100", "P2a-125", "P2a-150", "P2a-175", "P2a-200", "P2a-225", "P2b-00", "P2b-25", "P2b-50", "P2b-75", "P2b-100", "P2b-125", "P2b-150", "P2b-175",
  "P2c-00", "P2a-CT", "P2a-CM", "P2a-CB", "P2b-CT", "P2b-CM", "P2b-CB", "P2c-CT", "P2c-CM", "P2c-CB", "P2a1", "P2a2", "P2a3", "P2b1", "P2b2", "P2b3", "P1a-00", "P1a-25",
  "P1a-50", "P1a-75", "P1a-100", "P1b-00", "P1b-25", "P1b-50", "P1b-75", "P1b-100", "P1a-CT", "P1a-CM", "P1a-CB", "P1b-CT", "P1b-CM", "P1b-CB", "P1c-CT", "P1c-CM", "P1c-CB", "P4a1",
  "P4a2", "P4a3", "P4a4", "P3a-CT", "P3a-CM", "P3a-CB", "P3b-CT", "P3b-CM", "P3b-CB", "P3c-CT", "P3c-CM", "P3c-CB", "P1a1", "P1a2", "P1a3", "P1b1", "P1b2", "P1b3",
  "P1c1", "P1c2", "P1c3", "P1d1", "P1d2", "P1d3", "P1e1", "P1e2", "P1e3", "C1a1", "C1a2", "C1a3", "C1b1", "C1b2", "C1b3", "C1c1", "C1c2", "C1c3", "EXT_BL","PCR_BL1","H20_BL1",
  "T1a1", "T1a2", "T1a3", "T1b1", "T1b2", "T1b3", "T1c1", "T1c2", "T1c3", "T1d1", "T1d2", "T1d3", "T1e1", "T1e2", "T1e3", "T1f1", "T1f2", "T1f3", "T1g1", "T1g2", "T1g3", "T1h1", "T1h2", "T1h3", "T1i1", "T1i2", "T1i3",
  "T2a1", "T2a2", "T2a3", "T2b1", "T2b2", "T2b3", "T2c1", "T2c2", "T2c3", "T2d1", "T2d2", "T2d3", "T2e1", "T2e2", "T2e3", "T2f1", "T2f2", "T2f3", "T2g1", "T2g2", "T2g3", "T2h1", "T2h2", "T2h3", "T2i1", "T2i2", "T2i3",
  "T3a1", "T3a2", "T3a3", "T3b1", "T3b2", "T3b3", "T3c1", "T3c2", "T3c3", "T3d1", "T3d2", "T3d3", "T3e1", "T3e2", "T3e3", "T3f1", "T3f2", "T3f3", "T3g1", "T3g2", "T3g3", "T3h1", "T3h2", "T3h3", "T3i1", "T3i2", "T3i3",
  "EXT_BL6","PCR_BL3", "H20_BL3"
)

#change rownames here
rownames(v9.seqtab.merge) <- v9.new.rownames


#--- v4 ---#
#re-check sample rowname order
rownames(v4.seqtab.merge)
#define new row names here
v4.new.rownames <- c(
  #CAW.22.12
  "N1a1", "N1a2", "N1a3", "N1b1", "N1b2", "N1b3", "N1c1", "N1c2", "N1c3", "N1d1", "N1d2", "N1d3", "N1e1", "N1e2", "N1e3", "N1f1", "N1f2", "N1f3", "N1g1", "N1g2", "N1g3", "N1h1", "N1h2", "N1h3", "N1i1", "N1i2", "N1i3",
  "N2a1", "N2a2", "N2a3", "N2b1", "N2b2", "N2b3", "N2c1", "N2c2", "N2c3", "N2d1", "N2d2", "N2d3", "N2e1", "N2e2", "N2e3", "N2f1", "N2f2", "N2f3", "N2g1", "N2g2", "N2g3", "N2h1", "N2h2", "N2h3", "N2i1", "N2i2", "N2i3",
  "N3a1", "N3a2", "N3a3", "N3b1", "N3b2", "N3b3", "N3c1", "N3c2", "N3c3", "N3d1", "N3d2", "N3d3", "N3e1", "N3e2", "N3e3", "N3f1", "N3f2", "N3f3", "N3g1", "N3g2", "N3g3", "N3h1", "N3h2", "N3h3", "N3i1", "N3i2", "N3i3",
  "N4a1", "N4a2", "N4a3", "N4b1", "N4b2", "N4b3", "N4c1", "N4c2", "N4c3", "N4d1", "N4d2", "N4d3", "v4_ExtBL1", "v4_PCRBL1", "v4_H20BL1" , "N4e1", "N4e2", "N4e3", "N4f1", "N4f2", "N4f3", "N4g1", "N4g2", "N4g3", "N4h1", "N4h2", "N4h3", "N4i1", "N4i2", "N4i3", 
  "N5a1", "N5a2", "N5a3", "N5b1", "N5b2", "N5b3", "N5c1", "N5c2", "N5c3", "N5d1", "N5d2", "N5d3", "N5e1", "N5e2", "N5e3", "N5f1", "N5f2", "N5f3", "N5g1", "N5g2", "N5g3", 
  "N5h1", "N5h2", "N5h3", "N5i1", "N5i2", "N5i3", "N6a1", "N6a2", "N6a3", "N6b1", "N6b2", "N6b3", "N6c1", "N6c2", "N6c3", "N6d1", "N6d2", "N6d3", "N6e1", "N6e2", "N6e3", "N6f1", "N6f2", "N6f3", "N6g1", "N6g2", "N6g3", "N6h1", "N6h2", "N6h3", "N6i1", "N6i2", "N6i3", 
  "H01", "H03", "H25", "H26", "H32", "H37", "H50", "H53", "H84", "H95", "H96", "H101", "H102", "H105", "H108", "H115", "H116", "H119", "H120", "H121", "H124", "H159", "H165", "H171", "ExtBL2", "PCRBL2", "H20BL2",
  "H172", "H178", "H216", "H252", "21_AP1", "21_AP2", "21_AP3", "21_AP8", "21_AP9", "21_AP10", "21_AP11", "H167", "ExtBL3", "PCRBL3", "H20BL3",
  #CAW.23.22
  "P2a-00", "P2a-25", "P2a-50", "P2a-75", "P2a-100", "P2a-125", "P2a-150", "P2a-175", "P2a-200", "P2a-225", "P2b-00", "P2b-25", "P2b-50", "P2b-75", "P2b-100", "P2b-125", "P2b-150", "P2b-175",
  "P2c-00", "P2a-CT", "P2a-CM", "P2a-CB", "P2b-CT", "P2b-CM", "P2b-CB", "P2c-CT", "P2c-CM", "P2c-CB", "P2a1", "P2a2", "P2a3", "P2b1", "P2b2", "P2b3", "P1a-00", "P1a-25",
  "P1a-50", "P1a-75", "P1a-100", "P1b-00", "P1b-25", "P1b-50", "P1b-75", "P1b-100", "P1a-CT", "P1a-CM", "P1a-CB", "P1b-CT", "P1b-CM", "P1b-CB", "P1c-CT", "P1c-CM", "P1c-CB", "P4a1",
  "P4a2", "P4a3", "P4a4", "P3a-CT", "P3a-CM", "P3a-CB", "P3b-CT", "P3b-CM", "P3b-CB", "P3c-CT", "P3c-CM", "P3c-CB", "P1a1", "P1a2", "P1a3", "P1b1", "P1b2", "P1b3",
  "P1c1", "P1c2", "P1c3", "P1d1", "P1d2", "P1d3", "P1e1", "P1e2", "P1e3", "C1a1", "C1a2", "C1a3", "C1b1", "C1b2", "C1b3", "C1c1", "C1c2", "C1c3", "EXT_BL","PCR_BL1","H20_BL1",
  "T1c1", "T1c2", "T1c3", "T1d1", "T1d2", "T1d3", "T1e1", "T1e2", "T1e3", "T1f1", "T1f2", "T1f3", "T1g1", "T1g2", "T1g3", "T1h1", "T1h2", "T1h3", "T1i1", "T1i2", "T1i3",
  "T2a1", "T2a2", "T2a3", "T2b1", "T2b2", "T2b3", "T2c1", "T2c2", "T2c3", "T2d1", "T2d2", "T2d3", "T2e1", "T2e2", "T2e3", "T2f1", "T2f2", "T2f3", "T2g1", "T2g2", "T2g3", "T2h1", "T2h2", "T2h3", "T2i1", "T2i2", "T2i3",
  "T3a1", "T3a2", "T3a3", "T3b1", "T3b2", "T3b3", "T3c1", "T3c2", "T3c3", "T3d1", "T3d2", "T3d3", "T3e1", "T3e2", "T3e3", "T3f1", "T3f2", "T3f3", "T3g1", "T3g2", "T3g3", "T3h1", "T3h2", "T3h3", "T3i1", "T3i2", "T3i3",
  "EXT_BL6","PCR_BL3", "H20_BL3", "T1a1", "T1a2", "T1a3", "T1b1", "T1b2", "T1b3"
)
#change rownames here
rownames(v4.seqtab.merge) <- v4.new.rownames

#--- Merging data ---#
#combine metadata, sample information and metabarcoding files
v9 <- phyloseq(otu_table(v9.seqtab.merge, taxa_are_rows = FALSE), sample_data(v9.metadata), tax_table(v9.taxonomy))
v4 <- phyloseq(otu_table(v4.seqtab.nochim, taxa_are_rows = FALSE), sample_data(v4.metadata), tax_table(v4.taxonomy))


#---------------------------------------#
##### Subtracting negative controls #####
#---------------------------------------#

# This code segment subtracts sequences found in the negative control samples
# from all samples to account for contamination.


#--- V9 region ---#
# subsets neg control samples in merged data file
v9.negatives = subset_samples(v9, type =="control")
v9.negatives

# calculates the column sums for each variable in the `v9.negatives` OTU table and stores them in the `v9.negSums` vector.
v9.negSums <- colSums(otu_table(v9.negatives))

# Converts `v9.negSums` to a vector
v9.negSums_vec <- as.vector(v9.negSums)

# converts the object `otu_table(v9.negatives)` into a matrix and assigns it to the variable `v9.table`
v9.table = as(otu_table(v9.negatives), "matrix")

# converts the matrix object `v9.table` into a data frame and assigns it to the variable `v9.tabledf`
v9.tabledf = as.data.frame(v9.table)
v9.tabledf[,1:length(v9.tabledf)] <- sweep(v9.tabledf[,1:length(v9.tabledf)], 2,v9.negSums_vec)

# changes neg controls to zero
v9.tabledf <- replace(v9.tabledf, v9.tabledf < 0,0)

# removes zero samples (or negative controls)
v9.tabledf_noneg <- v9.tabledf[rowSums(v9.tabledf)!=0,]

# Creating a new phyloseq object with negative controls removed from all samples
# Then combining it with the appropriate metadata and taxonomy information
v9.cleaned <- phyloseq(otu_table(v9.tabledf_noneg, taxa_are_rows = FALSE),
                       sample_data(v9.metadata),
                       tax_table(v9.taxonomy))

v9.cleaned


#--------------------------------------------------#
##### Comparing original and clean sample sums #####
#--------------------------------------------------#

# calculates the total count sums for each sample in two different phyloseq objects (v9 and v9.cleaned) and 
# stores them in separate data frames. It then performs a left join operation to merge the sample sums from 
# both objects into a single data frame (control.rem.sums) for comparison.


org.ss <- as.data.frame(sample_sums(v9)) # calculate total count sums for each sample
org.ss$Names <- rownames(org.ss) # adds new column to data frame called names containing the sample names

new.ss <- as.data.frame(sample_sums(v9.cleaned)) # calculate total count sums for each sample
new.ss$Names <- rownames(new.ss) # adds new column to data frame called names containing the sample names

control.rem.sums <- dplyr::left_join(org.ss, new.ss, by="Names") # merging data frames based on matching `Name` values

colnames(control.rem.sums) <- c("Original", "Names", "New") # defining column names

control.rem.sums <- control.rem.sums %>%
  mutate(perc = New/Original*100) # finding the percent change between original and cleaned

control.rem.sums



#--- Save files ---#
#Saving clean data set to working directory
saveRDS(v9.cleaned,"C:/Users/user/Documents/R/v9_cleaned.rds")


#------------------------------------------------------------------------------------------------------------------------------------------------#
# track sequences that have been removed at each step

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaF.pseudo, getN), sapply(dadaR.pseudo, getN), sapply(mergers, getN))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names
track <- data.frame(track)
track$Perc <- track$merged/track$input*100

write.csv(track, "track.csv")

#------------------------------------------------------------------------------------------------------------------------------------------------#
# THE END - Now you can do your community and divesrity analysis, woo #
