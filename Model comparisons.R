###########################
## PRS model comparisons ##
###########################

## 4 comparisons
# 1. Standard PRS (default clumping) optimised for best p-value threshold
# 2. Stacked clumping and thresholding
# 3. Best clumping and thresholding combination
# 4. Lassosum


library(bigsnpr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)
library(lassosum)
library(methods)
library(parallel)






#################################################
## Load in sample data and relevant variables  ##
#################################################


## Data cleaning
## Remove sample_ID with missing height data
## This step only needs to be performed once!
base.id <- read.table("EUR.QC.fam") #Sample ID
base.id <- base.id[, c(1,2)]
colnames(base.id) <- c("FID", "IID")
tgt.data <- read.table("EUR.height", header=T) #Sample phenotype (height)
valid <- merge(tgt.data, base.id, by=c("FID", "IID"))
valid.ind <- which(base.id$FID %in% valid$FID)
## Create bigSNP file rds from bed file
snp_readBed2("EUR.QC.bed", ind.row = valid.ind)









# Relevant data and variables
base.bigSNP <- snp_attach("EUR.QC.rds")
G <- base.bigSNP$genotypes #[0,1,2] for each SNP, for each person (matrix)
CHR <- base.bigSNP$map$chromosome #Chr of each SNP (vector)
POS <- base.bigSNP$map$physical.pos #Pos of each SNP (vector)
NCORES <- nb_cores() 




## Load in sample heights
base.height <- read.table("EUR.height", header=T) 
# Merge with base data and remove samples with missing height
pheno <- merge(base.bigSNP$fam[c(1,2)], base.height, by=c(1,2))
colnames(pheno)[1:2] <- c("FID", "IID")




## Load in summary statistics. Merge with base data
sumstats <- read.table(gzfile("Height.QC.gz"), header=T)
sumstats <- sumstats[, c(1, 3, 2, 4, 5, 9, 8)]
# Convert odds ratio to beta
sumstats$OR <- log(sumstats$OR)
# Update and normalise names
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "p")
map <- base.bigSNP$map[,-(2:3)] #base map data
names(map) <- c("chr", "pos", "a0", "a1") #rename columns to standardise
info_snp <- snp_match(sumstats, map)




## Set up training/testing split datasets
set.seed(42)
# From 559 individuals, randomly choose 330 (70/30 split)
ind.train <- sample(nrow(G), 330)
# Testing data is everyone not in training data
ind.test <- setdiff(rows_along(G), ind.train)



