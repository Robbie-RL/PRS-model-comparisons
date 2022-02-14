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


## Create bigSNP file from QC Base plink files
# Only needs to be performed once!
snp_readBed("EUR.QC.final.bed") #create rds



# Relevant data and variables
base.bigSNP <- snp_attach("EUR.QC.stacking.rds")
G <- base.bigSNP$genotypes #[0,1,2] for each SNP, for each person (matrix)
CHR <- base.bigSNP$map$chromosome #Chr of each SNP (vector)
POS <- base.bigSNP$map$physical.pos #Pos of each SNP (vector)
NCORES <- nb_cores() 