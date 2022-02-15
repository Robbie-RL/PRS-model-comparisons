###########################
## PRS model comparisons ##
###########################

## 4 comparisons
# 1. Standard PRS (default clumping) optimised for best p-value threshold
# 2. Stacked clumping and thresholding
# 3. Best clumping and thresholding combination (No stacking)
# 4. Lassosum


library(bigsnpr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)
library(lassosum)
library(methods)
library(parallel)






########################################
## Data cleaning (ONLY PERFORM ONCE)  ##
########################################


## Remove sample_ID with missing height data
## This step only needs to be performed once!
system("unzip -j EUR.QC.zip")
base.id <- read.table("EUR.QC.fam") #Sample ID
base.id <- base.id[, c(1,2)]
colnames(base.id) <- c("FID", "IID")
tgt.data <- read.table("EUR.height", header=T) #Sample phenotype (height)
valid <- merge(tgt.data, base.id, by=c("FID", "IID"))
valid.ind <- which(base.id$FID %in% valid$FID)
## Create bigSNP file rds from bed file
snp_readBed2("EUR.QC.bed", ind.row = valid.ind)




##################################################
## Define important variables and load in data  ##
##################################################


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




######################################################
## Generate covariate data for training/testing set ##
######################################################


## Pruning and clumping of SNPs. Will output a list of SNPs that are relatively
## independent
pruned.list <- snp_clumping(G, CHR, ind.row=ind.train,
                            ncores=NCORES, exclude=snp_indLRLDR(CHR, POS))




## Perform PCA using bigsnpr SVD function and convert to PC 
## eigenvectors using prcomp().
svd <- big_randomSVD(G, snp_scaleBinom(), ncores=NCORES, ind.row=ind.train, 
                     ind.col=pruned.list)
eigenvec.train <- predict(svd, G, ind.row=ind.train, ind.col=pruned.list)[,1:6]




## Merge with sex data to form a covariate table for training set
covar.train <- data.frame(base.bigSNP$fam[ind.train, 1:2],
                          eigenvec.train, base.bigSNP$fam$sex[ind.train])
colnames(covar.train) <- c("FID", "IID","PC1", "PC2", 
                           "PC3", "PC4", "PC5", "PC6", "SEX")




## Generate principal component eigenvectors for testing data set
eigenvec.test <- predict(svd, G, ind.row=ind.test, ind.col=pruned.list)[,1:6]




## Merge with sex data to form a covariate table for testing set
covar.test <- data.frame(base.bigSNP$fam[ind.test, 1:2],
                         eigenvec.test, base.bigSNP$fam$sex[ind.test])
colnames(covar.test) <- c("FID", "IID", "PC1", "PC2", 
                          "PC3", "PC4", "PC5", "PC6", "SEX")




########################
## Standard PRS model ##
########################


## Perform clumping using standard hyper-parameter (LD r2 = 0.2,
## window size = 250kb , all p-values).
# Original summary stats file contains Odds ratio and must be converted
# to Beta score by log transformation
# Need to preserve original Plink format
sumstat.p <- read.table(gzfile("Height.QC.gz"), header=T)
sumstat.p$BETA <- log(sumstat.p$OR)
write.table(sumstat.p, "Height.QC.Transformed", quote=F, row.names=F)
# Perform clumping using plink
system(paste0("./plink \\",
              "--bfile EUR.QC \\",
              "--clump-p1 1 \\",
              "--clump-r2 0.2 \\",
              "--clump-kb 250 \\",
              "--clump Height.QC.Transformed \\",
              "--clump-snp-field SNP \\",
              "--clump-field P \\",
              "--out EUR"))
# Extract SNP ID
system("awk 'NR!=1{print $3}' EUR.clumped >  EUR.valid.snp")
# Extract p-value
system("awk '{print $3,$8}' Height.QC.Transformed > SNP.pvalue")




## Calculate PRS for separate p-value thresholds since optimum
## p-value is unknown
# Create separate files for p-value thresholds
system("echo \"0.001 0 0.001\" > range_list")
system("echo \"0.05 0 0.05\" >> range_list")
system("echo \"0.1 0 0.1\" >> range_list")
system("echo \"0.2 0 0.2\" >> range_list")
system("echo \"0.3 0 0.3\" >> range_list")
system("echo \"0.4 0 0.4\" >> range_list")
system("echo \"0.5 0 0.5\" >> range_list")




## Create txt file of training data for PRS calculation
write.table(base.bigSNP$fam[ind.train,1:2], file="train.txt", sep=" ", 
            row.names=FALSE, col.names=FALSE, quote=FALSE)




## Calculate PRS for each threshold using training dataset
system(paste0("./plink \\",
              "--bfile EUR.QC \\",
              "--keep train.txt \\",
              "--score Height.QC.Transformed 3 4 12 header \\",
              "--q-score-range range_list SNP.pvalue \\",
              "--extract EUR.valid.snp \\",
              "--out EUR"))




## Get best p-value threshold
p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
# Create training data frame for linear modelling
train.data <- merge(covar.train, pheno, by=c("FID", "IID"))
# Model using no PRS scores for baseline comparison
null.r2 <- summary(lm(Height~., data=train.data[-c(1,2)]))$r.squared
prs.result <- NULL
for (i in p.threshold) {
  score <-  read.table(paste0("EUR.",i,".profile"), header=T)[,c(1,2,6)]
  pthres.data <- merge(train.data, score, by=c("FID", "IID"))
  model <- summary(lm(Height~., pthres.data[,-c(1,2)]))
  model.r2 <- model$r.squared
  print(model.r2)
  prs.r2 <- model.r2-null.r2
  prs.coef <- model$coeff["SCORE",]
  prs.result <- rbind(prs.result, 
                      data.frame(Threshold=i, R2=prs.r2, 
                                 P=as.numeric(prs.coef[4]), 
                                 BETA=as.numeric(prs.coef[1]),
                                 SE=as.numeric(prs.coef[2])))
}
prs.result[which.max(prs.result$R2),] 
#p-value threshold of 0.4 is optimum
train.data.0.4 <- merge(train.data, 
                        read.table("EUR.0.4.profile", header=TRUE)[,c(1,2,6)],
                        by=c("FID", "IID"))




## Get PRS score using optimum p-value threshold
# Create text file of test ID
write.table(base.bigSNP$fam[ind.test,1:2], file="test.txt", sep=" ", 
            row.names=FALSE, col.names=FALSE, quote=FALSE)
# Calculate PRS using optimal p-value threshold of 0.4
system("echo \"0.4 0 0.4\" >> test_list")
system(paste0("./plink \\",
              "--bfile EUR.QC \\",
              "--keep test.txt \\",
              "--score Height.QC.Transformed 3 4 12 header \\",
              "--q-score-range test_list SNP.pvalue \\",
              "--extract EUR.valid.snp \\",
              "--out EUR.test"))
# PRS stored on EUR.test.0.4.profile
test.data.0.4 <- merge(covar.test, 
                       read.table("EUR.test.0.4.profile", header=T)[,c(1,2,6)],
                       by=c("FID", "IID"))




# Predict using optimum training linear model
train.prs <- lm(Height~., data=train.data.0.4[,-c(1,2)])
test.pred.prs <- predict(train.prs, test.data.0.4[,-c(1,2)])




## Calculate coefficient of determiniation for SCT model
SSR.PRS <- sum((pheno$Height[ind.test] - test.pred.prs)^2)
SST <- sum((pheno$Height[ind.test] - mean(pheno$Height[ind.test]))^2)
R2.PRS <- 1 - (SSR.PRS / SST)




#####################
## Stacking model ##
####################


## Create BETA and p-value (log) vectors
beta <- rep(NA, ncol(G))
beta[info_snp$`_NUM_ID_`] <- info_snp$beta
lpval <- rep(NA, ncol(G))
lpval[info_snp$`_NUM_ID_`] <- -log10(info_snp$p)



## LONGEST STEP
## Perform clumping
# 22 chromosomes x 28 clumping parameter combinations
all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train,
                              lpS = lpval, exclude = which(is.na(lpval)),
                              ncores = NCORES)
attr(all_keep, "grid")




## Perform threshold and PRS calculations using threshold/clumping sets
#backing file
loc <- paste0(getwd(),"/tmp")
multi_PRS <- snp_grid_PRS(G, all_keep, beta, lpval, ind.row=ind.train,
                          backingfile=loc, n_thr_lpS = 50, ncores = NCORES)

# 22 CHR x 28 clumping combinations x 50 p-val thresholds = 30800 PRS
# multi_PRS = 330 x 30800




## Perform linear regression to learn PRS weights
## Covariate penalisation is not performed (pf.covar)
final.mod.train <- snp_grid_stacking(multi_PRS, 
                                     pheno$Height[ind.train],
                                     ncores=NCORES, K=4,
                                     covar.train=covar_from_df(covar.train[,-c(1,2)]),
                                     pf.covar=c(rep(0, 7)))
#Update with new beta scores from derived weights
train.beta <- final.mod.train$beta.G #New beta scores
train.beta.covar <- final.mod.train$beta.covar #Covariate weightings
train.ind <- which(train.beta != 0) #Skip any beta scores with value of zero




## Visualise BETA score changes
ggplot(data.frame(y = train.beta, x = beta)[ind.train, ]) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_abline(slope = 0, intercept = 0, color = "blue") +
  geom_point(aes(x, y), size = 0.6) +
  theme_bigstatsr() +
  labs(x = "Effect sizes from GWAS", y = "Non-zero effect sizes from SCT")




## Use fitted model from training to predict test phenotype
# First create FBM of combined testing genotype and covariates
# Exclude training rows and genotypes with weighting of zero
test.data.sct <- as_FBM(cbind(G[ind.test, train.ind], covar.test[,-c(1,2)]))
test.pred.sct <- final.mod.train$intercept + big_prodVec(test.data.sct,
                                                         c(train.beta[train.ind],
                                                           train.beta.covar))




## Calculate coefficient of determiniation for SCT model
SSR.SCT <- sum((pheno$Height[ind.test] - test.pred.sct)^2)
R2.SCT <- 1 - (SSR.SCT / SST)




################################################
## Best clumping and thresholding combination ##
################################################

# Note that there is no stacking involved

#multi_PRS = 330 rows x 30800 columns
#          = 330 rows x [22 chr x 1400 (C+T combination)]




## Convert to tidyverse data frame and add two columns: lpval and id.number
grid2 <- attr(all_keep, "grid") %>%
  mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), id = row_number()) %>%
  tidyr::unnest(cols = "thr.lp")
s <- nrow(grid2)




## Get R2 performance of each 1400 C+T combination and add to grid
grid2$r2 <- big_apply(multi_PRS, a.FUN=function(X, ind, s, train){
  #Sum PRS score for one C+T combination over the 22 chromosomes
  single_PRS <- rowSums(X[,ind + s * (0:21)])
  dat <- cbind(train, PRS=single_PRS) 
  summary(lm(Height~., data=dat[,-c(1,2)]))$r.squared
}, ind=1:s, s=s, train=train.data,  
a.combine = 'c', block.size = 1, ncores = NCORES)




## Get SNP variants in top C+T combination
# Print top 10 C+T combinations for comparison and isolate top C+T
max_prs <- grid2 %>% arrange(desc(r2)) %>% slice(1:10) %>% print() %>% slice(1)
# Using the top C+T combination, get a list of the variants in that group
max_snp <- unlist(purrr::map(all_keep, max_prs$id))




## Create linear regression model from training
# Get PRS scores using optimal C+T parameters
maxCT.prs.train <- snp_PRS(G, beta[max_snp], ind.test=ind.train, 
                           ind.keep=max_snp,lpS.keep=lpval[max_snp], 
                           thr.list=max_prs$thr.lp)
# Combine with covariance data
maxCT.train <- cbind(train.data, maxCT.prs.train)
colnames(maxCT.train)[11] <- "PRS"
maxCT.model <- lm(Height~., data=maxCT.train[,-c(1,2)])


## Generate PRS score using optimal C+T 
maxCT.prs.test <- snp_PRS(G, beta[max_snp], ind.test=ind.test, ind.keep=max_snp,
                          lpS.keep=lpval[max_snp], thr.list=max_prs$thr.lp)
maxCT.test <- cbind(covar.test, maxCT.prs.test)
colnames(maxCT.test)[10] <- "PRS"
maxCT.pred <- predict(maxCT.model, maxCT.test[, -c(1,2)])




## Calculate coefficient of determiniation for SCT model
SSR.maxCT <- sum((pheno$Height[ind.test] - maxCT.pred)^2)
R2.maxCT <- 1 - (SSR.maxCT / SST)





####################
## Lassosum model ##
####################


## Define important variables
# Two threads for multi-thread processing
cl <- makeCluster(2)
# Bed/Bim/Fam file prefix
bfile <- "EUR.QC"
# Use LD regions defined in Berisa and Pickrell (2015) for the European 
# population and the hg19 genome.
ld.file <- "EUR.hg19"




## Summary statistic processing
# Read in summary statistics
ss <- read.table("Height.QC.gz", header=T)
# Remove p-values that equal zero
which(ss$p == 0) #There are none
# Convert p-values to correlation via t-statistic
cor <- p2cor(p=ss$P, n=ss$N, sign=log(ss$OR))




## Lassosum pipeline
out <- lassosum.pipeline(
  cor = cor,
  chr = ss$CHR,
  pos = ss$BP,
  A1 = ss$A1,
  A2 = ss$A2,
  ref.bfile=bfile,
  keep.ref=as.data.frame(covar.train[,c(1,2)]), #Train data ID
  test.bfile=bfile,
  keep.test=as.data.frame(covar.train[,c(1,2)]), #Train data ID
  LDblocks = ld.file, 
  cluster=cl
)




## Training performance model
v.train <- validate(out, 
                    pheno=as.data.frame(pheno[ind.train,]), 
                    covar=as.data.frame(covar.train))
max(v.train$validation.table$value)^2




## Perform on test data
# Get new beta scores
out2 <- subset(out, s=v.train$best.s, lambda=v.train$best.lambda)
v.test <- validate(out2, 
                   covar=as.data.frame(covar.test),
                   pheno=as.data.frame(pheno[ind.test,]),
                   test.bfile=bfile,
                   keep=as.data.frame(covar.test[,c(1,2)]))
R2.Lasso <- v.test$validation.table$value^2




############################
## Final table of results ##
############################

results <- data.frame("PRS Model"=c("Standard PRS","Stacked C+T", 
                                    "Max C+T","Lassosum"),
                      "R2" = c(R2.PRS, R2.SCT, R2.maxCT, R2.Lasso))
print(results)



