###################################################
#                                                 #
#      eQTL/pQTL & Mediation Examples             #
#          Updated June 15, 2016                  #
#              Steve Munger                       #
#           steven.munger@jax.org                 #
#                                                 #
###################################################

#### This script runs best in RStudio (especially the interactive plots)  www.rstudio.com
#### Please forgive my clunky code. You can teach an old dog new tricks, but they won't necessarily be pretty tricks.

### First, install some R packages and our data
options(stringsAsFactors = F) #This will save headaches when working with data frames that contain character strings...
install.packages("devtools")
library(devtools)

install_github("dmgatti/DOQTL")
install_github("simecek/intermediate")
install_github("kbroman/qtlcharts")

library(DOQTL)
library(intermediate)
library(qtlcharts)

setwd("~/") #Set working directory
load("~/ChickMungeretal2016_DiversityOutbred.Rdata")  ###Load in the dataset

# Data objects include:

#### RNA data
# expr.rna.192 <- 192 samples (rows) x 21,454 genes (columns) Liver RNA-seq expression data. Data was upper quartile normalized and transformed to rank normal scores.
# annotations.rna.192 <- Gene information for the 21,454 genes. 21,454 rows (genes) x 6 columns (annotations). Order of rows corresponds to order of columns in expr.rna.192.
# covariates.rna.192 <- Experimental variable information for the 192 Diversity Outbred samples.
# X.rna <- Covariate matrix for eQTL mapping (includes additive effects and interaction of sex, diet, and batch)

#### Protein data
# expr.protein.192 <- 192 samples (rows) x 8,050 proteins (columns) Liver protein abundance data. Data was quantile normalized and transformed to rank normal scores.
# annotations.protein.192 <- Protein information for 8,050 proteins. 8,050 rows (proteins) x 11 columns (annotations). Order of rows corresponds to order of columns in expr.protein.192.
# covariates.protein.192 <- Experimental variable information (pertinent to protein dataset) for the 192 Diversity Outbred samples.
# X.protein <- Covariate matrix for pQTL mapping (includes additive effects and interaction of sex and diet)

#### Common objects used in e/pQTL mapping and mediation analysis
# K.LOCO.192 <- Kinship Matrix for QTL mapping. Constructed using "Leave One Chromosome Out (LOCO)" method. List with 20 elements corresponding to the 19 autosomes and X chromosome.
# probs.192 <- 8-state founder genotypes at 64,000 imputed markers for each of 192 samples. 192x8x64000 array.
# markers.64k <- Location information for the 64,000 imputed markers.
# samples.192 <- List of 192 Diversity Outbred sample names.



### SCAN for eQTL and pQTL
## Example gene = Tmem68

##################################################################################################
## 1) Scan for eQTL

my.gene <- "Tmem68"  ### Input your gene of interest

target.rna.index <- which(annotations.rna.192$Gene == my.gene) ### find the index number for your gene of interest
annotations.rna.192[target.rna.index,]  ### Show information for the gene of interest
scanone.rna <- scanone(expr.rna.192, pheno.col=target.rna.index,K=K.LOCO.192, probs=probs.192, snps=markers.64K, addcovar=X.rna[,-1]) ## Perform the eQTL scan
plot(scanone.rna) ## plot the eQTL LOD scores


# A little function to find the SNP maximizing LOD score (autosomes only) 
argmax.lod <- function(scanone.fit)
  scanone.fit$lod$A$marker[which.max(scanone.fit$lod$A$lod)[1]]


# Let's plot the founder strain coefficients for the autosome with the peak LOD marker
argmax.snp.rna <- argmax.lod(scanone.rna)
coefplot(scanone.rna, chr=markers.64K[argmax.snp.rna,"chr"])


coefplot(scanone.rna,chr=4)  # Run this if you want to manually input the chromosome to plot

###################################################################################################





###################################################################################################
## 2) Scan for pQTL
target.protein.index <- which(annotations.protein.192$Associated.Gene.Name == my.gene)
scanone.protein <- scanone(expr.protein.192, pheno.col=target.protein.index, K=K.LOCO.192,probs=probs.192, snps=markers.64K, addcovar=X.protein[,-1])
plot(scanone.protein)
# 
# effect plot for autosome with max. LOD
argmax.snp.protein <- argmax.lod(scanone.protein)
coefplot(scanone.protein, chr=markers.64K[argmax.snp.protein,"chr"])

coefplot(scanone.protein, chr=13)  # Run this if you want to manually input the chromosome to plot
###################################################################################################





###################################################################################################
# Mediation Scan

#Requires:
#### target - numeric vector with transcript/protein expression
#### mediator - matrix, each column is one transcript/protein's expression
#### annotation - data.frame with mediator annotation, must include columns "chr" and "pos"
#### qtl.geno - matrix, haplotype probabilities at the QTL we want to mediate
#### covar - additive covariates
#### method = c("ignore", "lod-diff", "double-lod-diff", "lod-ratio")  ### we prefer "double-lod-diff"


## 3) Mediation Scan - Condition distant pQTL on protein intermediates
y <- expr.protein.192[,target.protein.index]
geno.argmax.protein <- probs.192[,-1,argmax.snp.protein]

# trim annotation, calculate middle point
annot.protein <- annotations.protein.192[,c("Ensembl.Protein.ID", "Ensembl.Gene.ID", "Associated.Gene.Name")]
annot.protein$Chr <- annotations.protein.192$Chromosome.Name
annot.protein$Pos <- (annotations.protein.192$Gene.Start..bp. + annotations.protein.192$Gene.End..bp.)/2

med <- mediation.scan(target=y, mediator=expr.protein.192, annotation=annot.protein, 
                            covar=X.protein[,-1], qtl.geno=geno.argmax.protein,method="double-lod-diff")
kplot(med) #Interactive Plot - Hover over points to see gene symbols
plot(med)  #Static plot



## 4) Mediation Scan - Condition distant pQTL on transcript intermediates

# trim annotation, calculate middle point
annot.rna <- annotations.rna.192[,c("EnsemblID", "Gene", "Chr")]
colnames(annot.rna) = c("Ensemble.Gene.ID","Associated.Gene.Name","Chr")
annot.rna$Pos <- (annotations.rna.192$Start.Mbp + annotations.rna.192$End.Mbp)/2

med <- mediation.scan(target=y, mediator=expr.rna.192, annotation=annot.rna, 
                            covar=X.protein[,-1], qtl.geno=geno.argmax.protein)

kplot(med)  #Interactive Plot - Hover over points to see gene symbols
plot(med)   #Static plot

######################################################################################################



#### Other example proteins to scan

##  Ndufaf1
##  Mtr
##  Cct7
##  Glul
##  Xrcc6
##  Elp3
##  Aven
##  Klc4




###################################################################
### Optional code to set different covariates for RNA and Protein
# X.rna <- model.matrix(~Sex*Diet*Batch, covariates.rna.192)
# colnames(X.rna)[2] <- "sex" # DOQTL requirement

# X.protein <- model.matrix(~Sex*Diet, covariates.protein.192)
# colnames(X.protein)[2] <- "sex" # DOQTL requirement
###################################################################








