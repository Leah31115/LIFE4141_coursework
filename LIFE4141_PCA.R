rm(list=ls()) # Clear my environment

# Load in libraries
library(ggplot2)
library(adegenet)
library(vcfR)
library(pegas)
library(StAMPP)
library(ade4)
library(adegraphics)
library(tidyverse)
library(qqman)
library(ape)
library(ggrepel)


# This script will load vcf files, perform PCA, and make a Neighbour-joining tree

# -----------------------------------------------------------------------------
# Full credit goes to Filip Kolar 2017, Sian Bray 2018 and Levi Yant 2022, for this code to
# Load vcf, make genlight objects and run PCA

# Set working directory 
setwd("C:/Users/leahe/Documents/LIFE4141_Comparative_and_Evolutionary_Genomics/LIFE4141_coursework_resources/evo_data_files")

# Import SNP data from VCF
vcf <- read.vcfR("LAB_NEN_ODN.clean_BI.ann.3mbChr5.vcf.gz")

# convert to genlight 	
aa.genlight <- vcfR2genlight(vcf) # vcfR2genlight for diploid organisms is used instead of the vcfR2genlight.tetra modified function
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")   # add real SNP.names
pop(aa.genlight) <- substr(indNames(aa.genlight),1,3)             # add pop names: here pop names are first 3 chars of ind name

# Check
aa.genlight
indNames(aa.genlight)
ploidy(aa.genlight)

# Fast pca function
glPcaFast <- function(x,
                      center=TRUE,
                      scale=FALSE,
                      nf=NULL,
                      loadings=TRUE,
                      alleleAsUnit=FALSE,
                      returnDotProd=FALSE){
  
  if(!inherits(x, "genlight")) stop("x is not a genlight object")
  # keep the original mean / var code, as it's used further down
  # and has some NA checks..
  if(center) {
    vecMeans <- glMean(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }
  if(scale){
    vecVar <- glVar(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }
  # convert to full data, try to keep the NA handling as similar
  # to the original as possible
  # - dividing by ploidy keeps the NAs
  mx <- t(sapply(x$gen, as.integer)) / ploidy(x)
  # handle NAs
  NAidx <- which(is.na(mx), arr.ind = T)
  if (center) {
    mx[NAidx] <- vecMeans[NAidx[,2]]
  } else {
    mx[NAidx] <- 0
  }
  # center and scale
  mx <- scale(mx,
              center = if (center) vecMeans else F,
              scale = if (scale) vecVar else F)
  # all dot products at once using underlying BLAS
  # to support thousands of samples, this could be
  # replaced by 'Truncated SVD', but it would require more changes
  # in the code around
  allProd <- tcrossprod(mx) / nInd(x) # assume uniform weights
  ## PERFORM THE ANALYSIS ##
  ## eigenanalysis
  eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]
  ## scan nb of axes retained
  if(is.null(nf)){
    barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n = 1))
  }
  ## rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  ##res$matprod <- allProd # for debugging
  ## use: li = XQU = V\Lambda^(1/2)
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")
  ## GET LOADINGS ##
  ## need to decompose X^TDV into a sum of n matrices of dim p*r
  ## but only two such matrices are represented at a time
  if(loadings){
    if(scale) {
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # create empty matrix
    ## use: c1 = X^TDV
    ## and X^TV = A_1 + ... + A_n
    ## with A_k = X_[k-]^T v[k-]
    myPloidy <- ploidy(x)
    for(k in 1:nInd(x)){
      temp <- as.integer(x@gen[[k]]) / myPloidy[k]
      if(center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.na(temp)] <- 0
      }
      if(scale){
        temp <- temp/vecSd
      }
      res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
    }
    res$loadings <- res$loadings / nInd(x) # don't forget the /n of X_tDV
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }
  ## FORMAT OUTPUT ##
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }
  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
      rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }
  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }
  res$call <- match.call()
  class(res) <- "glPca"
  return(res)
}


# Run PCA
pca.1 <- glPcaFast(aa.genlight, nf=300)

# Loading plot
loadingplot(pca.1)

# proportion of explained variance by first three axes

pca.1$eig[1]/sum(pca.1$eig) # proportion of variation explained by 1st axis
pca.1$eig[2]/sum(pca.1$eig) # proportion of variation explained by 2nd axis
pca.1$eig[3]/sum(pca.1$eig) # proportion of variation explained by 3rd axis


# Colour pops
col <- funky(10)

# PCA plot
pdf("2024_12_26_PCA.pdf", width = 14, height = 7)
g1 <- s.class(pca.1$scores, pop(aa.genlight), xax = 1, yax = 2, 
              col = transp(col,.6), ellipseSize = 0, starSize = 0, 
              ppoints.cex = 4, paxes.draw = T, pgrid.draw = F, plot = FALSE, 
              xlab = "PC1 (23.3%)", ylab = "PC2 (7.8%)")

g2 <- s.label(pca.1$scores, xax = 1, yax = 2, ppoints.col = "red", 
              plabels = list(box = list(draw = FALSE), optim = TRUE), 
              paxes.draw = T, pgrid.draw = F, plabels.cex = 1, plot = FALSE, 
              xlab = "PC1 (23.3%)", ylab = "PC2 (7.8%)")

ADEgS(c(g1, g2), layout = c(1, 2))
dev.off()

# -----------------------------------------------------------------------------

# PCA with points coloured by soil Zinc levels
sample_info <- read.csv("metal_levels.csv")

# Extract PC1 and PC2 values from pca.1
pca_values <- data.frame(pca.1[["scores"]]) %>%
  select("PC1", "PC2")

# Make index column a data column
pca_values <- rownames_to_column(pca_values)
colnames(pca_values)[1] ="Individual"

# Merge pc1 and sample_info
pca_sample_data <- merge(pca_values, sample_info, by = "Individual")

# Plot PCA and colour point by Zn levels
ggplot(pca_sample_data, aes(PC1, PC2, colour = Zn)) +
  geom_point() +
  geom_text_repel(aes(label=Individual), size=3, force=20, color="black", max.overlaps = 20) +
  xlab("PC1 (23.3%)") +
  ylab("PC2 (7.8%)") +
  theme_bw()


# Make a neighbour-joining tree using the genlight objects
nj_tree <- nj(dist(as.matrix(aa.genlight)))

# phylogram
plot(nj_tree, typ = "phylogram", cex = 0.8)
title("Neighbour-joining tree of the Cochlearia pyrenaica \n individuals")


