#!/usr/bin/env Rscript
## Step 3 of methylation pre-processing - normalization

## Set-up ----
library("argparser")

p <- arg_parser("Perform normalization and optional per-sample blood cell proportion/age estimation. Requires R packages 'minfi, 'wateRmelon', and a manifest package e.g. 'IlluminaHumanMethylationEPICmanifest' or 'IlluminaHumanMethylation450kmanifest'.")
p <- add_argument(p, short="-i", "--input", help="Required argument - Path to .Rdata object containing output of prior processing step (2_probeQC.R).", nargs=1)
p <- add_argument(p, short="-W", "--estimate_WBC_props", help="Should blood cell type proportions be estimated? Default=FALSE", type="logical", default=FALSE, nargs=1)
p <- add_argument(p, short="-W", "--estimate_age", help="Should individual age be estimated? Default=FALSE", type="logical", default=FALSE, nargs=1)
p <- add_argument(p, short="-s", "--set_poor_CpGs_to_NA", help="Should beta values that fail the detection p-value threshold be set to NA? Default=TRUE", type="logical", default=TRUE, nargs=1)
p <- add_argument(p, short="-o", "--output", help="Prefix for output files.", default="3_normalization_output", nargs=1)

args <- parse_args(p)

input <- args$input
predictWBC <- args$estimate_WBC_props
predict.age <- args$estimate_age
set.fail.betas.to.NA <- args$set_poor_CpGs_to_NA
output <- args$output

if (is.na(input)) {
  stop("Input argument is missing")
}

library(minfi)
library(wateRmelon)

## Load data ----
cat("loading .RData input... ")
load(input)
cat("done\n")

## Estimate WBC composition (only interpret for blood samples) ----
ARRAY <- annotation(Dat.sample.probesQCed)[[1]]
if (predictWBC == T) {
  cat("estimating blood cell proportions...\n")
  # EPIC-specific version of this function (estimateCellCounts2) broken as of 03/19/2021
  cellProps <- estimateCellCounts.wmln(as.methylumi(preprocessRaw(Dat.sample.probesQCed)), referencePlatform=ARRAY)
  phenos.df <- cbind(phenos.df, cellProps[match(phenos.df$Sample_Name, rownames(cellProps)),])
  cat("done\n")
}

## Normalization & final QC ----
# Colour and background correction with noob, followed by functional normalization (between array method)
cat("starting functional normalization...\n")
Dat.funnorm <- preprocessFunnorm(Dat.sample.probesQCed)
cat("done\n")
# Type I and II probe scaling w/ BMIQ (within array method)
probeType <- getProbeType(Dat.funnorm);probeType[probeType=="I"] <- 1; probeType[probeType=="II"] <- 2
cat("starting BMIQ normalization...\n")
  library(parallel)
  CL <- makeCluster(4, type="FORK")
    Dat.funnorm.BMIQ <- parApply(CL, getBeta(Dat.funnorm), 2, FUN=BMIQ, design.v=probeType, plots=F)
  stopCluster(CL)
unlist.BMIQ <- function(dat) {
  betas <- dat$nbeta
  return(betas)
}
Dat.funnorm.BMIQ.df <- as.data.frame(sapply(Dat.funnorm.BMIQ, unlist.BMIQ, simplify="array"))
cat("done\n")
# Setting beta values that do not pass p-value threshold to NA
if (set.fail.betas.to.NA == T) {
  cat("setting poor beta-value measures to NA... ")
  detP.sorted <- detP[match(rownames(Dat.funnorm.BMIQ.df), rownames(detP)), match(colnames(Dat.funnorm.BMIQ.df), colnames(detP))]
  d.pval.fails <- which(detP.sorted > pval.thresh, arr.ind=T)
  for (ind in unique(d.pval.fails[,2])) {
    Dat.funnorm.BMIQ.df[names(which(d.pval.fails[,2]==ind)), ind] <- NA
  }
  cat("done\n")
}

## Predict age ----
if (predict.age == T) {
  cat("estimating epigenetic age...\n")
  predictedAge <- agep(Dat.funnorm.BMIQ.df)
  phenos.df$Horvath_DNAmAge <- predictedAge[match(phenos.df$Sample_Name, rownames(predictedAge)),]
  cat("done\n")
}

## Finishing up ----
cat("finished probe QC, saving RData object and writing out .csv files\n")
save(list=c("Dat.funnorm.BMIQ", "phenos.df"), file=paste(output,"RData",sep="."))
write.csv(phenos.df, paste(output,"_sampleInfo.csv",sep=""), quote=F, row.names=F)
write.csv(Dat.funnorm.BMIQ.df, paste(output,"_normBetas.csv",sep=""), quote=F)
