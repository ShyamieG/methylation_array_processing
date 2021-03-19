#!/usr/bin/env Rscript
## Step 2 of methylation pre-processing - probe QC

## Set-up ----
library("argparser")

p <- arg_parser("Perform probe QC. Requires R package 'minfi', and a manifest package e.g. 'IlluminaHumanMethylationEPICmanifest' or 'IlluminaHumanMethylation450kmanifest'.")
p <- add_argument(p, short="-i", "--input", help="Required argument - Path to .Rdata object containing output of prior processing step (1_sampleQC.R).", nargs=1)
p <- add_argument(p, short="-c", "--cross_reactive_probes", help="Required argument - Path to .RData object containing list of array-specific cross-reactive probes.", nargs=1)
p <- add_argument(p, short="-s", "--exclude_sex_chrs", help="Should probes on X and Y chromosomes be removed? Default = TRUE", type="logical", default=TRUE, nargs=1)
p <- add_argument(p, short="-P", "--probe_p_val_cutoff", help="Maximum proportion of a probe's CpGs that can fail detection p-value criteria and still be retained.", default=0.05, type="double", nargs=1)
p <- add_argument(p, short="-o", "--output", help="Prefix for RData output file.", default="2_probeQC_output", nargs=1)

args <- parse_args(p)

input <- args$input
cr.probes <- args$cross_reactive_probes
exclude.sex.chrs <- args$exclude_sex_chrs
probe.pval.cutoff <- args$probe_p_val_cutoff
output <- args$output

if (is.na(input) | is.na(cr.probes)) {
  stop("At least one required argument is missing")
}

library(minfi)

## Load data ----
cat("loading .RData input... ")
load(input)
cat("done\n")

## Filter CpG sites ----
probes.to.remove <- c()
# Extract sex chromosome probes
if (exclude.sex.chrs==T) {
  cat("flagging sex chromosome probes... ")
  XY.probes <- which(getAnnotation(Dat.sampleQCed)$chr %in% c("chrX", "chrY"))
  cat("done\n")
  probes.to.remove <- c(probes.to.remove, XY.probes)
}
# Extract CpGs with a detection P value above the threshold
detection.fail <- function(pvals, pval.thresh) {
    return(length(which(as.numeric(pvals) > pval.thresh)))
}
cat("flagging probes with poor detection p-values... ")
detP.fail.probes <- names(which(apply(detP, MARGIN=1, FUN=detection.fail, pval.thresh=pval.thresh) > floor(length(sampleNames(Dat.sampleQCed))*probe.pval.cutoff)))
probes.to.remove <- c(probes.to.remove, detP.fail.probes)
cat("done\n")
# Load published list of cross-reactive probes
cat("loading cross-reactive probes... ")
load(cr.probes)
probes.to.remove <- c(probes.to.remove, CR.probes)
cat("done\n")
# Filter probes
cat("filtering probes from data... ")
Dat.sample.probesQCed <- subsetByLoci(Dat.sampleQCed, excludeLoci=probes.to.remove, keepControls=T, keepSnps=F)
cat("done\n")

## Finishing up ----
cat("finished probe QC, saving RData object\n")
save(list=c("Dat.sample.probesQCed", "XY.probes", "detP.fail.probes", "CR.probes", "phenos.df", "detP", "pval.thresh", "drop_sex_mismatches"), file=paste(output,"RData",sep="."))
