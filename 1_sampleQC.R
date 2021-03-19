#!/usr/bin/env Rscript
## Step 1 of methylation pre-processing - sample QC

## Set-up ----
library("argparser")
options(stringsAsFactors=F)

p <- arg_parser("Perform sample QC. Requires R package 'minfi', and a manifest package e.g. 'IlluminaHumanMethylationEPICmanifest' or 'IlluminaHumanMethylation450kmanifest'.")
p <- add_argument(p, short="-I", "--idat_dir", help="Required argument - Path to directories containing .idat files. Directories can be nested.", nargs=1)
p <- add_argument(p, short="-D", "--sample_sheet_dir", help="Required argument - Path to directory containing sample sheet from U Miami.", nargs=1)
p <- add_argument(p, short="-S", "--sample_sheet_csv", help="Required argument - Name of sample sheet (.csv format) provided by U Miami. Must contain sample IDs in a column called 'Sample_Names'. Must also contain 'Slide' and 'Array' columns", nargs=1)
p <- add_argument(p, short="-M", "--sample_manifest", help="Required argument - Path to sample manifest (.csv format) provided to U Miami. Must have columns called 'sample_external_id', 'gender', and 'affection_status'.", nargs=1)
p <- add_argument(p, short="-d", "--drop_sex_mismatches", help="Should samples with non-matching estimated sex and reported gender be dropped? Default = TRUE", default=TRUE, type="logical", nargs=1)
p <- add_argument(p, short="-m", "--median_cutoff", help="Minimum value of median signal (methylated or unmethylated) for a sample to pass.", default=10.5, type="double", nargs=1)
p <- add_argument(p, short="-p", "--CpG_p_val_threshold", help="Maximum detection p-value for a given measurement to be considered reliable.", default=0.05, type="double", nargs=1)
p <- add_argument(p, short="-P", "--sample_p_val_cutoff", help="Maximum proportion of a sample's CpGs that can fail detection p-value criteria and still be retained.", default=0.1, type="double", nargs=1)
p <- add_argument(p, short="-o", "--output", help="Prefix for RData output file.", default="1_sampleQC_output", nargs=1)

args <- parse_args(p)

IDAT.dir <- args$idat_dir
sample.sheet.dir <- args$sample_sheet_dir
sample.sheet.file <- args$sample_sheet_csv
sample.manifest <- args$sample_manifest
drop_sex_mismatches <- args$drop_sex_mismatches
median.cutoff <- args$median_cutoff
pval.thresh <- args$CpG_p_val_threshold
sample.pval.cutoff <- args$sample_p_val_cutoff
output <- args$output

if (is.na(IDAT.dir) | is.na(sample.sheet.dir) | is.na(sample.sheet.file) | is.na(sample.manifest)) {
  stop("At least one required argument is missing")
}

library("minfi")

## Load data ----
cat("loading idats... ")
Dat <- read.metharray.exp(base=IDAT.dir, recursive=T, force=T)
cat("done\n")
sample.sheet <- read.metharray.sheet(base=sample.sheet.dir, pattern=sample.sheet.file)
# Rename samples
sampleNames(Dat) <- sample.sheet[match(sampleNames(Dat), paste(sample.sheet$Slide, sample.sheet$Array, sep="_")), "Sample_Name"]
manifest <- read.csv(sample.manifest, header=T)

## Flag sex mismatches ----
cat("calculating quality control metrics... ")
QCinfo <- minfiQC(preprocessRaw(Dat))
cat("done\n")
phenos.df <- sample.sheet[match(rownames(QCinfo$qc), sample.sheet[, "Sample_Name"]),c("Sample_Name","Slide","Array")]
phenos.df <- as.data.frame(cbind(phenos.df, QCinfo$qc[,c("mMed", "uMed", "predictedSex")]))
# Identify mismatches based on sex
phenos.df <- cbind(phenos.df, manifest[match(phenos.df$Sample_Name, manifest$hihg_sample_id),c("sample_external_id", "gender", "affection_status")])
sex.mismatches <- phenos.df[phenos.df[,"gender"] != phenos.df[,"predictedSex"], "Sample_Name"]
cat(paste(length(sex.mismatches), " samples with inconsistent sex\n"))
phenos.df$sex_mismatch <- "PASS"
if (length(sex.mismatches)>0) {
  phenos.df$sex_mismatch[match(sex.mismatches, phenos.df$Sample_Name)] <-"FAIL"
}

## Technical QC ----
# Extract detection p-values
cat("extracting detection p-values... ")
detP <- detectionP(Dat)
cat("done\n")
# Identify individuals who fail threshold
detection.fail <- function(pvals, pval.thresh) {
  return(length(which(as.numeric(pvals) > pval.thresh)))
}
# Flag samples with low median signal (either U or M) or a high proportion (10%+) of CpGs that fail the detection p-value filter
per.sample.pval.fails <- apply(detP, 2, detection.fail, pval.thresh=pval.thresh)
pval.failed.samples <- sampleNames(Dat)[which(per.sample.pval.fails > floor(nrow(detP)*sample.pval.cutoff))]
cat(paste(length(pval.failed.samples), " samples failed p-value cutoff\n"))
signal.failed.samples <- sampleNames(Dat)[which(phenos.df$mMed < median.cutoff | phenos.df$uMed < median.cutoff)]
cat(paste(length(signal.failed.samples), " samples failed signal cutoff\n"))
# Filter out all samples that fail QC (sex check, median intensity check)
if (drop_sex_mismatches==TRUE) {
  bad.samples <- unique(c(sex.mismatches, pval.failed.samples, signal.failed.samples))
} else {
  bad.samples <- unique(c(pval.failed.samples, signal.failed.samples))
}
if (length(bad.samples) > 0) {
    Dat.sampleQCed <- Dat[,-match(bad.samples, colnames(Dat))]
} else {
    Dat.sampleQCed <- Dat
}
# Annotated samples filtered out in phenotype dataframe
phenos.df$detectionP <- "PASS"
phenos.df$medianSignal <- "PASS"
if (length(pval.failed.samples)>0) {
  phenos.df$detectionP[match(pval.failed.samples, phenos.df$Sample_Name)] <-"FAIL"
}
if (length(signal.failed.samples)>0) {
  phenos.df$medianSignal[match(signal.failed.samples, phenos.df$Sample_Name)] <-"FAIL"
}

## Finishing up ----
cat("finished sample QC, saving RData object\n")
save(list=c("Dat.sampleQCed", "detP", "phenos.df", "pval.thresh", "drop_sex_mismatches"), file=paste(output,"RData",sep="."))
