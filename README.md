# methylation_array_processing
These scripts perform DNA methylation array pre-processing - filtering & normalization
Follows the first steps outlined in Wilhelm-Benartzi et al. 2013 (see Figure 1 and Table 1)

I recommend creating a new environment (i.e. using conda) to install all required software.
First install R version 3.5.1 (this older version avoids package conflicts), then install the following:

CRAN R packages:
RPMM

Bioconductor R packages:
minfi
wateRmelon
IlluminaHumanMethylation450kmanifest
IlluminaHumanMethylationEPICmanifest
IlluminaHumanMethylationEPICanno.ilm10b4.hg19
FlowSorted.Blood.450k
FlowSorted.Blood.EPIC

-- RUNNING THE PIPELINE --
1. Copy all scripts to a working directory
2. Specify all desired options for the R scripts within 'process_methylation_data.sh', including adding paths to .idats, sample sheets, etc.
3. Set up and activate the working environment as described above
4. Execute 'process_methylation_data.sh'
