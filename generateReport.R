# R script to generate dMLPA reports from filtered excel files
#
# ENSURE PACKAGES ARE INSTALLED
#
options(menu.graphics = FALSE)
options(repos = structure(c(CRAN="http://cran.ma.imperial.ac.uk/")))
if (!("knitr" %in% installed.packages())) install.packages("knitr",dependencies=TRUE)
if (!("kableExtra" %in% installed.packages())) install.packages("kableExtra",dependencies=TRUE)


#
# LOAD LIBRARIES
#
cat(paste('Loading Packages\n'))
require(warn.conflicts=FALSE,quietly=FALSE,package="GenomicRanges")
require(warn.conflicts=FALSE,quietly=FALSE,package="readxl")
require(warn.conflicts=FALSE,quietly=FALSE,package="xtable")
require(warn.conflicts=FALSE,quietly=FALSE,package="knitr")
require(warn.conflicts=FALSE,quietly=FALSE,package="kableExtra")
require(warn.conflicts=FALSE,quietly=FALSE,package="Rsamtools")
require(warn.conflicts=FALSE,quietly=FALSE,package="stringr")
require(warn.conflicts=FALSE,quietly=FALSE,package="dplyr")
source('ed2vcf.R')

#
# READ ARGS
#
cmd<-commandArgs(trailingOnly = FALSE)
runningScript<-unlist(strsplit(cmd[which(substr(cmd,1,7)=='--file=')],'='))[2]
scriptDirectory<-normalizePath(dirname(runningScript))

args<-commandArgs(trailingOnly = TRUE)
infile<-args[1]
outfile<-args[2]
annfile<-args[3]


# define version string (this is imported from the docker image if available)
versionfile <- paste(scriptDirectory, "VERSION", sep="/")
versionstr <- ifelse(file.exists(versionfile), readChar(versionfile, 7), "DEVELOP")

#
# Read input data
#
if (!file.exists(infile)) stop(paste("Input file", infile, "does not exist"))
metadata <- read_excel(infile, sheet = "Sample info")
results <- read_excel(infile, sheet = "Gene filtered results")
samplename <- colnames(metadata)[2]

#
# load annotations
#
annotations<-NA
if (file.exists(annfile)) {
  message('Reading annotations...')
  anntable<-read.table(annfile,header=FALSE)
  annotations<-with(anntable, GRanges(seqnames = V1, IRanges(start=V2+1, end=V3, names=V4),names=gsub("_"," ",V4)))
}

#
# report run configuration
#
message(paste(" ---------------------------------------------------"))
message(paste("        Version:", versionstr))
message(paste("          Input:", infile))
message(paste("         Output:", outfile))
message(paste("    Annotations:", annfile))
message(paste(" ---------------------------------------------------"))
message(paste("    Sample Name:", samplename))
message(paste(" ---------------------------------------------------"))


#stop('DEVSTOP')

#
# Build QC table
#
# get only control regions
results.controls <- results[which(results[["Probe type"]] == "CTRL"), ]

#
# Create summary table of CNVs (from results table, with own annotations)
#
# remove control regions
results.clinical <- results[which(results[["Probe type"]] != "CTRL"), ]
# group  by gene
# call CNA (for table summary)
# annotate table summary (if annotations supplied)

### this if from exomedepth, one might want to usinge the overlap functions of GenomicRanges
message("Annotating CNVs...")
overlap_frac <- 0.000000001  # 1bp/1Gb, basically any overlap
if (length(cnvs@CNV.calls)>0) {
  # add extra annotation
  if (all(is.na(annotations))) {
	message("No CNV annotaions available")
  } else {
    message("Adding extra annotations...")
    cnvs.annotated <- AnnotateExtra(x = cnvs.annotated,
                                    reference.annotation = annotations,
                                    min.overlap = overlap_frac,
                                    column.name = "annotation")
  }
  cnvs.annotated@annotations$name <- as.factor(
	sapply(strsplit(
	  as.character(cnvs.annotated@annotations$name), '_'), "[[", 1))
  results[[rs]] <- cnvs.annotated
} else {
  message(paste('No CNVs called with refset',toupper(rs)))
  results[[rs]] <- NA  # no CNVs called
}

#
# pepare data for plotting
#
# move non-plotting/printing from generateReport.Rnw here for clarity

#
# knit report (using grouped data & metadata)
#
# set working directory to target directory (for intermediary files)
setwd(dirname(args[2]))
# get knitr script
knitrScript <- paste(scriptDirectory, "generateReport.Rnw", sep="/")
# knit report
if (sub(".*[.]", "", outfile, perl = TRUE) == "pdf") {
  knit2pdf(knitrScript, output=sub("[.][^.]*$", ".tex", args[2], perl=TRUE))
}

#
# save workspace
#
save.image(sub("[.][^.]*$", ".RData", outfile, perl = TRUE))
