# R script to generate dMLPA reports from filtered excel files
#
# ENSURE PACKAGES ARE INSTALLED AND LOAD LIBRARIES
#
required_packages <- c(
  "GenomicRanges",
  "readxl",
  "xtable",
  "knitr",
  "kableExtra",
  "Rsamtools",
  "stringr",
  "dplyr",
  "tibble")

# load packages or install from CRAN
options(menu.graphics = FALSE)
options(repos = structure(c(CRAN = "http://cran.ma.imperial.ac.uk/")))
for (p in required_packages) {
  if (!require(p, warn.conflicts = FALSE, character.only = TRUE)) {
    install.packages(p, dependencies = TRUE)
    library(p, warn.conflicts = FALSE, character.only = TRUE)
  }
}

#
# READ ARGS
#
cmd <- commandArgs(trailingOnly = FALSE)
script <- unlist(strsplit(cmd[which(substr(cmd, 1, 7) == "--file=")], "="))[2]
scriptdir <- normalizePath(dirname(script))

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]
annfile <- args[3]


# define version string (this is imported from the docker image if available)
versionfile <- paste(scriptdir, "VERSION", sep="/")
versionstr <- ifelse(file.exists(versionfile), readChar(versionfile, 7), "DEV")

#
# Read input data
#
if (!file.exists(infile)) stop(paste("Input file", infile, "does not exist"))
metadata <- readxl::read_excel(infile, sheet = "Sample info")
results <- readxl::read_excel(infile, sheet = "Gene filtered results")
results <- results %>% mutate(Exon = as.integer(Exon))
samplename <- colnames(metadata)[2]
panel <- stringr::str_extract(samplename, "Pan[0-9]+")

#
# custom annotations
# (could be used to add own annotations, eg known CNVs, DGV, etc)
#
annotations <- NA
if (file.exists(annfile)) {
  message("Reading annotations...")
  anntable <- read.table(annfile, header = FALSE)
  annotations <- with(anntable,
    GRanges(seqnames = V1, IRanges(start = V2 + 1, end = V3, names = V4),
      names = gsub("_", " ", V4)))
}
 
#
# Create summary table of CNVs (from results table, with own annotations)
#
# remove control regions ang get gene list
results_clinical <- results[which(results[["Probe type"]] != "CTRL"), ]
genes <- unique(results_clinical[["Gene"]])
gene_list <- paste(genes, collapse = ", ")

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
message(paste("          Panel:", panel))
message(paste("          Genes:", gene_list))
message(paste(" ---------------------------------------------------"))

# define cutoffs for CNA ([x[)
cnchange <- list(
  down = c(0, 0.8),
  # neutral = c(0.8, 1.2),
  up = c(1.2, Inf)
  )

# split 05-123,232 into two columns (chrom, pos in kb)
results_coordinates <- data.frame(
  t(sapply(results_clinical[["Mapview (hg38) in kb"]],
      function(x) unlist(strsplit(x, "-")))))
colnames(results_coordinates) <- c("Chrom", "PosKb")
results_coordinates <- transform(results_coordinates,
  Chrom = sub("^0", "", Chrom), PosKb = as.numeric(sub(",", "", PosKb)))

# merge coordiante columns with results
results_clinical_rownames <- rownames(results_clinical)
results_clinical <- cbind(results_clinical, as_tibble(results_coordinates))
rownames(results_clinical) <- results_clinical_rownames

# # rename column with sample name
# results <- rename(results, all_of(c("Copy Number Change" = samplename)))

# calculate CNA with set cutoffs and save data per gene for plotting
cna_summary <- data.frame()
cna_detail <- list()
for (gene in genes) {
  # get results per gene
  results_gene <- results_clinical[which(results_clinical$Gene == gene),]
  # find CN up and down
  for (direction in names(cnchange)) {
    print(paste(gene, direction))
    # get probes with given CN change
    cnc_rows <- which(cnchange[[direction]][1] <= results_gene[, samplename] &
      results_gene[, samplename] < cnchange[[direction]][2])
    # skip aggregation if no rows found
    if (length(cnc_rows) == 0) next
    # find contiguous rows (probe groups representing a uniform CN segment)
    cnc_data <- results_gene[cnc_rows, ]
    cnv_groups <- cumsum(c(1, diff(cnc_data[["Probe order"]]) != 1))
    # split results per CN segment
    cnv_list <- split(results_gene[cnc_rows, ], cnv_groups)
    # summarise and save contiguous CNA segments
    for (cnv_index in names(cnv_list)) {
      cnv_lines <- as_tibble(cnv_list[[cnv_index]])
      start_exon <- min(cnv_lines[["Exon"]])
      end_exon <- max(cnv_lines[["Exon"]])
      start_pos <- min(cnv_lines[["PosKb"]])
      end_pos <- max(cnv_lines[["PosKb"]])
      copy_change_mean <- round(mean(cnv_lines[[samplename]]), 3)
      copy_change_sd <- round(sd(cnv_lines[[samplename]]), 3)
      copy_number_est <- round(mean(
        cnv_lines[["Normal copy number"]] *
        cnv_lines[[samplename]]))
      # create an aggregate line
      cnv_aggregate <- data.frame(
        Gene = unique(cnv_lines[["Gene"]]),
        Direction = direction,
        Exons = ifelse(start_exon != end_exon,
          paste(start_exon, end_exon, sep = "-"), start_exon),
        Chrom = unique(cnv_lines[["Chrom"]]),
        Coordinates.kb = paste(start_pos, end_pos, sep = "-"),
        Copy.Number.Change = paste(copy_change_mean, copy_change_sd, sep = "±"),
        Estimated.Copy.Number = copy_number_est,
        Supporting.Probes = nrow(cnv_lines)
      )
      cna_summary <- rbind(cna_summary, cnv_aggregate)
      # save detail per gene (for plotting/detail view)
      cnv_context <- cbind(results_gene, marked_cnv = 
        results_gene[["Probe order"]] %in% cnv_lines[["Probe order"]])
      cna_detail[[nrow(cna_summary)]] <- cnv_context
    }
  }
}
print(cna_summary)

#
# knit report (using grouped data & metadata)
#
setwd(dirname(args[2]))
if (sub(".*[.]", "", outfile, perl = TRUE) == "pdf") {
  knitrscript <- paste(scriptdir, "generateReport.Rnw", sep="/")
  knit2pdf(knitrscript, output = sub("[.][^.]*$", ".tex", args[2], perl = TRUE))
}

#
# save workspace
#
save.image(sub("[.][^.]*$", ".RData", outfile, perl = TRUE))
