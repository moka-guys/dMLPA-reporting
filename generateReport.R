# R script to generate dMLPA reports from filtered excel files
#
# ENSURE PACKAGES ARE INSTALLED AND LOAD LIBRARIES TODO check if GenomicRanges, xtable, kableExtra and Rsamtools packages are required
#
required_packages <- c(
  "readxl",
  "xtable",
  "knitr",
  "kableExtra",
  "stringr",
  "dplyr",
  "tibble",
  "tidyr")

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
infile <- args[1] # output from digital MLPA script, .xlsx file format
outfile <- args[2] # name of file to output TODO consider how this will be used when called from the app
annfile <- args[3] # MLPA probe information file


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
# get probe position information from MLPA probe information sheet
# extract chr, start and end positions from probe info sheet
# conver probe positions to integers
#
annotations <- NA
# if (file.exists(annfile)) {  #TODO remove this once the alternative below is tested
#   message("Reading MLPA probe information file...")
#   probe_file <- annfile
#   probe_info <- readxl::read_excel(probe_file, sheet = "probe_info")
#   probe_coords <- probe_info %>% separate("Mapview (hg38)a", c("chrom", "coords"), ":")
#   probe_coords <- probe_coords %>% separate("coords", c("pos_start", "pos_end"), "-") 
#   probe_coords <- probe_coords %>% mutate(pos_start = as.integer(pos_start))
#   probe_coords <- probe_coords %>% mutate(pos_end = as.integer(pos_end))
#   probe_coords <- probe_coords %>% select("Probe number", chrom, pos_start, pos_end)
# }
if (file.exists(annfile)) {
  message("Reading MLPA probe information file...")
  
  probe_coords <- readxl::read_excel(annfile, sheet = "probe_info") %>%
    separate("Mapview (hg38)a", c("chrom", "pos_start", "pos_end"), ":|-") %>%
    mutate(across(c(pos_start, pos_end), as.integer)) %>%
    select("Probe number", chrom, pos_start, pos_end)
}
#
# Create summary table of CNVs (from results table, with own annotations)
#
# combine probe position information with results
# remove control regions ang get gene list
results_probes <- left_join(results,probe_coords, by = "Probe number")
results_clinical <- results_probes[which(results[["Probe type"]] != "CTRL"), ]
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

# define cutoffs for MLPA DQ values. only autosomal genes considered ("normal" == 2 copies)
cnchange <- list(
  RefToScientist1 = c(0, 0.4),
  HetDel = c(0.4, 0.6), # theoretical value 0.5
  RefToScientist2 = c(0.6, 0.8),
#  normal = c(0.8, 1.2), # theoretical value 1
  RefToScientist3 = c(1.2, 1.4),
  HetDup = c(1.4, 1.6), # theoretical value 1.5
  RefToScientist4 = c(1.6, Inf)
  )

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
      start_pos <- min(cnv_lines[["pos_start"]])
      end_pos <- max(cnv_lines[["pos_end"]])
      # calc mean and report to 3 significant figures
      copy_change_mean <- round(mean(cnv_lines[[samplename]]), 3)
      copy_change_sd <- round(sd(cnv_lines[[samplename]]), 3)
      copy_number_est <- round(mean(
        cnv_lines[["Normal copy number"]] *
        cnv_lines[[samplename]]))
      # create an aggregate line
      cnv_aggregate <- data.frame(
        Gene = unique(cnv_lines[["Gene"]]),
        Copy.Number = direction,
        Exons = ifelse(start_exon != end_exon,
          paste(start_exon, end_exon, sep = "-"), start_exon),
        Chrom = unique(cnv_lines[["chrom"]]),
        Probe.Coordinates = paste(start_pos, end_pos, sep = "-"),
        Mean.Dosage.Quotient = paste(copy_change_mean, copy_change_sd, sep = "+/-"),
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
