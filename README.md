---
title: "dMPLA report generation"
author: "David Brawand <dbrawand@nhs.net>"
date: "2023/03/01"
output:
  pdf_document:
    toc: true
    highlight: zenburn
---

# digital MPLA report generation

## Usage

The scripts takes 2 arguments and can be invoked as follows:

`generateReport.R <input file> <output file>`

### Input file

This script reads the filtered dMLPA output file. It _must_ contain two sheets:
1. `Sample info` - Sample and QC metadata
2. `Gene filtered results` - dMLPA results (including control regions)

An example file is provided in `test/dmlpa.xlsx`.

### Output file

A single PDF format report is generated as output.
Intermediary files will be written into the same folder (e.g. figures as separate PDFs).

## Usage with docker

Build the docker container which contains all scripts with `make`.

The built docker image is tagged as `seglh/dmlpa:latest`.
All reports generated with the docker image contain the version number in the report header.

The docker image contains all required packages and scripts but no reference genome.

The entrypoint is set as the report generation script and invication only requires the input and output file names.
```
docker run -it \
	-v /path_to_input:/input \
	-v /path_to_output:/output \
	seglh/dmpla:latest \
	/input/dmpla_example.xlsx \
	/output/dmpla_example.pdf \
	/input/extra.bed
```

NB: The last argument is optional (see below).


## Integrated Quality Control

The report contains QC data from the input file. Additional QC is performed on the control region data.


### Additional annotations

Additional annotations of can be supplied as a BED file (third argument).

