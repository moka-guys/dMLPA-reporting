# Base Image
FROM r-base:4.2.1

# Metadata
LABEL base.image="r-base:4.2.1"
LABEL software="dMLPA-reporting"
LABEL description="dMLPA report generation"
LABEL maintainer="David Brawand <dbrawand@nhs.net>"
LABEL website="https://github.com/moka-guys/dMLPA-reporting"
LABEL documentation="https://github.com/moka-guys/dMLPA-reporting"
LABEL license="https://github.com/moka-guys/dMLPA-reporting"

# set working directory
WORKDIR /root
ENV HOME /root

# update base system and install R package dependencies and LateX
RUN apt-get -y update
RUN apt-get --yes install \
	libcurl4-openssl-dev \
	time \
	libpoppler-cpp-dev \
	libssl-dev \
	libfontconfig1-dev \
	libxml2-dev

# install livetex (APT is not always compatible with current remote repository)
RUN cd /tmp && \
	wget -O- https://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz | tar xvfz - && \
	cd install-tl-* && \
	perl ./install-tl --no-interaction --scheme=small
ENV PATH="/usr/local/texlive/2023/bin/x86_64-linux:${PATH}"

# LateX packages
RUN tlmgr init-usertree
RUN tlmgr update --self
RUN tlmgr install framed && \
	tlmgr install mdframed && \
	tlmgr install zref && \
	tlmgr install needspace && \
	tlmgr install dingbat && \
	tlmgr install booktabs && \
	tlmgr install makecell && \
	tlmgr install colortbl && \
	tlmgr install parskip && \
	tlmgr install fancyhdr && \
	tlmgr install lastpage && \
	tlmgr install xcolor && \
	tlmgr install fp && \
	tlmgr install ms && \
	tlmgr install pgf

# set CRAN mirror
RUN echo 'options(repos=structure(c(CRAN="http://cran.ma.imperial.ac.uk/")))' > /root/.Rprofile

# R packages (BioC)
RUN Rscript -e "install.packages('BiocManager');BiocManager::install(version = '3.15', ask = FALSE)" && \
	Rscript -e "BiocManager::install(c('Biostrings','IRanges','Rsamtools','GenomicRanges','GenomicAlignments'))"

# R packages (CRAN)
RUN Rscript -e "install.packages('optparse')" && \
    Rscript -e "install.packages('readxl')" && \
    Rscript -e "install.packages('stringr')" && \
    Rscript -e "install.packages('dplyr')" && \
    Rscript -e "install.packages('knitr')" && \
    Rscript -e "install.packages('kableExtra')" && \
    Rscript -e "install.packages('tinytex')" && \
    Rscript -e "install.packages('xtable')"

# add scripts and report template
ADD generateReport.R /root
ADD generateReport.Rnw /root
ADD VERSION /root

# Rscript as entrypoint
ENTRYPOINT ["/usr/bin/Rscript", "generateReport.R"]
