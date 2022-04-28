FROM bioconductor/bioconductor_docker:RELEASE_3_15
#FROM bioconductor/bioconductor_docker:devel

# install base packages
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}; \
    BiocManager::install(c("TxDb.Hsapiens.UCSC.hg38.knownGene", "TxDb.Hsapiens.UCSC.hg19.knownGene"))'
RUN Rscript -e 'install.packages(c("optparse", "R.utils")); \
    BiocManager::install(c("remotes", "raerose01/deconstructSigs"));'
RUN Rscript -e 'BiocManager::install(c("GenomicRanges", "IRanges", "DNAcopy", "Biostrings", "GenomicFeatures", "rtracklayer",\
"S4Vectors", "rhdf5", "VariantAnnotation", "Rsamtools", "BiocGenerics"))'

# patched PSCBS with support of interval weights
RUN Rscript -e 'BiocManager::install("lima1/PSCBS", ref="add_dnacopy_weighting")'

RUN apt-get update \
    && apt-get install -y --no-install-recommends apt-utils python-is-python3 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# RUN apt-get install -y --no-install-recommends \
#     texlive \
#     texlive-latex-extra \
#     texlive-fonts-extra \
#     texlive-bibtex-extra \
#     texlive-science \
#     texi2html \
#     texinfo \
#     && apt-get clean \
#     && rm -rf /var/lib/apt/lists/*

# install GenomicsDB
ENV GENOMICSDB_PATH=/opt/GenomicsDB
RUN mkdir $GENOMICSDB_PATH
ENV INSTALL_PREFIX=$GENOMICSDB_PATH
ENV PREREQS_ENV=$GENOMICSDB_PATH/genomicsdb_prereqs.sh

WORKDIR /tmp
RUN git clone -b develop https://github.com/GenomicsDB/GenomicsDB && \
    cd GenomicsDB && \
    git checkout f945f3db507c32e460033c284addc801a05b6919
RUN cd GenomicsDB/scripts/prereqs && \
    ./install_prereqs.sh && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
RUN chmod +x $GENOMICSDB_PATH/genomicsdb_prereqs.sh && \
    $GENOMICSDB_PATH/genomicsdb_prereqs.sh && \
    cmake -DCMAKE_INSTALL_PREFIX=$GENOMICSDB_PATH ./GenomicsDB && \
    make && make install && \
    rm -rf /tmp/GenomicsDB

# install GenomicsDB R bindings
RUN Rscript -e 'library(remotes);\
remotes::install_github("nalinigans/GenomicsDB-R", ref="master", configure.args="--with-genomicsdb=/opt/GenomicsDB/")'

# install PureCN
RUN Rscript -e 'BiocManager::install("PureCN", dependencies = TRUE)'
#RUN Rscript -e 'BiocManager::install("lima1/PureCN", dependencies = TRUE)'
ENV PURECN=/usr/local/lib/R/site-library/PureCN/extdata

# add symbolic link and paths
ENV PATH $GENOMICSDB_PATH/bin:$PATH
WORKDIR /opt
RUN ln -s $PURECN /opt/PureCN

# install GATK4
ENV GATK_VERSION="4.2.6.1"
RUN wget --no-verbose https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip && \
    unzip gatk-${GATK_VERSION}.zip -d /opt && \
    rm gatk-${GATK_VERSION}.zip

ENV PATH /opt/gatk-${GATK_VERSION}:$PATH

CMD ["/bin/bash"]
