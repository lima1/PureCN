FROM bioconductor/bioconductor_docker:RELEASE_3_19
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

RUN apt update \
    && apt install -y --no-install-recommends apt-utils python-is-python3 \
    openjdk-17-jre-headless \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# tex support for building vignettes
# RUN apt update \
#     && apt install -y --no-install-recommends \
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
ENV GENOMICSDB_BRANCH=master
RUN mkdir $GENOMICSDB_PATH
ENV INSTALL_PREFIX=$GENOMICSDB_PATH
ENV PREREQS_ENV=$GENOMICSDB_PATH/genomicsdb_prereqs.sh
#ARG TARGETPLATFORM
#RUN if [ "$TARGETPLATFORM" = "linux/amd64" ]; then JAVA_HOME="/usr/lib/jvm/java-17-openjdk-amd64"; else JAVA_HOME=/usr/lib/jvm/java-17-openjdk-arm64; fi
#ENV MAVEN_VERSION=3.9.5

WORKDIR /tmp

RUN git clone --recursive --branch $GENOMICSDB_BRANCH https://github.com/GenomicsDB/GenomicsDB.git && \
    cd GenomicsDB/scripts/prereqs && \
    ./install_prereqs.sh && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN chmod +x $PREREQS_ENV && \
    $PREREQS_ENV && \
    cmake -DCMAKE_INSTALL_PREFIX=$GENOMICSDB_PATH -DCMAKE_BUILD_TYPE=Release ./GenomicsDB && \
    make && make install && \
    rm -rf /tmp/GenomicsDB
  
# install GenomicsDB R bindings
RUN Rscript -e 'library(remotes);\
remotes::install_github("nalinigans/GenomicsDB-R", ref="master", configure.args="--with-genomicsdb=/opt/GenomicsDB/")'

# install PureCN
RUN Rscript -e 'BiocManager::install("PureCN", dependencies = TRUE)'
#RUN Rscript -e 'BiocManager::install("lima1/PureCN", ref = "RELEASE_3_19", dependencies = TRUE)'
ENV PURECN=/usr/local/lib/R/site-library/PureCN/extdata

# add symbolic link and paths
ENV PATH $GENOMICSDB_PATH/bin:$PATH
WORKDIR /opt
RUN ln -s $PURECN /opt/PureCN

# install GATK4
ENV GATK_VERSION="4.5.0.0"
RUN wget --no-verbose https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip && \
    unzip gatk-${GATK_VERSION}.zip -d /opt && \
    rm gatk-${GATK_VERSION}.zip

ENV PATH /opt/gatk-${GATK_VERSION}:$PATH

CMD ["/bin/bash"]
