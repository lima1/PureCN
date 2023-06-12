FROM bioconductor/bioconductor_docker:RELEASE_3_17
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
    texlive \
    openjdk-17-jre-headless \
    texlive-latex-extra \
    texlive-fonts-extra \
    texlive-bibtex-extra \
    texlive-science \
    texi2html \
    texinfo \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# install GenomicsDB
ENV GENOMICSDB_PATH=/opt/GenomicsDB
RUN mkdir $GENOMICSDB_PATH
ENV INSTALL_PREFIX=$GENOMICSDB_PATH
ENV PREREQS_ENV=$GENOMICSDB_PATH/genomicsdb_prereqs.sh
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
ENV MAVEN_VERSION=3.9.2

RUN ls $JAVA_HOME

WORKDIR /tmp

RUN git clone --recursive --branch master https://github.com/GenomicsDB/GenomicsDB.git && \
    cd GenomicsDB/scripts/prereqs && \
    ./install_prereqs.sh && \
    rm -rf /var/lib/apt/lists/*

RUN chmod +x $PREREQS_ENV && \
    $PREREQS_ENV && \
    cmake -DCMAKE_INSTALL_PREFIX=$GENOMICSDB_PATH -DCMAKE_BUILD_TYPE=Release ./GenomicsDB && \
    make && make install && \
    rm -rf /tmp/GenomicsDB
  
# install GenomicsDB R bindings
RUN Rscript -e 'library(remotes);\
remotes::install_github("nalinigans/GenomicsDB-R", ref="develop", configure.args="--with-genomicsdb=/opt/GenomicsDB/")'

# install PureCN
#RUN Rscript -e 'BiocManager::install("PureCN", dependencies = TRUE)'
RUN Rscript -e 'BiocManager::install("lima1/PureCN", ref = "RELEASE_3_17", dependencies = TRUE)'
ENV PURECN=/usr/local/lib/R/site-library/PureCN/extdata

# add symbolic link and paths
ENV PATH $GENOMICSDB_PATH/bin:$PATH
WORKDIR /opt
RUN ln -s $PURECN /opt/PureCN

# install GATK4
ENV GATK_VERSION="4.4.0.0"
RUN wget --no-verbose https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip && \
    unzip gatk-${GATK_VERSION}.zip -d /opt && \
    rm gatk-${GATK_VERSION}.zip

ENV PATH /opt/gatk-${GATK_VERSION}:$PATH

CMD ["/bin/bash"]
