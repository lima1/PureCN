[![R-CMD-check-bioc](https://github.com/lima1/PureCN/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/lima1/PureCN/actions/workflows/check-bioc.yml)
[![BioC status](http://www.bioconductor.org/shields/build/release/bioc/PureCN.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/PureCN)
[![Platforms](http://www.bioconductor.org/shields/availability/release/PureCN.svg)](https://www.bioconductor.org/packages/release/bioc/html/PureCN.html#archives)
[![Coverage](https://img.shields.io/codecov/c/github/lima1/PureCN.svg)](https://codecov.io/gh/lima1/PureCN)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic%202.0-0298c3.svg)](https://opensource.org/licenses/Artistic-2.0) 

# PureCN

A tool developed for tumor-only diagnostic sequencing using hybrid-capture
protocols. It provides copy number adjusted for purity and ploidy and can
classify mutations by somatic status and clonality. It requires a pool of
process-matched normals for coverage normalization and artifact filtering.
PureCN was parameterized using large collections of diverse samples, ranging
from low coverage whole-exome to ultra-deep sequenced plasma gene-panels.

## Installation

To install this package, start R and enter:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("PureCN")
```

If your R/Bioconductor version is
[outdated](https://bioconductor.org/about/release-announcements/), this will
install an old and unsupported version.

For outdated R/Bioconductor versions, you can try backporting the latest stable
version (this should work fine for Bioconductor 3.3 and later):

```
BiocManager::install("lima1/PureCN", ref = "RELEASE_3_19")
```

If you want the latest and greatest from the developer branch:

```
BiocManager::install("lima1/PureCN")
```

To get the lastest stable version from
[Conda](https://anaconda.org/bioconda/bioconductor-purecn) (unstable is
currently only available from GitHub directly):

```
conda install -c bioconda bioconductor-purecn=2.8.1
```

A [Dockerhub](https://hub.docker.com/r/markusriester/purecn) image of the
latest stable version with recommended dependencies such as
[GenomicsDB](https://github.com/GenomicsDB/GenomicsDB) and
[GATK 4](https://github.com/broadinstitute/gatk) pre-installed:

```
docker pull markusriester/purecn:latest
```


## Tutorials

To get started:

```
vignette("Quick", package = "PureCN")
```

For the R package and more detailed information:

```
vignette("PureCN", package = "PureCN")
```

These tutorials are also available on the Bioconductor project page
([devel](https://bioconductor.org/packages/devel/bioc/html/PureCN.html),
[stable](https://doi.org/doi:10.18129/B9.bioc.PureCN)).

## Bugs

Before [posting](https://github.com/lima1/PureCN/issues) a bug report:

* update to the latest version 
* confirm with sessionInfo() that the latest version is used
* if this is a first PureCN attempt, closely follow the Quick vignette 
([devel](https://bioconductor.org/packages/devel/bioc/vignettes/PureCN/inst/doc/Quick.html),
[stable](https://bioconductor.org/packages/release/bioc/vignettes/PureCN/inst/doc/Quick.html))
* make sure that the issue is not covered in the Support section of the main
  vignette

## Papers

* Main paper describing the likelihood model:

    Riester M, Singh A, Brannon A, Yu K, Campbell C, Chiang D and Morrissey M
    (2016). “PureCN: Copy number calling and SNV classification using targeted
    short read sequencing.” _Source Code for Biology and Medicine_, **11**, pp. 13.
    doi: [10.1186/s13029-016-0060-z](https://doi.org/10.1186/s13029-016-0060-z).

* Validation paper, including description of novel additions, such as off-target
  support, tangent normalization and tweaks to the likelihood model:

    Oh S, Geistlinger L, Ramos M, Morgan M,  Waldron L, Riester M (2020).
    Reliable analysis of clinical tumor-only whole exome sequencing data.
    _JCO Clinical Cancer Informatics_. doi: [10.1200/CCI.19.00130](https://doi.org/10.1200/CCI.19.00130);  
    _bioRxiv_. doi: [10.1101/552711](https://doi.org/10.1101/552711)

## Selected citations

Pereira et al. (2021). "Cell-free DNA captures tumor heterogeneity and driver
alterations in rapid autopsies with pre-treated metastatic cancer". _Nature
Communications_. doi:
[10.1038/s41467-021-23394-4](https://doi.org/10.1038/s41467-021-23394-4).

Dummer et al. (2020). "Combined PD-1, BRAF and MEK inhibition in advanced
BRAF-mutant melanoma: safety run-in and biomarker cohorts of COMBI-i". _Nature
Medicine_. doi: [10.1038/s41591-020-1082-2](https://doi.org/10.1038/s41591-020-1082-2).

Bertucci et al. (2019). "Genomic characterization of metastatic breast cancers".
_Nature_. doi: [10.1038/s41586-019-1056-z](https://doi.org/10.1038/s41586-019-1056-z).

Dagogo-Jack et al. (2018). "Tracking the evolution of resistance to ALK tyrosine kinase
inhibitors through longitudinal analysis of circulating tumor DNA". _JCO
Precision Oncology_. doi:
[10.1200/PO.17.00160](https://doi.org/10.1200/PO.17.00160).

Orlando et al. (2018). "Genetic mechanisms of target antigen loss in CAR19 therapy of
acute lymphoblastic leukemia". _Nature Medicine_.
doi: [10.1038/s41591-018-0146-z](https://doi.org/10.1038/s41591-018-0146-z).

Pal et al. (2018). "Efficacy of BGJ398, a fibroblast growth factor receptor 1-3
inhibitor, in patients with previously treated advanced urothelial carcinoma
with FGFR3 alterations". _Cancer Discovery_. doi:
[10.1158/2159-8290.CD-18-0229](https://doi.org/10.1158/2159-8290.CD-18-0229).

Pitt et al. (2018). "Characterization of Nigerian breast cancer reveals
prevalent homologous recombination deficiency and aggressive molecular
features". _Nature Communications_. doi:
[10.1038/s41467-018-06616-0](https://doi.org/10.1038/s41467-018-06616-0).
