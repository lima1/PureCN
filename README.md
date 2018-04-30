[![Build Status](https://travis-ci.org/lima1/PureCN.svg?branch=master)](https://travis-ci.org/lima1/PureCN)
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
source("https://bioconductor.org/biocLite.R")
biocLite("PureCN")
```

If your R/Bioconductor version is outdated, this will install an old and
unsupported version.

For outdated R/Bioconductor versions, you can try backporting the latest stable
version (this should work fine for Bioconductor 3.3 and later):

```
biocLite("lima1/PureCN", ref="RELEASE_3_7")
```

If you want the latest and greatest from the developer branch:

```
biocLite("lima1/PureCN")
```


## Tutorials

To get started:

```
vignette("Quick", package="PureCN")
```

For the R package and more detailed information:

```
vignette("PureCN", package="PureCN")
```

These tutorials are also available on the Bioconductor project page ([devel](https://bioconductor.org/packages/devel/bioc/html/PureCN.html), [stable](https://doi.org/doi:10.18129/B9.bioc.PureCN)).

## Paper

Riester M, Singh A, Brannon A, Yu K, Campbell C, Chiang D and Morrissey M
(2016). “PureCN: Copy number calling and SNV classification using targeted
short read sequencing.” _Source Code for Biology and Medicine_, **11**, pp. 13.
doi: [10.1186/s13029-016-0060-z](http://doi.org/10.1186/s13029-016-0060-z).

---

[![PureCN_example.png](https://s9.postimg.org/6emxz4f5b/Pure_CN_example.png)](https://postimg.org/image/yer1jeiln/)

Copy number profile obtained by a 500-gene panel. 
