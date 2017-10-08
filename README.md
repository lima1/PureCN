[![Build Status](https://travis-ci.org/lima1/PureCN.svg?branch=master)](https://travis-ci.org/lima1/PureCN)

[![PureCN_example.png](https://s9.postimg.org/6emxz4f5b/Pure_CN_example.png)](https://postimg.org/image/yer1jeiln/)

## Installation

To install this package, start R and enter:

```
source("https://bioconductor.org/biocLite.R")
biocLite("PureCN")
```

Note that if your R/Bioconductor version is outdated, this will install an old
and unsupported version.

For outdated R/Bioconductor versions, you can try backporting the latest stable
version (this should work fine for Bioconductor 3.3 and 3.4):

```
biocLite("lima1/PureCN", ref="RELEASE_3_5")
```

If you want the latest and greatest from the developer branch (since 1.7.44
recommended over stable 1.6.3):

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

## Paper

Riester M, Singh A, Brannon A, Yu K, Campbell C, Chiang D and Morrissey M
(2016). “PureCN: Copy number calling and SNV classification using targeted
short read sequencing.” _Source Code for Biology and Medicine_, **11**, pp. 13.
doi: [10.1186/s13029-016-0060-z](http://doi.org/10.1186/s13029-016-0060-z).


## License

[Artistic-2.0](https://opensource.org/licenses/Artistic-2.0).
