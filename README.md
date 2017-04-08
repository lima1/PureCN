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
biocLite("Bioconductor-mirror/PureCN", ref="release-3.5")
```

If you want the latest and greatest from the developer branch:

```
biocLite("Bioconductor-mirror/PureCN")
```
