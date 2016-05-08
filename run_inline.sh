#!/bin/sh

#/usr/local/bin/R --vanilla <<RSCRIPT
R --no-save <<RSCRIPT
library(inlinedocs);
package.skeleton.dx('.')
RSCRIPT

cd ..

# do the full thing only for final release
if [ 1 == 0 ]; then
    R CMD build --no-build-vignettes PureCN
else
    /usr/local/bin/R CMD build PureCN
    rm PureCN.pdf
    /usr/local/bin/R CMD Rd2pdf PureCN
fi    
cd -


