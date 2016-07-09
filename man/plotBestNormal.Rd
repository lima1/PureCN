\name{plotBestNormal}
\alias{plotBestNormal}
\title{Plot the PCA of tumor and its best normal(s)}
\description{This function can be used to understand how a best normal is chosen
by the \code{\link{findBestNormal}} function. It can be also used 
to tune the best normal selection by finding good parameter values for
num.normals and pcs.}
\usage{plotBestNormal(gatk.normal.files, gatk.tumor.file, 
    normalDB, x = 1, y = 2, col.tumor = "red", col.best.normal = "blue", 
    col.other.normals = "black", ...)}
\arguments{
  \item{gatk.normal.files}{GATK coverage file of normal files, typically identified via 
\code{\link{findBestNormal}}.}
  \item{gatk.tumor.file}{GATK coverage file of a tumor sample.}
  \item{normalDB}{Database of normal samples, created with 
\code{\link{createNormalDatabase}}.}
  \item{x}{PC to be plotted on x-axis.}
  \item{y}{PC to be plotted on y-axis.}
  \item{col.tumor}{Color of tumor in plot.}
  \item{col.best.normal}{Color of best normals in plot.}
  \item{col.other.normals}{Color of best normals in plot.}
  \item{\dots}{Arguments passed to the plot function.}
}

\value{Returns NULL.}

\author{Markus Riester}



\seealso{\code{\link{createNormalDatabase} \link{findBestNormal}}}
\examples{
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.normal2.file <- system.file("extdata", "example_normal2.txt",
     package="PureCN")
gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file)
normalDB <- createNormalDatabase(gatk.normal.files)

gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
gatk.best.normal.file <- findBestNormal(gatk.tumor.file, normalDB)
plotBestNormal(gatk.best.normal.file, gatk.tumor.file, normalDB)

# Display sample sex. The first point in the plot is always tumor.
plotBestNormal(gatk.best.normal.file, gatk.tumor.file, normalDB,  
    pch=c(1,ifelse(normalDB$sex=="F", 1, 2)))
}
