#' Read file containing segmentations
#'
#' Read segmentation files produced by DNAcopy, CNVkit or GATK4.
#'
#' @param seg.file File with segmentation
#' @param sampleid Sampleid, for segmentation files containing multiple samples
#' @param model.homozygous Unless \code{TRUE}, checks for very small log2-ratios
#' that cannot happen in samples with normal contamination
#' @param format File format. If missing, derived from the file
#' extension. Currently DNAcopy, and GATK4
#' (ModelSegments) format supported. CNVkit uses DNAcopy format.
#' @param zero Start position is 0-based. Default is \code{FALSE}.
#' @param verbose Verbose output.
#' @return A \code{data.frame}.
#' @author Markus Riester
#' @examples
#'
#' seg.file <- system.file("extdata", "example_seg.txt",
#'     package = "PureCN")
#' seg <- readSegmentationFile(seg.file, "Sample1")
#'
#' @export readSegmentationFile
readSegmentationFile <- function(seg.file, sampleid, model.homozygous = FALSE,
    format, zero = FALSE, verbose = TRUE) {
    if (is.null(seg.file)) return(NULL)
    seg <- read.delim(seg.file, comment.char = "@", stringsAsFactors = FALSE)
    if (missing(format)) format <- .getSegFormat(seg)
    if (format == "GATK4") seg <- .convertSegGATK4(seg, sampleid)
    if (verbose) flog.info("Loaded provided segmentation file %s (format %s).",
        basename(seg.file), format)
    .checkSeg(seg, sampleid, model.homozygous, verbose)
}
.getSegFormat <- function(seg) {
    if (colnames(seg)[1] == "CONTIG") return("GATK4")
    return("DNAcopy")
}
.convertSegGATK4 <- function(seg, sampleid) {
    lr_col_id <- "LOG2_COPY_RATIO_POSTERIOR_50"
    if (is.null(seg[[lr_col_id]])) {
        lr_col_id <- "MEAN_LOG2_COPY_RATIO"
    }
    data.frame(
        ID = sampleid,
        chrom = seg$CONTIG,
        loc.start = seg$START,
        loc.end = seg$END,
        num.mark = seg$NUM_POINTS_COPY_RATIO,
        seg.mean = seg[[lr_col_id]]
    )
}
.checkSegFormatting <- function(seg) {
    required.colnames <- c("ID", "chrom", "loc.start", "loc.end", "num.mark",
        "seg.mean")
    required.colnames2 <- c("ID", "chromosome", "start", "end", "num_probes",
        "mean")
    if (ncol(seg) > length(required.colnames)) {
        seg <- seg[, seq_along(required.colnames)]
    }
    if (identical(colnames(seg), required.colnames2)) {
        colnames(seg) <- required.colnames
    }

    if (!identical(as.character(suppressWarnings(as.numeric(seg$chrom))), as.character(seg$chrom))) {
        flog.warn("Expecting numeric chromosome names in seg.file, assuming file is properly sorted.")
        seg$chrom <- .strip.chr.name(seg$chrom, .getChrHash(seg$chrom))
    }

    if (any(is.na(seg$chrom) | is.na(seg$loc.start) | is.na(seg$loc.end))) {
        flog.warn("Coordinates in seg.file contain missing values.")
    }

    if (!identical(colnames(seg), required.colnames)) {
        .stopUserError(paste("Segmentation file expected with colnames",
                paste(required.colnames, collapse = ", ")))
    }
    seg
}
    
.checkSeg <- function(seg, sampleid, model.homozygous, verbose = TRUE) {
    # first check basic file formatting and attempt to standardize
    seg <- .checkSegFormatting(seg)
    
    # now do some sanity checks and pick the correct sample in
    # multi sample files

    # The smallest possible log-ratio is about 8
    # for 0.99 purity and high ploidy.
    # remove artifacts with lower log-ratio
    if (!model.homozygous && min(seg$seg.mean, na.rm = TRUE) < -8) {
        nBefore <- nrow(seg)
        seg <- seg[which(seg$seg.mean >= -8 | seg$num.mark >= 4), ]
        if (verbose) flog.warn("Removing %i short segments with log-ratio < -8.",
            nBefore - nrow(seg))
    }

    segs <- split(seg, seg$ID)
    matchedSeg <- match(make.names(sampleid), make.names(names(segs)))

    if (length(segs) == 1) {
        if (!is.null(sampleid) && is.na(matchedSeg)) {
            flog.warn("Provided sampleid (%s) does not match %s found in %s",
                      sampleid, names(segs)[1], "segmentation.")
        }
        matchedSeg <- 1
    } else if (is.null(sampleid)) {
        .stopUserError("seg.file contains multiple samples and sampleid missing.")
    } else if (is.na(matchedSeg)) {
        .stopUserError("seg.file contains multiple samples and sampleid does not match any.")
    } else {
        seg <- segs[[matchedSeg]]
    }
    offset <- weighted.mean(seg$seg.mean, seg$num.mark, na.rm = TRUE)
    if (abs(offset) > 0.001 && verbose) {
        flog.info("Re-centering provided segment means (offset %.4f).", offset)
    }
    seg$seg.mean <- seg$seg.mean - offset
    seg
}
