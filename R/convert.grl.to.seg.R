# Helper function. Convert GRangesList to list-as-a-matrix
# @description \code{\link{convert.grl.to.seg}} Function that coverts GRangesList to the
# list-as-matrix format
#' @import GenomicRanges
#' @import S4Vectors
#' @import GenomeInfoDb
#' @import IRanges
# @param grl Segment copy number information in the form of a GRangesList object
# @return Segment copy number information in list-as-a-matrix format
# @examples
# #library(GenomicRanges)
# #data(grl.data) #Segment copy number information in the form of a GRangesList object
# #out.seg <- convert.grl.to.seg(grl = grl.seg)

convert.grl.to.seg <- function(grl=NULL) {
    #browser()

    unlisted_grl <- unlist(grl, use.names=FALSE)
    all_segs <- data.frame(start=start(unlisted_grl),
                           end=end(unlisted_grl),
                           value=mcols(unlisted_grl)[ , "value"])
    f_levels <- paste(rep.int(seq_along(grl), length(seqlevels(grl))),
                      rep(seq_along(seqlevels(grl)), each=length(grl)),
                      sep="-")
    f <- factor(paste(rep.int(seq_along(grl), elementNROWS(grl)),
                      as.character(seqnames(unlisted_grl)),
                      sep="-"),
                levels=f_levels)
    seg <- lapply(split(all_segs, f),
                  function(df) {
                      m <- as.matrix(df)
                      rownames(m) <- NULL
                      m
                  })
    dim(seg) <- c(length(grl), length(seqlevels(grl)))
    dimnames(seg) <- list(names(grl), seqlevels(grl))
    seg
}

