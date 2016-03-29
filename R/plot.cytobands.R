# Helper function - plot cytoband on heatmap
# @description Plots cytobands on a heatmap. Helper function called by cytobands.cin.draw()
#' @import bitops
#' @importFrom graphics abline par plot polygon rect segments text
# @param cb A data frame consists of columns 'chrom', 'start', 'end', 'name' and 'stain'.
# Each row corresponds to a cytoband;
# @param bot The bottom position of a horizontally drawn chromosome;
# @param top The top position of the chromosome.
# @return Plots cytobands on a heatmap
# Note that the unit of x-axis must be nucliotide base.

"plot.cytobands" <- function(cb, bot, top) {

    # require('bitops') devtools::use_package('bitops')

    get.stain.col <- function(stain) {
        switch(stain, gneg = "gray100", gpos25 = "gray90", gpos50 = "gray70", gpos75 = "gray40", gpos100 = "gray0", gvar = "gray100",
            stalk = "brown3", acen = "brown4", "white")
    }

    cb <- cbind(cb, shape = rep(0L, nrow(cb)))

    # determine shapes of bands
    I <- grep("^p", cb$name)
    cb[I[1], "shape"] <- 1L
    cb[tail(I, 1), "shape"] <- 2L

    I <- grep("^q", cb$name)
    cb[I[1], "shape"] <- 1L
    cb[tail(I, 1), "shape"] <- 2L

    for (i in which(cb$stain == "stalk")) {
        cb[i, "shape"] <- as.integer(NA)
        if (i > 1)
            cb[i - 1, "shape"] <- bitops::bitOr(cb[i - 1, "shape"], 2L)
        if (i < nrow(cb))
            cb[i + 1, "shape"] <- bitops::bitOr(cb[i + 1, "shape"], 1L)
    }

    # convert stain to colors
    cb$stain <- sapply(cb$stain, get.stain.col)

    # draw cytobands
    for (i in 1:nrow(cb)) {
        start <- cb$start[i]
        end <- cb$end[i]
        color <- cb$stain[i]
        if (is.na(cb$shape[i])) {
            d <- (end - start)/3
            segments(start + d, top, start + d, bot, col = color, xpd = NA)
            segments(end - d, top, end - d, bot, col = color, xpd = NA)
        } else {
            d <- (end - start)/4  # horizontal and vertical sizes of corners
            d.y <- (top - bot)/4
            if (cb$shape[i] == 3L)
                polygon(c(start, start + d, end - d, end, end, end - d, start + d, start), c(top - d.y, top, top, top - d.y,
                  bot + d.y, bot, bot, bot + d.y), col = color, xpd = NA) else if (cb$shape[i] == 1L)
                polygon(c(start, start + d, end, end, start + d, start), c(top - d.y, top, top, bot, bot, bot + d.y), col = color,
                  xpd = NA) else if (cb$shape[i] == 2L)
                polygon(c(start, end - d, end, end, end - d, start), c(top, top, top - d.y, bot + d.y, bot, bot), col = color,
                  xpd = NA) else if (cb$shape[i] == 0L)
                polygon(c(start, end, end, start), c(top, top, bot, bot), col = color, xpd = NA) else stop("Invalid cytoband shape.")
        }
    }

    # draw cytoband names
    text((cb$start + cb$end)/2, top, cb$name, srt = 90, cex = 1, adj = c(-0.15, 0.5), xpd = NA)
}
