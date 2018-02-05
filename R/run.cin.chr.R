#' Calculate chromosome CIN
#' @export
#' @import GenomicRanges
#' @description \code{\link{run.cin.chr}} calculates chromosome level CIN for the following default thresholds
#' (with and without normalization): (a) gain threshold 2.5 and loss threshold 1.5 (b) gain threshold 2.25
#' and loss threshold 1.75 (c) gain threshold 2.10 and loss threshold 1.90. For each of these threshold
#' settings, this function will calculate CIN for gains, losses, and a combination of gains and losses
#' (referred to as 'sum' or 'overall' CIN). This will allow user to examine and select the best setting
#' of gain and loss threshold for their data. More details and tutorial are given in the accompanying
#' vignette.
#'
#' @param grl.seg The result of any segmentation algorithm such as CBS,FMR. Should be a data frame of 3
#' column-lists or matrix of three-column lists
#' @param out.folder.name Name of output folder, where the CIN ojbects for each setting will be created
#' @param thr.gain  A numeric list that contains values set as threshold gain
#' @param thr.loss  A numeric list that contains values set as threshold loss
#' @param V.mode A vector that has 3 options: 'sum', 'amp' and 'del'
#' @param V.def An integer vector that has different CIN definitions (2 means normalized, 3 means
#' un-normalized)
#'
#' @return Creates a dataMatrix R object for each setting that contains CIN values
#'
#' @examples
#' # Run chromosome level CIN calculation for all thresholds. This is how command should be run:
#' # A number of RData objects will be created in 'output_chr' folder.
#' \dontrun{
#' run.cin.chr(grl.seg = grl.data)
#' }
#'
#' #For this example, we run this function for one threshold only
#'
#' data("grl.data")
#' run.cin.chr(grl.seg = grl.data, thr.gain=2.25, thr.loss=1.75, V.def=3, V.mode="sum")
#'
#' # Next step: Plot chromosome level heatmap \code{\link{comp.heatmap}}
#' # More details and tutorial are given in the accompanying vignette
#' @seealso See accompanying vignette for end-to-end tutorial

run.cin.chr <- function(grl.seg, out.folder.name = "output_chr_cin",
                        thr.gain = c(2.5, 2.25, 2.1), thr.loss = c(1.5, 1.75, 1.9),
                        V.def = 2:3, V.mode = c("sum", "amp", "del")) {
    #browser()

    #checking inputs
    if(is(grl.seg, "GRangesList")) {


        if (length(thr.gain) == length(thr.loss)) {
            unlink(out.folder.name)
            dir.create(out.folder.name)

            # convert GRrangesList into format for analysis
            seg <- convert.grl.to.seg(grl = grl.seg)

            for (i.thr in 1:length(thr.gain)) for (i.mode in 1:length(V.mode)) for (i.def in 1:length(V.def)) {
                cin.standard(seg, out.folder.name, thr.gain[i.thr], thr.loss[i.thr], V.def[i.def], V.mode[i.mode])
            }

        } else {
            stop("Length of thr.gain and thr.loss must be the same")
        }
    } else {
        stop("Input 'grl.seg' must be a GRangesList object")
    }
}
