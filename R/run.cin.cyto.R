#' Calculate cytoband CIN
#' @export
#' @import GenomicRanges
#' @description \code{run.cyto.chr} calculates cytoband level CIN for the following default thresholds
#' (with and without normalization): (a) gain threshold 2.5 and loss threshold 1.5 (b) gain threshold 2.25
#' and loss threshold 1.75 (c) gain threshold 2.10 and loss threshold 1.90. For each of these threshold
#' settings, this function will calculate CIN for gains, losses, and a combination of gains and losses
#' (referred to as 'sum' or 'overall' CIN). This will allow user to examine and select the best setting
#' of gain and loss threshold for their data.
#' More details and tutorial are given in the accompanying vignette.
#' @param grl.seg The result of any segmentation algorithm such as CBS,FMR. Should be a GRangesList
#' @param cnvgr Probe annotation info for the copy number probes - GRanges object
#' @param snpgr Probe annotation info for the SNP probes - GRanges object
#' @param genome.ucsc A Reference genome
#' @param out.folder.name Name of output folder, where the CIN objects for each setting will be created
#' @param thr.gain A numeric list that contains values set as threshold gain
#' @param thr.loss A numeric list that contains values set as threshold loss
#' @param V.mode A vector that has 3 options: 'sum', 'amp' and 'del'
#' @param V.def An integer vector that has 2 different CIN definitions - normalized (value=2) and un-normalized (value=3)
#' @param chr.num Number of chromosomes in input. Typically 22.
#' @return Creates a dataMatrix and cytobands.cin R objects for each setting that contains CIN values
#' @examples
#' #### For this example, we run cytoband CIN calculation for one setting on chromosome 1 only
#' data("grl.data") #need segment level data
#'
#' #getting genome reference file
#' data("hg18.ucsctrack")
#' hg18.ucsctrack.chr <- subset(hg18.ucsctrack, seqnames(hg18.ucsctrack) %in% "chr22")
#'
#' #get probe annotation information
#' data("cnvgr.18.auto")
#'
#' #Call function to run cytoband CIN
#' run.cin.cyto(grl.seg = grl.data, cnvgr=cnvgr.18.auto, snpgr=NULL,
#' genome.ucsc = hg18.ucsctrack.chr, thr.gain = 2.25,thr.loss = 1.75,
#' V.def = 3, V.mode="sum",chr.num = 22)
#'
#' #Run cytoband level CIN calculation for all thresholds. This is how command should be run:
#' \dontrun{
#' run.cin.cyto(grl.seg = grl.data, cnvgr=cnvgr.18.auto, snpgr=snpgr.18.auto,
#' genome.ucsc = hg18.ucsctrack)
#' }
#' # A number of RData objects will be created in 'output_cyto' folder.
#'
#' @seealso Accompanying vignette for complete end-to-end tutorial

run.cin.cyto <- function(grl.seg, cnvgr=NULL, snpgr=NULL, genome.ucsc,
                         out.folder.name = "output_cyto_cin", thr.gain = c(2.5, 2.25, 2.1),
                         thr.loss = c(1.5,1.75, 1.9), V.def = 2:3,
                         V.mode = c("sum", "amp", "del"), chr.num = 22) {
    #browser()

    if(class(grl.seg)== "GRangesList") {

        if((class(cnvgr)=="GRanges" || class(snpgr)=="GRanges") &&
           (class(genome.ucsc)=="GRanges" )) {

            if (length(thr.gain) == length(thr.loss))
            {
                # load('seg2profiles.Rdata')

                unlink(out.folder.name)
                dir.create(out.folder.name)
                for (i.thr in 1:length(thr.gain)) for (i.mode in 1:length(V.mode)) for (i.def in 1:length(V.def)) {

                    #convert GRangesList to Segment Object
                    out.seg <- convert.grl.to.seg(grl = grl.seg)

                    #convert probe annptation from GenomeicRanges into required format
                    probe.info <- process.probe.anno(cnvgr1 = cnvgr, snpgr1 = snpgr)

                    ## convert reference GRanges object into required format
                    ref.info <- process.reference.genome(genome = genome.ucsc)

                    cin.cytobands(out.seg, probe.info, ref.info, out.folder.name, thr.gain[i.thr], thr.loss[i.thr], V.def[i.def],
                                  V.mode[i.mode], chr.num)

                }
            }  else {
                stop("Length of thr.gain and thr.loss must be the same") }
        } else {
            stop("Inputs 'cnvgr', 'snpgr' and 'genome.ucsc' must be GRanges object") }
    } else {
        stop("Input 'grl.seg' must be a GRangesList object") }
}


