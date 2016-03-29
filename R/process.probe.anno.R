# Process the Copy number andor SNP Probe annotation file into a specific format for CIN analysis. Helper function
#' @import stringr
#' @importFrom dplyr select
# @description \code{\link{process.probe.anno}} Process the Copy number / SNP Probe annotation file into a specific format for CIN analysis. The copy number probes span from a start position to end position. We calculate the mid point of this span and input into CIN anlaysis. No change is made to the SNP probe information. The the resulting ojbect is sorted based on chr, and then location.
# @param cnvgr1 GRanges object with copy number probe information: includes probe name, chromosome number
# chromosome start and end location
# @param snpgr1 GRanges object with SNP probe information: includes probe name, chromosome number and physical location
# @return Processed probe information ready for CIN analysis

process.probe.anno <- function(cnvgr1, snpgr1) {
    #browser()

    if(is.null(cnvgr1)) {
        stop("Input copy number data. This cannot be NULL")
    } else {
        cnv <- as.data.frame(cnvgr1) #convert GRange into data frame

        temp = (as.numeric(as.character(cnv$start)) + as.numeric(as.character(cnv$end)))/2 #calculate midpoint
        #cnv$midpoint <- round(temp, digits=0) #round to nearest number
        cnv$midpoint <- ceiling(temp)
        cnv.select <- dplyr::select(cnv, ID, seqnames, midpoint) #select required columns

        colnames(cnv.select) <- c("probe","chr","location") #rename columns

    }

    if(is.null(snpgr1)) {
        probe <- cnv.select
    } else {
        snp <- as.data.frame(snpgr1)
        snp.select <- dplyr::select(snp, ID, seqnames, start) #select required columns
        colnames(snp.select) <- c("probe","chr","location") #rename columns

        probe <- (rbind(cnv.select, snp.select)) # save information
    }

    ### remove the term "chr" from chromosome name
    #chrName <- sapply(strsplit(as.character(probe$chr), "chr"), function(x) x[[2]])
    chrName <- stringr::str_replace_all(probe$chr,"chr", "")

    #use 23 instead of X and 24 instead of Y - everything has to be numeric
    chrName <- stringr::str_replace_all(chrName,"X", "23")
    chrName <- stringr::str_replace_all(chrName,"Y","24")
    chrName <- stringr::str_replace_all(chrName,"MT","25")

    probe2 <- cbind(as.numeric(as.character(chrName)), as.numeric(as.character(probe$location)))
    row.names(probe2) <- probe$probe
    colnames(probe2) <- c("chr","location")

    ###sorting based on chromosome number and then location

    probe3 <- probe2[order(probe[,"chr"], probe[,"location"]),]


    return(probe3)
}
