# Process reference genome file a specific format for CIN analysis. Helper function
# @description \code{\link{process.reference.genome}} Process reference genome file a specific format for CIN analysis. This file includes length of each chromosome of the genome,  the start and end positions, and cytoband names
# @param genome GRanges object containing reference information (incl cytobands) from UCSC
# @return Processed reference information ready for CIN analysis (list)
#' @importFrom dplyr select
# @examples
# #data(hg18.ucsctrack)
# #ref.info <- process.reference.genome(hg18.ucsctrack)

process.reference.genome <- function(genome) {
    #browser()
    ref <- {}
    temp1 <- seqlengths(genome) #save length of each chromosome
    temp2 <- as.data.frame(genome) #save GRanges object as data frame
    colnames(temp2) <- c("chrom","start","end","width","strand","name","stain")
    temp2$chrom <- as.character(temp2$chrom)
    temp2$name <- as.character(temp2$name)
    temp2$stain <- as.character(temp2$stain)

    ref$len <- temp1
    ref$cytobands <- dplyr::select(temp2, chrom, start,
                                   end, name,
                                   stain ) #select only required columns

    return(ref)
}
