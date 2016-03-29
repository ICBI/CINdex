#' Given an input of cytobands, it outputs a list of genes that are present in the cytoband regions
#' @export
#' @param folder.name  Name of output folder
#' @param cyto.cin4heatmapObj  Output of the cytoband T test results
#' @param genome.ucsc  Reference sequence
#' @param gene.annotations  Information about CDS start and end positions, Gene names
#' @return Output files: The genes names present in the cytoband regions
#' @seealso  See accompanying vignette for an end-to-end tutorial
#' @description Once the user has a list of cytobands of interest, one downstream application could
#' be to find the list of genes present in the cytoband regions. This \code{extract.genes.in.cyto.regions}
#' function can be used for this purpose. The following steps should be run before this function can
#' be called:
#' #Step 1 : Run cytoband CIN - using \code{run.cin.chr()}
#' #Step 2: Plot cytoband level heatmap - using \code{comp.heatmap()}
#' #Step 3: Go through heatmaps as select one appropriate threshold. Load the file.
#' #Step 4: Perform T test to find differentially expressed cytobands - using \code{ttest.cyto.cin.heatmap()}
#' #Step 5: Call this funtion to extract genes located in cytoband regions
#' #More details and tutorial are given in the accompanying vignette
#' @examples
#' #For this example, we load example T test output object
#'data("cyto.cin4heatmap")
#'data("hg18.ucsctrack") #load Hg 18 reference annotation file
#'data("geneAnno") #load Gene annotations file
#'extract.genes.in.cyto.regions(cyto.cin4heatmapObj =cyto.cin4heatmap,
#'genome.ucsc = hg18.ucsctrack, gene.annotations = geneAnno)
#######################################################################################################


"extract.genes.in.cyto.regions" <- function(cyto.cin4heatmapObj = NULL,
                                            genome.ucsc = NULL,
                                            gene.annotations = NULL,
                                            folder.name = "output_genename") {
    #browser()

    if (is.null(cyto.cin4heatmapObj) || is.null(genome.ucsc) || is.null(gene.annotations)) {
        stop("One of the inputs is NULL. Please enter input in the correct format")
    }
    if(is.matrix(gene.annotations)) {

        ## convert reference GRanges object into required format
        ref.info <- process.reference.genome(genome.ucsc)
        ref <- ref.info$cytobands

        cyto.list = rownames(cyto.cin4heatmapObj)
        gene.chrs = gene.annotations[, 1]
        gene.start = as.numeric(gene.annotations[, 3])
        gene.end = as.numeric(gene.annotations[, 4])
        mark = rep(FALSE, nrow(gene.annotations))
        cytoband_genes <- {
        }
        for (i in 1:length(cyto.list)) {
            info = strsplit(cyto.list[i], "_")[[1]]  #KB changed - to _ , changed cyto.list[1] to cyto.list[i]
            chr = paste(info[1], info[2], sep = "")
            cyto = info[3]
            idx = (ref[, 1] == chr) & (ref[, 4] == cyto)
            cur.cyto.start = ref[idx, 2]
            cur.cyto.end = ref[idx, 3]
            idx.chr = (gene.annotations[, 1] == chr)
            idx.start = gene.start >= cur.cyto.start
            idx.end = gene.end <= cur.cyto.end
            idx.genes = idx.chr & idx.start & idx.end
            mark = mark | idx.genes  # KB changed || to |

            # added by KB
            geneNames <- unlist(gene.annotations[idx.genes, 5])
            tempM <- cbind(cyto.list[i], geneNames)
            cytoband_genes <- rbind(cytoband_genes, tempM)
            colnames(cytoband_genes) <- c("Cytoband", "geneNames")
        }
        # browser()
        out.genes = gene.annotations[mark, ]

        dir.create(folder.name, showWarnings = FALSE)
        filename1 <- paste(folder.name, "/", "out.genes.csv", sep = "")
        filename2 <- paste(folder.name, "/", "cytoband_genes.csv", sep = "")
        filename3 <- paste(folder.name, "/", "cytoband_genes.Rdata", sep = "")
        write.csv(out.genes, filename1, row.names = FALSE)
        write.csv(cytoband_genes, filename2, row.names = FALSE)
        save(cytoband_genes, file = filename3)

    } else {
        stop("Input 'gene.annotations' must be a matrix")
    }
}
