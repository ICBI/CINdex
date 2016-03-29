#' A comprehensive heatmap function that plots Chromosome and Cytoband heatmaps
#' @export
#' @import png
#' @import gridExtra
#' @description When the \code{run.cin.chr} and \code{run.cyto.chr} functions are
#' called, we get Chromosome and Cytoband CIN values for various gain/loss threshold settings.
#' This \code{comp.heatmap} function can be used to pick the best threshold for the input data.
#' It plots heatmaps for two groups of interest (case and control) for all the input gain/loss threshold
#' settings. By visually checking the heatmaps, the user can pick the threshold/setting that shows the best
#' contrast between two groups of interest.
#' Steps:
#' #Step 1: Run cytoband CIN or chromosome CIN - using \code{run.cin.chr()} or \code{run.cin.cyto()}
#' #Step 2: Call this function to create chromosome or cytoband level heatmaps. Pick gain/loss threshold
#' appropriate for data.
#' See vignette for more details.
#'
#' @param R_or_C The value'Regular' plots chromosome level heatmap and 'Cytobands' plots cytoband level
#' heatmaps
#' @param clinical.inf An n*2 matrix, the 1st column is 'sample name', the second is 'label'
#' @param genome.ucsc A Reference genome
#' @param in.folder.name Name of folder where the Chromsome CIN or Cytoband CIN objects are present
#' @param out.folder.name Name of folder where the Chromosome heatmaps or Cytoband heatmaps will be saved
#' @param plot.choice A choice of whether the heatmaps should be .png or .pdf format
#' @param base.color A choice of 'black' or 'white' base color for the heatmap (indicating no instability)
#' @param thr.gain A threshold above which will be set as gain
#' @param thr.loss A threshold below which will be set as loss
#' @param V.mode There are 3 options: 'sum', 'amp' and 'del'
#' @param V.def There are 2 different CIN definitions - normalized (value=2) and un-normalized (value=3)
#' @return No value returned. If R_or_C='Regular', it will genearte chromosome level heatmap,
#' If R_or_C='Cytobands',it will generate cytoband level heatmap
#' @examples
#' ###### Example 1 - Chromosome level
#'
#' ## Step 1: Run chromosome CIN
#' # This is how command should be run:
#' \dontrun{
#' run.cin.chr(grl.seg = grl.data)
#' }
#' # For this example, we run chr CIN on one threshold only
#' data("grl.data")
#' run.cin.chr(grl.seg = grl.data, thr.gain=2.25, thr.loss=1.75, V.def=3, V.mode="sum")
#'
#' ## Step 2: Plot chromosome level heatmap
#' # This is how the command must be called:
#' \dontrun{
#' comp.heatmap(R_or_C="Regular", clinical.inf=clin.crc, genome.ucsc=hg18.ucsctrack, thr.gain = 2.25,
#' thr.loss = 1.75,V.def = 3,V.mode = "sum")
#'}
#' # For this example, we run chr heatmap on one threshold only
#' comp.heatmap(R_or_C='Regular', clinical.inf=clin.crc, genome.ucsc=hg18.ucsctrack, thr.gain = 2.25,
#' thr.loss = 1.75,V.def = 3,V.mode = "sum")
#'
#' ###### Example 2 - Cytoband level
#'
#' ## Step 1 : Run cytoband CIN
#' # This is how command should be run:
#' \dontrun{
#' run.cin.cyto(grl.seg = grl.data,cnvgr=cnvgr.18.auto, snpgr=snpgr.18.auto,
#' genome.ucsc = hg18.ucsctrack)
#'
#' ## Step 2: Plot cytoband level heatmap
#'
#' comp.heatmap(R_or_C="Cytobands", clinical.inf=clin.crc, genome.ucsc=hg18.ucsctrack,
#' thr.gain=2.25, thr.loss=1.75,V.def=3,V.mode="sum")
#' }
#'
#' @seealso See accompanying vignette for end-to-end tutorial


comp.heatmap <- function(R_or_C = "Regular", clinical.inf = NULL, genome.ucsc = NULL,
                         in.folder.name = "output_chr_cin", out.folder.name = "output_chr_plots",
                         plot.choice = "png", base.color = "black", thr.gain = c(2.5, 2.25, 2.1),
                         thr.loss = c(1.5, 1.75, 1.9), V.def = 2:3, V.mode = c("sum","amp", "del")) {

    if(is.matrix(clinical.inf)) {

        if(class(genome.ucsc)=="GRanges") {

            color.sum = "yellow"
            color.amp = "red"
            color.del = "blue"
            color.set = c(color.sum, color.amp, color.del)

            ## convert reference GRanges object into required format
            ref <- process.reference.genome(genome.ucsc)

            if (length(thr.gain) == length(thr.loss)) {

                if (R_or_C == "Regular") {
                    unlink(out.folder.name)
                    dir.create(out.folder.name, showWarnings = FALSE)

                    for (i.thr in 1:length(thr.gain)) for (i.mode in 1:length(V.mode)) {

                        if (V.mode[i.mode] == "del")
                            cin.max.set = 2 else cin.max.set = 8


                            for (i.def in 1:length(V.def)) {
                                def = V.def[i.def]

                                if (def == 2) {
                                    file_name_p = paste("dataMatrix_", thr.gain[i.thr], "_", thr.loss[i.thr], "_", "normalized", "_", V.mode[i.mode],
                                                        sep = "")

                                } else if (def == 3) {
                                    file_name_p = paste("dataMatrix_", thr.gain[i.thr], "_", thr.loss[i.thr], "_", "unnormalized", "_", V.mode[i.mode],
                                                        sep = "")

                                }
                                file_name = paste(file_name_p, ".RData", sep = "")
                                load(paste(in.folder.name, "/", file_name, sep = ""))
                                file_name_p <- paste(out.folder.name, "/", file_name_p, sep = "")
                                heatmap.draw(dataMatrix, clinical.inf, plot.choice, base.color, file_name_p, cin.max.set, color.set[i.mode])



                            }

                    }

                    #dev.off()
                    #consolidating all images into one pdf file
                    if (plot.choice == "png") {
                        #reading all the files and putting it into one pdf file
                        #browser()

                        path1 <- paste(getwd(),"/",out.folder.name,"/",sep="")
                        fileNames <-  dir(path =path1)
                        fullPath <- paste(path1,fileNames, sep="")

                        pdf("all.chr.plots.pdf")
                        rl = lapply(sprintf(fullPath), png::readPNG)
                        gl = lapply(rl, grid::rasterGrob)
                        do.call(gridExtra::grid.arrange, c(gl, ncol=3))
                        #dev.off()
                        graphics.off()

                    }

                }# end of Regular


                if (R_or_C == "Cytobands") {

                    unlink(out.folder.name)
                    dir.create(out.folder.name, showWarnings = FALSE)

                    for (chr in 1:22) {
                        cat(chr, " ")
                        for (i.thr in 1:length(thr.gain)) for (i.mode in 1:length(V.mode)) {
                            if (V.mode[i.mode] == "del")
                                cin.max.set = 2 else cin.max.set = 8
                                for (i.def in 1:length(V.def)) {
                                    def = V.def[i.def]

                                    if (def == 2) {
                                        file_name_p = paste("dataMatrix_", thr.gain[i.thr], "_", thr.loss[i.thr], "_", "normalized", "_", V.mode[i.mode],
                                                            sep = "")
                                        folder_name = paste(out.folder.name, "/", thr.gain[i.thr], "_", thr.loss[i.thr], "_", "normalized",
                                                            "_", V.mode[i.mode], sep = "")

                                    } else if (def == 3) {
                                        file_name_p = paste("dataMatrix_", thr.gain[i.thr], "_", thr.loss[i.thr], "_", "unnormalized", "_",
                                                            V.mode[i.mode], sep = "")
                                        folder_name = paste(out.folder.name, "/", thr.gain[i.thr], "_", thr.loss[i.thr], "_", "unnormalized",
                                                            "_", V.mode[i.mode], sep = "")

                                    }
                                    unlink(folder_name)
                                    dir.create(folder_name, showWarnings = FALSE)
                                    file_name = paste(file_name_p, ".RData", sep = "")
                                    file_path = paste(folder_name, "/chr_", chr, "_", file_name_p, sep = "")
                                    #browser()
                                    if(file.exists(paste(in.folder.name, "/", file_name, sep = ""))) {
                                        load(file = paste(in.folder.name, "/", file_name, sep = ""))
                                        cytobands.cin.draw(dataMatrix, clinical.inf, chr,
                                                           ref$cytobands, plot.choice,
                                                           base.color, file_path,
                                                           cin.max.set,
                                                           color.set[i.mode])
                                    }
                                }
                        }
                    } #end of for - chr

                }# end of Cytobands



            } else {
                stop("Length of thr.gain and thr.loss must be the same")
            }

        } else {
            stop("Input 'genome.ucsc' has to be a GRanges object")
        }

    } else {
        stop("Input Clinical must be a n*2 matrix, the 1st column is 'Sample', the second is 'Label'.Please re-format the input")
    }

} #end of function


