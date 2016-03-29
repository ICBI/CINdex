#' Performs T test on cytoband level CIN data, and plots heatmap
#' @export
#' @description \code{ttest.cyto.cin.heatmap} to perform T test to find differentially expressed cytobands.
#' It also plots a heatmap after performing heirarchical clustering. When to use this function:
#' #Step 1: Run cytoband CIN - using \code{run.cin.chr()}.
#' #Step 2: Plot cytoband level heatmap - using \code{comp.heatmap()}.
#' #Step 3: Go through heatmaps as select one appropriate threshold. Load the file.
#' #Step 4: Call this function.
#' More details and tutorial are given in the accompanying vignette
#' @importFrom gplots heatmap.2
#' @importFrom utils write.csv write.table
#' @importFrom stats t.test
#' @import som
#' @param cytobands.cin.obj (eg. cytobands.cin_2.25_1.75_unnormalized_amp.Rdata), a list in which each cell is chromosome cin matrix
#' @param clinical.inf In a clinical.inf.Rdata is a two columns array, the 1st column is samplename, the 2nd is the label
#' @param genome.ucsc  Reference sequence
#' @param file.ext Provide a meaningful file name extension. Ideally include the gain, loss threshold settings
#' @param folder.name Name of folder where the output files will be generated
#' @param combine.cyto.flag  Whether or not to save the combine cytobands as a uni array rather than a list
#' @return #Outputs:
#'      1. cyto.cin.uni.file.ext.Rdata (eg. cyto.cin.uni.gainT_lossT_unnormalized.Rdata)
#'      2. Heatmaps: eg. CIN relapse-free VS relapse for gainT_lossT_unnormalized_dendrogram.pdf
#'      3. Raw CIN array for the corresponding heatmap:
#'      #ttest.cyto.cin4heatmap.gainT_lossT_unnormalized.csv
#'      #ttest.cyto.cin4heatmap.gainT_lossT_unnormalized.Rdata
#'      4. T test results for all cytobands on the whole genome
#'      #ttest.cytobands.cin.gainT_lossT_unnormalized.txt
#' @examples
#' #For this example, we load an example cytoband CIN data
#' data("cytobands.cin")
#' data("clin.crc") # sample names with group information
#' data("hg18.ucsctrack") #hg18 reference file
#' ttest.cyto.cin.heatmap(cytobands.cin.obj = cytobands.cin,
#' clinical.inf = clin.crc, genome.ucsc = hg18.ucsctrack)
#' @seealso See accompaying vignette for a detailed end to end workflow tutorial

###################################################################################################

"ttest.cyto.cin.heatmap" <- function(cytobands.cin.obj, clinical.inf, genome.ucsc,
                                     file.ext = "gainT_lossT_unnorm", folder.name = "output_ttest",
                                     combine.cyto.flag = FALSE) {
    #browser()
    #library('gplots')

    if(is.matrix(clinical.inf)) {

        if(class(genome.ucsc)=="GRanges") {



            labels = unique(clinical.inf[, 2])
            samples.label.1 = clinical.inf[clinical.inf[, 2] == labels[1], 1]
            samples.label.2 = clinical.inf[clinical.inf[, 2] == labels[2], 1]

            n.chr = length(cytobands.cin.obj)

            cyto.cin.uni = NULL

            ## convert reference GRanges object into required format
            ref.info <- process.reference.genome(genome.ucsc)
            ref <- ref.info$cytobands

            dir.create(folder.name, showWarnings = FALSE)

            combine.cyto.Rdata.name = paste(folder.name, "/", "cyto.cin.uni.", file.ext, ".RData", sep = "")
            ttest.file.name = paste(folder.name, "/", "ttest.result.cyto.cin.", file.ext, ".csv", sep = "")
            unlink(ttest.file.name)

            out_line = c("chr", "cur.cyto.name", "start.position", "end.position", "p.value", "fold.change", "flag.variation")
            # write(out_line,file=ttest.file.name, ncolumns=7, append=TRUE,sep='\t')
            finalM <- {
            }
            for (chr in 1:n.chr) {
                cur.chr.cyto = cytobands.cin.obj[[chr]]
                if (chr == 23)
                    hg.chr = paste("chr", "X", sep = "") else hg.chr = paste("chr", chr, sep = "")
                    cur.chr.ref = ref[ref[, 1] == hg.chr, ]
                    cur.chr.cyto.tmp = cur.chr.cyto
                    colnames(cur.chr.cyto.tmp) = paste("chr_", chr, "_", colnames(cur.chr.cyto.tmp), sep = "")
                    cyto.cin.uni = cbind(cyto.cin.uni, cur.chr.cyto.tmp)
                    for (i.cyto in 1:ncol(cur.chr.cyto)) {

                        cin.label.1 = cur.chr.cyto[samples.label.1, i.cyto]
                        cin.label.2 = cur.chr.cyto[samples.label.2, i.cyto]

                        cur.cyto.name = colnames(cur.chr.cyto)[i.cyto]
                        idx.hg.chr = cur.chr.ref[, "name"] == cur.cyto.name
                        start.position = cur.chr.ref[idx.hg.chr, "start"]
                        end.position = cur.chr.ref[idx.hg.chr, "end"]

                        test.res = t.test(cin.label.1, cin.label.2)
                        pvalue = NA
                        folder.change = NA
                        flag.variation = NA
                        if (!is.na(test.res$p.value)) {
                            pvalue = sprintf("%.4f", test.res$p.value)
                            folder.change = sprintf("%.4f", test.res$estimate[2]/test.res$estimate[1])

                            if (folder.change > 1)
                                flag.variation = "+" else flag.variation = "-"
                        }

                        out_line = c(chr, cur.cyto.name, start.position, end.position, pvalue, folder.change, flag.variation)
                        colnames(out_line)

                        # write(out_line,file=ttest.file.name, ncolumns=7, append=TRUE,sep='\t')
                        finalM <- rbind(finalM, out_line)
                        colnames(finalM) <- c("chr", "cur.cyto.name", "start.position", "end.position", "p.value", "fold.change", "flag.variation")
                    }


            }
            #browser()
            write.table(finalM, file = ttest.file.name, sep = ",", row.names = FALSE)

            # browser()
            if (combine.cyto.flag) {
                # cyto.cin.uni=t(cyto.cin.uni)
                save(cyto.cin.uni, file = combine.cyto.Rdata.name)
                # write.csv(cyto.cin.uni[idx,],file= paste(folder.name,'/','cyto.cin.uni.sig.csv', sep=''))

            }
            # browser() ttest.res = read.table(file = ttest.file.name, header = TRUE, sep = '\t')

            ttest.res <- finalM

            idx.non.na = !is.na(ttest.res[, "p.value"])
            idx.significant = as.double(ttest.res[, "p.value"]) < 0.05
            # browser()
            idx = idx.non.na & idx.significant

            cyto.cin.uni = t(cyto.cin.uni)  ###added by kb, else code would give error

            cyto.cin4heatmap = cbind(cyto.cin.uni[idx, samples.label.1], cyto.cin.uni[idx, samples.label.2])
            save(cyto.cin4heatmap, file = paste(folder.name, "/", "cyto.cin4heatmap.", file.ext, ".RData", sep = ""))
            output.label = clinical.inf[match(colnames(cyto.cin4heatmap), clinical.inf[, 1]), 2]

            cyto.cin4heatmap.withLabel <- t(cbind(output.label, t(cyto.cin4heatmap)))
            # save(cyto.cin4heatmap.withLabel, file=paste(folder.name,'/','cyto.cin4heatmapWithLabel.', file.ext,'.RData',sep='') )

            write.csv(t(cbind(output.label, t(cyto.cin4heatmap))), file = paste(folder.name, "/", "cyto.cin4heatmap.", file.ext,
                                                                                ".csv", sep = ""))

            # require(gplots) require(som)
            cyto.cin4heatmap2 = som::normalize(cyto.cin4heatmap, byrow = TRUE)
            dimnames(cyto.cin4heatmap2) = dimnames(cyto.cin4heatmap)
            cyto.cin4heatmap = cyto.cin4heatmap2
            col.cin = c(rep("blue", length(samples.label.1)), rep("red", length(samples.label.2)))

            dev = pdf(file = paste(folder.name, "/", "CIN_", labels[1], "_VS_", labels[2], "_for_", file.ext, "_dendrogram.pdf",
                                   sep = ""))
            main.name = paste("CIN ", labels[1], " VS ", labels[2], " for ", strsplit(file.ext, "_")[[1]][1], "\n", "blue: ", labels[1],
                              "\n", "red: ", labels[2], sep = "")
            par(mar = c(7, 4, 4, 2) + 0.1)
            gplots::heatmap.2(cyto.cin4heatmap, ColSideColors = col.cin, col = "greenred", dendrogram = "both", Colv = TRUE, Rowv = TRUE,
                              scale = "row", trace = "none", main = main.name, labRow = rownames(cyto.cin4heatmap), cexRow = 0.75, cexCol = 0.8,
                              margins = c(12, 8))

            #dev.off()
            #graphics.off()

        } else {
            stop("Input 'genome.ucsc' must be a GRanges object")
        }

    } else {
        stop("Input Clinical must be a n*2 matrix, the 1st column is 'sample name', the second is 'label'.Please re-format the input")
    }
}
