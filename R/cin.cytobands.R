# Helper function - help calculate cytoband CIN
# @description This function calculates chromosome instability index on cytobands level. Called by run.cin.cyto()
# @param segments The result of any segmentation algorithm such as CBS,FMR, etc . Should be a data
# frame or matrixconsists of three columns: 'begin': the begin position of
#  segments; 'end': the end position of segments; 'value': the magnitude of segments.
# @param pos A dataframe, with 2 colums chromosome and position
# @param ref A Reference sequence
# @param out.folder.name Name of output folder, where the CIN objects for each setting will be created
# @param thr.gain A threshold above which will be set as gain
# @param thr.loss Athreshold below which will be set as loss
# @param mode There are 3 options: 'sum', 'amp' and 'del'
# @param def There are 2 different CIN definitions - normalized (value=2) and un-normalized (value=3)
# @param chr.num A number indicating which chromosomes to be processed
# @return chromosome instability index value on cytobands level


cin.cytobands <- function(segments, pos, ref, out.folder.name, thr.gain = 2.5, thr.loss = 1.5, def = 2, mode = "sum", chr.num = 22) {
    #browser()

    cytobands = ref$cytobands

    n.samples = nrow(segments)  # no. of samples

    cytobands.cin = vector("list", chr.num)
    for (i in 1:length(cytobands.cin)) {
         cat("Chr:",i, '\n')
        cur.chr = paste("chr", i, sep = "")
        idx.chr = (cytobands[, "chrom"] == cur.chr)
        cur.chr.cyto = cytobands[idx.chr, ]
        ncyto.cur.chr = nrow(cur.chr.cyto)
        cur.cin.cyto = matrix(0, n.samples, ncyto.cur.chr)
        colnames(cur.cin.cyto) = cur.chr.cyto[, "name"]
        rownames(cur.cin.cyto) = rownames(segments)

        idx.pos = (pos[, 1] == i)
        probe.pos = pos[, 2][idx.pos]

        for (j in 1:ncyto.cur.chr) {
            cat(j, " ")
            start = cur.chr.cyto[j, "start"]
            end = cur.chr.cyto[j, "end"]
            for (k in 1:n.samples) {

                cur.chr.signal = segments.to.profile(segments[[k, i]])

                w.start = which(probe.pos > start)
                idx.start = w.start[1]
                w.end = which(probe.pos < end) #checking if hte probe lies within the cytoband
                idx.end = tail(w.end, 1)
                cur.cyto.signal = NULL
                if ((length(w.start) > 0) & (length(w.end) > 0))
                  cur.cyto.signal = cur.chr.signal[idx.start:idx.end]

                if (!is.null(cur.cyto.signal)) {
                  cur.cyto.seg = profile.to.segments(cur.cyto.signal)
                  cur.cin.cyto[k, j] = cin.core(cur.cyto.seg, thr.gain, thr.loss, def, mode)
                }
            } #end of k loop
        } # end of j loop

        cytobands.cin[[i]] = cur.cin.cyto
        cat("\n")
    }
    # browser()
    cytobands.cin.bindall = double()
    for (i in 1:chr.num) cytobands.cin.bindall = cbind(cytobands.cin.bindall, cytobands.cin[[i]])


    if (def == 2) {
        file_name = paste("cytobands.cin_", thr.gain, "_", thr.loss, "_", "normalized", "_", mode, ".RData", sep = "")

    } else if (def == 3) {
        file_name = paste("cytobands.cin_", thr.gain, "_", thr.loss, "_", "unnormalized", "_", mode, ".RData", sep = "")

    }

    save(cytobands.cin, file = paste(out.folder.name, "/", file_name, sep = ""))
    dataMatrix = t(cytobands.cin.bindall)

    if (def == 2) {
        file_name = paste("dataMatrix_", thr.gain, "_", thr.loss, "_", "normalized", "_", mode, ".RData", sep = "")
    } else if (def == 3) {
        file_name = paste("dataMatrix_", thr.gain, "_", thr.loss, "_", "unnormalized", "_", mode, ".RData", sep = "")

    }

    save(dataMatrix, file = paste(out.folder.name, "/", file_name, sep = ""))


}

