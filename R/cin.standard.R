# Helper function - help calculate chromosome CIN
# @description Calculate chromosome instability index on chromosome level. Called by run.cin.chr()
# @param segments The result of any segmentation algorithm such as CBS,FMR, etc.Should be a data frame or matrix that consists of three columns: 'begin': the begin position of
# segments; 'end': the end position of segments; 'value': the magnitude of segments
# @param out.folder.name Name of output folder, where the CIN ojbects for each setting will be created
# @param thr.gain  A threshold above which will be set as gain
# @param thr.loss  A threshold below which will be set as loss
# @param mode There are 3 options: 'sum', 'amp' and 'del'
# @param def There are 2 different CIN definitions -normalized (value=2) and un-normalized (value=3)
# @return cin chromosome instability index value



cin.standard <- function(segments, out.folder.name, thr.gain = 2.5, thr.loss = 1.5, def = 2, mode = "sum") {
    #browser()
    n = nrow(segments)  # no. of samples
    m = ncol(segments)  # no. of chromosomes
    cin = matrix(0, n, m, dimnames = dimnames(segments))  # amplitudes of segments
    for (i in 1:n) for (j in 1:m) {

        cin[i, j] = cin.core(segments[[i, j]], thr.gain, thr.loss, def, mode)
    }
    dataMatrix = t(cin)
    if (def == 2) {
        file_name = paste("dataMatrix_", thr.gain, "_", thr.loss, "_", "normalized", "_", mode, ".RData", sep = "")
    } else if (def == 3) {
        file_name = paste("dataMatrix_", thr.gain, "_", thr.loss, "_", "unnormalized", "_", mode, ".RData", sep = "")

    }

    save(dataMatrix, file = paste(out.folder.name, "/", file_name, sep = ""))
}
