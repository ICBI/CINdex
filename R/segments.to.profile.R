# Helper function - Construct a profile from segments
# @description Construct a profile from segments. Helper function called by cin.cytobands()
# @param seg If seg is an 3-element array the recovered profile only
# contains one segment. If seg is an n by 3 matrix, then
# seg[i,1], seg[i,2], and seg[i,3] are the begin indices,
# end indices, and intensity of the i-th segment.
#
# @param sig.prf If `sig.prf' is not NULL, the intensity of a segment is
# recalculated as the mean of the signals in sig.prf within
# the segment, i.e. mean(sig.prf[seg[i,1]:seg[i,2]])
# overrides the constant in seg[i,3] if `seg' has 3 columns.
# If `sig.prf' is not null `seg' can have only two columns
# indicating the begin and end positions of segments.
# @return The reconstructed copy number profile (a vector of length max(seg[,2]))
#
# Remarks:
#   It is safe to pass a segment array with m > 3 columns to the function as
#   long as the first three columns have the require meanings.
#
#
segments.to.profile <- function(seg, sig.prf = NULL) {
    #browser()
    if (!(is.numeric(seg) || is.data.frame(seg))) {
        stop("Wrong representation of segments.")
    }

    if (is.vector(seg))
        seg = t(seg)

    y = rep(as.double(NA), max(seg[, 2]))

    if (is.null(sig.prf)) {
        for (i in 1:nrow(seg)) y[seg[i, 1]:seg[i, 2]] = seg[i, 3]
    } else {
        for (i in 1:nrow(seg)) y[seg[i, 1]:seg[i, 2]] = mean(sig.prf[seg[i, 1]:seg[i, 2]])
    }
    return(y)
}
