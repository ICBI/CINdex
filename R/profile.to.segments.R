# Helper function - Converts a copy number profile to segments
# @description Converts a copy number profile to segments. Helper function Called by cin.cytobands()
# @param  y The profile to be converted to segmental representation.
# @return A data frame consists of three columns: 'begin': the begin position of
# segments; 'end': the end position of segments; 'value': the magnitude of segments

"profile.to.segments" <- function(y) {
    I <- which(diff(y) != 0)
    segments <- data.frame(begin = as.integer(c(1, I+1)), end = as.integer(c(I, length(y))), value = as.double(0))
    segments$value <- y[segments$begin]

    return(segments)
}
