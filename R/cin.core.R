# Helper function - help calculate CIN
# @description The core function is used to calculate the CIN for the input segments. Helper function
# called by cin.standard() and cin.cytobands()
# @param segments A matrix with 3 colums: start, end, value
# @param threshold.gain A threshold above which will be set as gain
# @param threshold.loss A threshold below which will be set as loss
# @param DEF There are 2 different CIN definitions:
#      Definition 2 - sum of all abnormal probes' distances divided by the # of the abnomral probes
#      Definition 3 - sum of all abnormal segments' distances
# @param MODE There are 3 options: 'sum', 'amp' and 'del'
# @return chromosome instability index value

cin.core <- function(segments, threshold.gain = 2.5, threshold.loss = 1.5, DEF = 2, MODE = "sum") {

    cin = 0



    if (DEF == 2) {
        y = segments.to.profile(segments)
        gains = y[y > threshold.gain] - 2
        loss = 2 - y[y < threshold.loss]
        amp = sum(gains)
        del = sum(loss)
        if (MODE == "sum") {
            L = length(gains) + length(loss)
            if (L > 0)
                cin = (amp + del)/L
        } else {
            if (MODE == "amp") {
                L = length(gains)
                if (L > 0)
                  cin = amp/L
            } else {
                L = length(loss)
                if (L > 0)
                  cin = del/L
            }
        }

    }

    if (DEF == 3) {
        intensity = segments[, 3]
        gains = intensity[intensity > threshold.gain] - 2
        loss = 2 - intensity[intensity < threshold.loss]
        amp = sum(gains)
        del = sum(loss)
        if (MODE == "sum")
            cin = amp + del else {
            if (MODE == "amp")
                cin = amp else cin = del
        }
    }


    return(cin)

}
