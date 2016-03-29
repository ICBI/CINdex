# Draws heatmaps at chromosome level. PDF or PNG files will be created for each setting in the working directory. Helper function called by comp.heatmap()
#' @import grDevices
# @param cin  Regular CIN matrix, each row is a chromosome, each column is a sample(patient)
# @param clinical.inf A n*2 matrix, the 1st column is 'sample name', the second is 'label'
# @param plot.choice A choice of whether the heatmaps should be .png or .pdf format
# @param base.color A choice of 'black' or 'white' base color for the heatmap (indicating no instability)
# @param heatmap.title Title for the heatmap
# @param cin.max.set Set the threshold for maximum cin value, above which will be fix to the highest color in heatmap
# @param color.set Color for the heatmap
# @return No value returned.


"heatmap.draw" <- function(cin, clinical.inf, plot.choice, base.color = "black", heatmap.title = "Chromosome instability index",
    cin.max.set = 5, color.set = "red") {
    ## re-arrange the cin and labels browser()
    labels = matrix(clinical.inf[, 2])
    samplenames = matrix(clinical.inf[, 1])
    labels.sort = sort(labels)
    idx.sort = order(labels)
    samplenames = samplenames[idx.sort]
    cin = t(cin)
    cin = cin[samplenames, ]

    n.samp = nrow(cin)
    row.names = rownames(cin)


    if (plot.choice == "png") {
        filename = paste(heatmap.title, ".png", sep = "")
        dev = png(filename, pointsize = 9)
        main.title = heatmap.title
        chr.len = ncol(cin)
        plot(x = c(-1, chr.len), y = c(-1, n.samp + 1), xaxt = "n", yaxt = "n", ann = FALSE, xpd = NA, bty = "n", type = "n")
    } else {
        dev = pdf(file = paste(heatmap.title, ".pdf", sep = ""), width = 12, height = 15, pointsize = 9)
        main.title = heatmap.title
        chr.len = ncol(cin)
        plot(x = c(-1, chr.len), y = c(-1, n.samp + 1), xaxt = "n", yaxt = "n", ann = FALSE, xpd = NA, bty = "n", type = "n",
            main = heatmap.title)
    }

    if (base.color == "black") {
        ramp <- grDevices::colorRamp(c("black", color.set))
    } else {
        ramp <- grDevices::colorRamp(c("white", color.set))
    }
    palette = rgb(ramp(seq(0, 1, length = 100)), maxColorValue = 255)


    if (cin.max.set == -1) {
        cin.max = max(cin)
        cin.min = min(cin)
        if (cin.max == cin.min)
            cin.max == cin.min + 1e-04
    } else {
        cin.max = cin.max.set
        cin.min = 0
    }
    cin.mid = (cin.max + cin.min)/2

    whole.seq = seq(cin.min, cin.max, length = 100)

    for (i in 1:n.samp) {
        for (j in 1:ncol(cin)) {
            text(x = j, y = n.samp + 1, labels = j, srt = 90, cex = 1, xpd = NA)
            cin.value = cin[i, j]
            idx.color = tail(which(cin.value >= whole.seq), 1)
            rect(xleft = j - 0.5, ybottom = i - 1, xright = j + 0.5, ytop = i, col = palette[idx.color], border = NA, ljoin = 1)
        }
        text(x = -chr.len * 0.02, y = (i - 1) + 0.5, labels = row.names[i], srt = 0, cex = 1, xpd = NA)
    }

    type.classID = unique(labels.sort)
    n.subclass = rep(0, length(type.classID))
    for (i in 1:length(type.classID)) {
        n.subclass[i] = sum(labels.sort == type.classID[i])
        if (i == 1)
            mark.pos = n.subclass[i]/2 else mark.pos = (sum(n.subclass[1:i]) + sum(n.subclass[1:(i - 1)]))/2
        text(x = chr.len + 1, y = mark.pos, labels = type.classID[i], srt = 90, cex = 1.5, xpd = NA)
        if (i != (length(type.classID)))
            abline(h = sum(n.subclass[1:i]), col = "green", lwd = 2, xpd = FALSE)
    }


    # draw color bar

    max.pos = chr.len
    nlevels = length(palette)
    x <- seq(1, max.pos, len = nlevels + 1)
    pal.top = -0.25
    pal.bot = -0.45
    rect(x[1:nlevels], pal.bot, x[-1], pal.top, col = palette, border = NA, ljoin = 1, xpd = NA)
    segments(x0 = c(1, max.pos/2, max.pos), y0 = pal.bot, x1 = c(1, max.pos/2, max.pos), y1 = pal.top, col = "black", lend = 1,
        lwd = 2, xpd = NA)  # draw ticks
    text(c(1, max.pos/2, max.pos), pal.bot - 0.05, c(sprintf("%.2f", cin.min), sprintf("%.2f", cin.mid), sprintf("%.2f",
        cin.max)), adj = c(0.5, 1), xpd = NA)

    text(max.pos/2, pal.bot - 2, heatmap.title, cex = 1.5, xpd = NA)
    dev.off()
}



