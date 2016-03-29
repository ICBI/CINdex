# Helper function - help plot cytoband level heatmap
# @description Plots cytoband level heatmap based on group infomation. Can create pdf of png. Helper function
# called by comp.heatmap().
# @param cin A cytobands CIN matrix, each row is a cytoband, each column is a sample
# @param clinical.inf A n*2 matrix, the 1st column is 'sample name', the second is 'label'
# @param chr Chromosome number
# @param annot The annotation dataframe for cytobands including 5 columns: 'chrom', 'start', 'end',
# 'name',and 'stain'. Each row corresponds to a cytoband; which is hg18_annot.Rdata
# @param plot.choice A choice of whether the heatmaps should be .png or .pdf format
# @param base.color A choice of 'black' or 'white' base color for the heatmap (indicating no instability)
# @param title_text The saved file name and title
# @param cin.max.set Set the threshold for maximum cin value, above which will be fix to the highest
# color in heatmap
# @param color.set color set for heatmap
# @return Plots cytoband level heatmap based on group infomation

cytobands.cin.draw <- function(cin, clinical.inf, chr, annot, plot.choice, base.color, title_text = "text", cin.max.set = 5,
    color.set = "red") {
    n.cur.chr = 0
    sn.cur.chr = 0
    for (i in 1:chr) {
        if (i == 23)
            cur.chr = paste("chr", "X", sep = "") else cur.chr = paste("chr", i, sep = "")
        n.cur.chr = sum(annot[, "chrom"] == cur.chr)
        sn.cur.chr = sn.cur.chr + n.cur.chr
    }
 #browser()
    end.row.idx = sn.cur.chr
    start.row.idx = (sn.cur.chr - n.cur.chr + 1)

    if (start.row.idx > nrow(cin))
        return(NULL)

    cin = t(cin[start.row.idx:end.row.idx, ])
    if (i == 23)
        cur.chr = paste("chr", "X", sep = "") else cur.chr = paste("chr", i, sep = "")
    annot = annot[annot[, "chrom"] == cur.chr, ]

    ## re-arrange the cin and labels
    labels = matrix(clinical.inf[, 2])
    samplenames = matrix(clinical.inf[, 1])
    labels.sort = sort(labels)
    idx.sort = order(labels)
    samplenames = samplenames[idx.sort]
    cin = cin[samplenames, ]

    n.samp = nrow(cin)
    row.names = rownames(cin)

    # dev.new() width=12, height=15
    if (plot.choice == "png") {
        dev = png(filename = paste(title_text, ".png", sep = ""), pointsize = 9)
        main.title = title_text
        chr.len = annot[nrow(annot), "end"]
        plot(x = c(-1, chr.len), y = c(-1, n.samp + 1), xaxt = "n", yaxt = "n", ann = FALSE, xpd = NA, bty = "n", type = "n")
    } else {
        dev = pdf(file = paste(title_text, ".pdf", sep = ""), width = 12, height = 15, title = title_text, pointsize = 9)
        main.title = title_text
        chr.len = annot[nrow(annot), "end"]
        plot(x = c(-1, chr.len), y = c(-1, n.samp + 1), xaxt = "n", yaxt = "n", ann = FALSE, xpd = NA, bty = "n", type = "n")
    }
    # browser() palette set
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
            cin.max = cin.min + 1
    } else {
        cin.max = cin.max.set
        cin.min = 0
    }
    cin.mid = (cin.max + cin.min)/2

    whole.seq = seq(cin.min, cin.max, length = 100)

    for (i in 1:n.samp) {
        for (n.cyto in 1:ncol(cin)) {
            start = annot[n.cyto, "start"]
            end = annot[n.cyto, "end"]
            cin.value = cin[i, n.cyto]
            idx.color = tail(which(cin.value >= whole.seq), 1)
            rect(xleft = start, ybottom = i - 1, xright = end, ytop = i, col = palette[idx.color], border = NA, ljoin = 1)
        }
        text(x = -chr.len * 0.02, y = (i - 1) + 0.5, labels = row.names[i], srt = 0, cex = 1, xpd = NA)
    }

    plot.cytobands(annot, bot = n.samp + 0.05, top = n.samp + 0.2 + 0.05)

    type.classID = unique(labels.sort)
    n.subclass = rep(0, length(type.classID))
    for (i in 1:length(type.classID)) {
        n.subclass[i] = sum(labels.sort == type.classID[i])
        if (i == 1)
            mark.pos = n.subclass[i]/2 else mark.pos = (sum(n.subclass[1:i]) + sum(n.subclass[1:(i - 1)]))/2
        text(x = chr.len * 1.02, y = mark.pos, labels = type.classID[i], srt = 90, cex = 1.5, xpd = NA)
        if (i != (length(type.classID)))
            abline(h = sum(n.subclass[1:i]), col = "green", xpd = FALSE)
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
    if (chr == 23)
        chr = "X"
    text(max.pos/2, pal.bot - 2, paste("chromosome", chr, "cytobands CIN overview"), cex = 1.5, xpd = NA)
    dev.off()

}
