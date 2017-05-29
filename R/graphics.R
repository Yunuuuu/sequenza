plotWindows <- function(seqz.window, m.lty = 1, m.lwd = 3,
                         m.col = "black", q.bg = "lightblue",
                         log2.plot = FALSE, n.min = 1, xlim, ylim,
                         add = FALSE, ...) {
    if (log2.plot) {
        seqz.window[, c(3, 4, 5)] <- log2(seqz.window[, c(3, 4, 5)])
    }
    if (!add) {
        if (missing(xlim))
            xlim <- c(seqz.window$start[1], seqz.window$end[nrow(seqz.window)])
        if (missing(ylim))
           ylim <- c(min(seqz.window$q0, na.rm = TRUE),
               max(seqz.window$q1, na.rm = TRUE))
        plot(xlim, ylim, type = "n", ...)
    }
    seqz.window <- seqz.window[seqz.window$N >= n.min, ]
    rect(xleft = seqz.window$start, ybottom = seqz.window$q0,
        xright = seqz.window$end, ytop = seqz.window$q1,
        col = q.bg, border = NA)
    segments(y0 = seqz.window$mean, x0 = seqz.window$start,
        x1 = seqz.window$end, lty = m.lty, lwd = m.lwd, col = m.col)

}

plot_gc <- function(gc_mat, ...) {
    colorgram(x = as.numeric(rownames(gc_mat)),
        y = as.numeric(colnames(gc_mat)),
        z = gc_mat,  ...)
}

cp.plot <- function (cp.table, xlab = "Ploidy", ylab = "Cellularity",
    zlab = "Scaled rank LPP",
    colFn = colorRampPalette(c("white", "lightblue")), ...) {
    z <- matrix(rank(cp.table$lpp), nrow = nrow(cp.table$lpp)) /
        length(cp.table$lpp)
    map <- makecmap(c(0, 1), colFn = colFn, include.lowest = TRUE)
    colorgram(x = cp.table$ploidy, y = cp.table$cellularity, z = z,
        map = map, las = 1, xlab = xlab, ylab = ylab, zlab = zlab, ...)
}

chromosome.view <- function(baf.windows, ratio.windows, mut.tab = NULL,
    segments = NULL,  min.N.baf = 1, min.N.ratio = 1e4, main = "",
    vlines = FALSE, legend.inset = c(-20 * strwidth("a", units = 'figure'), 0),
    BAF.style = "lines", CNn = 2, cellularity = NULL, ploidy = NULL,
    avg.depth.ratio = NULL, model.lwd = 1, model.lty = "24", model.col = 1,
    x.chr.space = 10) {
    if (is.null(segments)) {
      data.model <- NULL
    } else {
        if ("CNt" %in% colnames(segments)) {
            if (length(c(cellularity, ploidy, avg.depth.ratio)) != 3) {
                data.model <- NULL
            } else {
                data.model <- list()
                CNt.max <- max(segments$CNt, na.rm = TRUE) + 1
                CNt.min <- 0
                data.model$baf <- expected.baf(sd = mean(segments$sd.BAF,
                    na.rm = TRUE), CNn = CNn, CNt = CNt.max,
                    cellularity = cellularity)
                if (CNn == 2) {
                    data.model$baf <- rbind(c(0 ,0,
                        max(data.model$baf$BAF), 0), data.model$baf)
                } else {
                    data.model$baf <- rbind(c(0, 0, 1, 0), data.model$baf)
                }
                types <- types.matrix(CNt.min = CNt.min,
                    CNt.max = CNt.max, CNn = CNn)
                data.model$muf <- cbind(types,
                    model.points(cellularity = cellularity,
                        ploidy = ploidy, types = types,
                        avg.depth.ratio = avg.depth.ratio))
            }
        } else {
            data.model <- NULL
        }
    }
    if (is.null(mut.tab)) {
        par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0),
            mfcol = c(2, 1), xaxt = "n")
        min.x <- min(c(min(baf.windows$start), min(ratio.windows$start)))
        max.x <- max(c(max(baf.windows$end), max(ratio.windows$end)))
        xlim <- c(min.x, max.x)
    } else {
        min.x <- min(c(min(baf.windows$start),
            min(ratio.windows$start), min(mut.tab$position)))
        max.x <- max(c(max(baf.windows$end),
            max(ratio.windows$end), max(mut.tab$position)))
        xlim <- c(min.x, max.x)
        par(mar = c(0, 4, 0, 10), oma = c(5, 0, 4, 0),
            mfcol = c(3, 1), xaxt = "n", xpd = TRUE)
        mutation.colors <- c(
            "A>C" = rgb(red = 0, green = 178, blue = 238,
                alpha = 120, maxColorValue = 255),
            "T>G" = rgb(red =   0, green = 178, blue = 238,
                alpha = 120, maxColorValue = 255),
            "A>G" = rgb(red = 255, green = 64, blue = 64,
                alpha = 120, maxColorValue = 255),
            "T>C" = rgb(red = 255, green = 64, blue = 64,
                alpha = 120, maxColorValue = 255),
            "A>T" = rgb(red =  34, green = 139, blue = 34,
                alpha = 120, maxColorValue = 255),
            "T>A" = rgb(red =  34, green = 139, blue = 34,
                alpha = 120, maxColorValue = 255),
            "C>A" = rgb(red = 139, green = 90, blue = 0,
                alpha = 120, maxColorValue = 255),
            "G>T" = rgb(red = 139, green = 90, blue = 0,
                alpha = 120, maxColorValue = 255),
            "C>G" = rgb(red = 127, green =   0, blue = 255,
                alpha = 120, maxColorValue = 255),
            "G>C" = rgb(red = 127, green =   0, blue = 255,
                alpha = 120, maxColorValue = 255),
            "C>T" = rgb(red = 255, green = 215, blue = 0,
                alpha = 120, maxColorValue = 255),
            "G>A" = rgb(red = 255, green = 215, blue = 0,
                alpha = 120, maxColorValue = 255))
        plot(x = mut.tab$position, y = mut.tab$F,
            ylab = "Mutant allele frequency", las = 1, pch = 19,
            col = c(mutation.colors, "NA" = NA)[as.character(mut.tab$mutation)],
            ylim = c(min(mut.tab$F, na.rm = TRUE), 1), xlim = xlim)
        unique.colors <- unique(mutation.colors)
        labels <- sapply(unique.colors, function(a) {
            paste(names(mutation.colors)[mutation.colors == a],
                collapse = ", ")
            })
        legend(y = "center", x  = "right", legend = labels,
            inset = legend.inset, pch = 19, col = unique.colors,
            pt.bg = unique.colors, border = NA, bty = "n")
        if (!is.null(segments)){
            if (vlines) {
                abline(v = segments$end.pos, lwd = 1, lty = 2)
            }
            if (!is.null(data.model)) {
                for (i in 1:nrow(segments)) {
                segments(x0 = segments$start.pos[i],
                    x1 = segments$end.pos[i],
                    y0 = unique(data.model$muf$mufreqs[
                        data.model$muf$CNt == segments$CNt[i]]),
                    lwd = model.lwd, lty = model.lty, col = model.col)
                }
            }
        }
    }
    if (!is.null(segments)){
        plot(ylab = "B allele frequency", type = "n",
            x = xlim, y = c(0, 0.5), las = 1)
        plotWindows(baf.windows, ylab = "B allele frequency",
                  xlim = xlim, ylim = c(0, 0.5), las = 1,
                  n.min = min.N.baf, add = TRUE)
        if (vlines) {
            abline(v = segments$end.pos, lwd = 1, lty = 2)
        }
        segments(x0 = segments$start.pos, y0 = segments$Bf,
            x1 = segments$end.pos, y1 = segments$Bf, col = "red", lwd = 3)
        if (!is.null(data.model)) {
            for (i in 1:nrow(segments)) {
               segments(x0 = segments$start.pos[i], x1 = segments$end.pos[i],
                        y0 = unique(data.model$baf$BAF[
                            data.model$baf$CNt == segments$CNt[i]]),
                        lwd = model.lwd, lty = model.lty, col = model.col)
            }
        }
    } else {
        plotWindows(baf.windows, ylab = "B allele frequency",
            xlim = xlim, ylim = c(0, 0.5), las = 1, n.min = min.N.baf)
    }
    plotWindows(ratio.windows, ylab = "Depth ratio",
        las = 1, n.min = min.N.ratio, ylim = c(0, 2.5))
    if (!is.null(segments)){
        if (vlines) {
            abline(v = segments$end.pos, lwd = 1, lty = 2)
        }
        segments(x0 = segments$start.pos, y0 = segments$depth.ratio,
            x1 = segments$end.pos, y1 = segments$depth.ratio,
            col = "red", lwd = 3)
        if (!is.null(data.model)) {
            ratios.theoric <- unique(data.model$muf[, c("CNt", "depth.ratio")])
            segments(x0 = rep(min(segments$start.pos, na.rm = TRUE),
                    times = nrow(ratios.theoric)),
                x1 = rep(max(segments$end.pos, na.rm = TRUE),
                    times = nrow(ratios.theoric)),
                y0 = ratios.theoric$depth.ratio, lwd = model.lwd,
                lty = model.lty, col = model.col)

            axis(labels = as.character(ratios.theoric$CNt), side = 4,
                line = 0, las = 1, at = ratios.theoric$depth.ratio)
            mtext(text = "Copy number", side = 4, line = 2,
                cex = par("cex.lab") * par("cex"))
        }
    }
    par(xaxt = "s")
    axis(labels = as.character(round(seq(xlim[1] / 1e6, xlim[2] / 1e6,
            by = x.chr.space), 0)),
        side = 1, line = 0, at = seq(xlim[1], xlim[2], by = 1e6 * x.chr.space),
        outer = FALSE, cex = par("cex.axis") * par("cex"))
    mtext("Position (Mb)", side = 1, line = 3, outer = FALSE,
        cex = par("cex.lab") * par("cex"))
    mtext(main, 3, outer = TRUE, cex = par("cex.main") * par("cex"), line = 2)
}
