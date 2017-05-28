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
