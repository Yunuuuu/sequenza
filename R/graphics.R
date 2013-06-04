cp.plot <- function(cp.table, map = makecmap(seq(from = median(cp.table[,3], na.rm = TRUE), 
                                                 to = max(cp.table[,3], na.rm = TRUE), by = 0.1), n = 10),
                    outlier = "white", ...) {
   require(squash)
   cell.steps   <- length(unique(cp.table[, 2]))
   color.mat    <- cp.table[, 3][seq(from = 1, to = cell.steps,by = 1)]
   cp.row.mat   <- seq(from = cell.steps + 1, to = nrow(cp.table), by = cell.steps)
   last.but.one <- length(cp.row.mat) - 1
   for (i in 1:last.but.one) {
      start     <- as.numeric(cp.row.mat[i])
      n         <- as.numeric(i + 1)
      end       <- as.numeric(cp.row.mat[n] - 1)
      color.mat <- rbind(color.mat, cp.table[, 3][seq(from = start[1], to = end[1], by = 1)])
   }
   colorgram(x=unique(x = cp.table[, 1]), y = unique(cp.table[, 2]), z = color.mat,
             colFn = jet, las = 1, xlab= "DNA-content", ylab = "cellularity", map = map, outlier = outlier, ...)
   L.max <- cp.table[which.max(cp.table[, 3]),]
   points(x = L.max[1], y = L.max[2], pch = 18)
}


# plot.fit.model <- function(mufreq.tab, cellularity, dna.content, chr23 = "XY",
#                            cn.ratio.range = c(0.5:2), avg.depth.ratio = avg.depth.ratio,
#                            cex.m = 1, cex.d = 1, ...) {
#    xy.index   <- mufreq.tab$chr == "chrX" | mufreq.tab$chr == "chrY"
#    plot(x = mufreq.tab$F, y = mufreq.tab$adjusted.ratio,
#         xlab = "mutation frequency", ylab = "depth.ratio",
#         las = 1, type="n", ...)
#    points(x = mufreq.tab$F[!xy.index], y = mufreq.tab$adjusted.ratio[!xy.index],
#           pch = 19, col = "blue", cex = cex.d)
#    types      <- types.matrix(cn.ratio.range = cn.ratio.range, chr23 = chr23)
#    if (length(which(xy.index)) >= 1) {
#       if (chr23 == "XY") {
#          points(x = mufreq.tab$F[xy.index], y = mufreq.tab$adjusted.ratio[xy.index], pch = 19, col = "green")
#       } else {
#          points(x = mufreq.tab$F[xy.index], y = mufreq.tab$adjusted.ratio[xy.index], pch = 19, col = "blue")
#       }
#    }
#    points.fit <-model.points(cellularity = cellularity, dna.content = dna.content,
#                              types = types, avg.depth.t = avg.depth.t,
#                              avg.depth.r = avg.depth.r)
#    points(points.fit,pch = 19, col = "red", cex = cex.m)
#    if (length(which(xy.index)) >= 1 & chr23 == "XY") {
#       legend("bottomright", c("Male chr X/Y",paste(paste("C", cellularity,sep = " : "), paste("P", dna.content,sep = " : "), sep = "; ")), pch = 19, col = c("green","red"))
#    } else {
#       legend("bottomright", paste(paste("C", cellularity,sep = " : "), paste("P", dna.content,sep = " : "), sep = "; "), pch = 19, col = "red")
#    }
# }

# plot.gc.depth <- function(d.values, gc.contents, ...) {
#    require(beeswarm)
#    gc.values <- sort(unique(gc.contents))
#    d.list <- list()
#    #sizevect <- rep(0,length(gc.values))
#    for ( i in 1:length(gc.values)) {
#       d.list[[i]] <- d.values[gc.contents == gc.values[i]]
#       #sizevect[i]    <- length(d.list[[i]])
#    }
#    bxplot(d.list, names = gc.values, ...)
# }

plotWindows <- function(abf.window, m.lty = 1, m.lwd = 3,
                         m.col = "black", q.bg = "lightblue", log2.plot = FALSE,
                         n.min = 1, xlim = NULL, ylim = NULL, add = FALSE, ...) {
   if (log2.plot == TRUE) {
      abf.window[, c(3, 4, 5)] <- log2(abf.window[, c(3, 4, 5)]) 
   }
   x.min        <- abf.window[1, 1]
   x.max        <- abf.window[nrow(abf.window), 2]
   y.min        <- min(abf.window[, 4])
   y.max        <- max(abf.window[, 5])
   if (add == FALSE) { 
      if (!is.null(xlim) & !is.null(ylim)) {
         plot(xlim = xlim, ylim = ylim, type = "n", x = NULL, ...)
      } else if (!is.null(xlim) & is.null(ylim)){   
         plot(xlim = xlim, ylim = c(y.min, y.max), type = "n", x = NULL, ...)
      } else if (is.null(xlim) & !is.null(ylim)) {
         plot(xlim = c(x.min, x.max), ylim = ylim, type = "n", x = NULL, ...)
      } else {
         plot(xlim = c(x.min, x.max), ylim = c(y.min, y.max), type = "n", x = NULL, ...)
      }
   }
   abf.window <- abf.window[abf.window[, 6] >= n.min, ]
   rect(xleft = abf.window[, 1], ybottom = abf.window[, 4], 
        xright = abf.window[, 2], ytop = abf.window[, 5],
        col = q.bg, border = NA)
   segments(y0 = abf.window[, 3], x0 = abf.window[, 1] , x1 = abf.window[, 2], lty = m.lty, lwd = m.lwd, col = m.col)

}

mutation.colors <- function(mut.type = NULL, l.pos = "topright") {
   col1 <- rgb(red = 0, green = 178, blue = 238, alpha = 120, max = 255)
   col2 <- rgb(red = 255, green = 64, blue = 64, alpha = 120, max = 255)
   col3 <- rgb(red = 34, green = 139, blue = 34, alpha = 120, max = 255)
   col4 <- rgb(red = 139, green = 90, blue = 0, alpha = 120, max = 255)
   col5 <- rgb(red = 127, green = 0, blue = 255, alpha = 120, max = 255)
   col6 <- rgb(red = 255, green = 215, blue = 0, alpha = 120, max = 255)
   g1   <- c('A>C', 'T>G')
   g2   <- c('A>G', 'T>C')
   g3   <- c('A>T', 'T>A')
   g4   <- c('C>A', 'G>T')
   g5   <- c('C>G', 'G>C')
   g6   <- c('C>T', 'G>A')
   colors.vect <- c(rep(x = col1, times = 2),
                    rep(x = col2, times = 2),
                    rep(x = col3, times = 2),
                    rep(x = col4, times = 2),
                    rep(x = col5, times = 2),
                    rep(x = col6, times = 2))
   all.types   <- c(g1, g2, g3, g4, g5, g6)
   if (is.null(mut.type)) {
      g1 <- paste(g1, collapse = ", ")
      g2 <- paste(g2, collapse = ", ")
      g3 <- paste(g3, collapse = ", ")
      g4 <- paste(g4, collapse = ", ")
      g5 <- paste(g5, collapse = ", ")
      g6 <- paste(g6, collapse = ", ")
      legend(x = l.pos, legend = c(g1, g2, g3, g4, g5, g6), col = c(col1, col2, col3, col4, col5, col6), pch = 15, bty = "n")
   } else {
      mut.to.col   <- setNames(colors.vect, all.types)
      as.character(mut.to.col[mut.type])
   }
}

chromosome.view <- function(baf.windows, ratio.windows, mut.tab = NULL,
                            segments = NULL, main = "", CNr = 2,
                            cellularity = NULL, dna.content = NULL, avg.depth.ratio = NULL) {
   if (is.null(segments)) {
      data.model <- NULL
      } else {
         if ("CNt" %in% colnames(segments)) {
            if (length(c(cellularity, dna.content, avg.depth.ratio)) != 3) {
               data.model <- NULL
               } else {
            data.model     <- list()
            CNt.max        <- max(segments$CNt, na.rm = TRUE)
            CNt.min        <- 0
            data.model$baf <- theoric.baf(CNr = 2, CNt = CNt.max, cellularity = cellularity)
            types          <- types.matrix(CNt.min = CNt.min, CNt.max = CNt.max, CNr = CNr)
            data.model$muf <- cbind(types, model.points(cellularity = cellularity, dna.content = dna.content,
                                                   types = types, avg.depth.ratio = avg.depth.ratio))
            }
         } else {
            data.model <- NULL
         }
      }
   if (is.null(mut.tab)) {
      par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0), mfcol = c(2,1), xaxt='n')
      min.x <- min(c(min(baf.windows$start), min(ratio.windows$start)))
      max.x <- max(c(max(baf.windows$end), max(ratio.windows$end)))
      xlim <- c(min.x, max.x)
   } else {
      min.x <- min(c(min(baf.windows$start), min(ratio.windows$start), min(mut.tab$n.base)))
      max.x <- max(c(max(baf.windows$end), max(ratio.windows$end), max(mut.tab$n.base)))
      xlim <- c(min.x, max.x)      
      par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0), mfcol = c(3,1), xaxt='n')

      plot(x = mut.tab$n.base, y = mut.tab$F, 
           ylab = "mutation frequency", las = 1, pch = 19,
           col = mutation.colors(mut.tab$mutation),
           ylim = c(min(mut.tab$F, na.rm = TRUE), 1), xlim = xlim)
      mutation.colors()
      if (!is.null(segments)){
         abline(v = segments$end.pos, lwd = 0.7, lty = 2)
         if (!is.null(data.model)) {
            for (i in 1:nrow(segments)) {
                segments(x0 = segments$start.pos[i], x1 = segments$end.pos[i], 
                         y0 = unique(data.model$muf$mufreqs[data.model$muf$CNt == segments$CNt[i]]), lwd = 0.4, lty = "24")
            }
         }
      }

   }

   plotWindows(baf.windows, ylim = c(0, 0.5),
                     ylab = "B allele frequency", xlim = xlim, las = 1)
   if (!is.null(segments)){
      abline(v = segments$end.pos, lwd = 0.7, lty = 2)
      segments(x0 = segments$start.pos, y0 = segments$Bf, x1=segments$end.pos, y1 = segments$Bf, col = "red", lwd = 3)
      if (!is.null(data.model)) {
         for (i in 1:nrow(segments)) {
            segments(x0 = segments$start.pos[i], x1 = segments$end.pos[i], 
                     y0 = unique(data.model$baf$BAF[data.model$baf$CNt == segments$CNt[i]]), lwd = 0.4, lty = "24")
         }
      }
   }
   plotWindows(ratio.windows, ylab = "depth ratio", 
                      las = 1, n.min = 1e4, ylim = c(0, 2))
   if (!is.null(segments)){
      abline(v = segments$end.pos, lwd = 0.7, lty = 2)
      segments(x0 = segments$start.pos, y0 = segments$depth.ratio, x1=segments$end.pos, y1 = segments$depth.ratio, col = "red", lwd = 3)
      if (!is.null(data.model)) {
         ratios.theoric <- unique(data.model$muf[,c('CNt', 'depth.ratio')])

         segments(x0 = rep(min(segments$start.pos, na.rm =TRUE), times = nrow(ratios.theoric)),
                  x1 = rep(max(segments$end.pos, na.rm = TRUE), times = nrow(ratios.theoric)), 
                  y0 = ratios.theoric$depth.ratio, lwd = 0.4, lty = "24")
         text(x = rep(min(segments$start.pos, na.rm =TRUE), times = nrow(ratios.theoric)),
              y = ratios.theoric$depth.ratio, labels = ratios.theoric$CNt, adj = c(0, -0.5), cex = 0.8)
         
      }
   }
   mtext(main, 3, outer = TRUE, cex = par("cex.main"), line = 2)
   mtext("Base Pairs 1e+06", 1, outer = TRUE, cex = par("cex.main"), line = 3)

   mtext(at = seq(xlim[1], xlim[2], by = 1e7), text = round(seq(xlim[1]/1e6, xlim[2]/1e6, by = 10), 0), side = 1, cex = 0.6)
}