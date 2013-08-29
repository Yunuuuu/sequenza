cp.plot <- function(cp.table, map = makecmap(seq(from = median(cp.table$L, na.rm = TRUE), 
                                                 to = max(cp.table$L, na.rm = TRUE), by = 0.1), n = 10),
                    outlier = "white", ...) {
   require(squash)
   z <- tapply(cp.table[, 'L'], list(cp.table[, 'dna.index'], cp.table[, 'cellularity']), mean)
   x <- as.numeric(rownames(z))
   y <- as.numeric(colnames(z))
   colorgram(x, y, z,
             colFn = jet, map = map, outlier = outlier, las = 1, 
             xlab= "DNA-index", ylab = "Cellularity", 
             zlab = "log-likelihood", ...)
   L.max <- cp.table[which.max(cp.table$L),]
   points(x = L.max$dna.index, y = L.max$cellularity, pch = 18)
}

cp.plot.contours <- function(cp.table, likThresh = c(0.5, 0.9, 0.99, 0.999), 
                             col = palette(), legend.pos = 'bottomright', ...) {
   require(squash)
   z <- tapply(cp.table[, 'L'], list(cp.table[, 'dna.index'], cp.table[, 'cellularity']), mean)
   x <- as.numeric(rownames(z))
   y <- as.numeric(colnames(z))
   max.lik <- max(cp.table[, 3])
   LogSumLik <- log2(sum(2^(cp.table[, 3] - max.lik))) + max.lik
   znorm <- z - LogSumLik
   znormsort <- sort(znorm, decreasing = TRUE)
   znormcumLik <- cumsum(2 ^ znormsort)
   n <- sapply(likThresh, function(x) sum(znormcumLik < x) + 1)
   logLikThresh <- znormsort[n]
   names(logLikThresh) <- paste0(likThresh * 100, '%')
   
   contour(x, y, znorm, levels = znormsort[n], col = col,
           drawlabels = FALSE,
           xlab= "DNA index", ylab = "Cellularity", ...)
   if(!is.na(legend.pos)) {
      legend(legend.pos, legend = names(logLikThresh), 
             col = col, lty = 1, title = 'cumLik')
   } 
   invisible(logLikThresh)
}

# plot.fit.model <- function(mufreq.tab, cellularity, dna.index, chr23 = "XY",
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
#    points.fit <-model.points(cellularity = cellularity, dna.index = dna.index,
#                              types = types, avg.depth.t = avg.depth.t,
#                              avg.depth.r = avg.depth.r)
#    points(points.fit,pch = 19, col = "red", cex = cex.m)
#    if (length(which(xy.index)) >= 1 & chr23 == "XY") {
#       legend("bottomright", c("Male chr X/Y",paste(paste("C", cellularity,sep = " : "), paste("P", dna.index,sep = " : "), sep = "; ")), pch = 19, col = c("green","red"))
#    } else {
#       legend("bottomright", paste(paste("C", cellularity,sep = " : "), paste("P", dna.index,sep = " : "), sep = "; "), pch = 19, col = "red")
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
                         n.min = 1, xlim, ylim, add = FALSE, ...) {
   if (log2.plot == TRUE) {
      abf.window[, c(3, 4, 5)] <- log2(abf.window[, c(3, 4, 5)]) 
   }
   if(!add) {
      if(missing(xlim)) 
         xlim <- c(abf.window[1, 1], abf.window[nrow(abf.window), 2])
      if(missing(ylim))
         ylim <- c(min(abf.window[, 4], na.rm = TRUE), max(abf.window[, 5], na.rm = TRUE))
      plot(xlim, ylim, type = "n", ...)  
   }
   abf.window <- abf.window[abf.window[, 6] >= n.min, ]
   rect(xleft = abf.window[, 1], ybottom = abf.window[, 4], 
        xright = abf.window[, 2], ytop = abf.window[, 5],
        col = q.bg, border = NA)
   segments(y0 = abf.window[, 3], x0 = abf.window[, 1] , x1 = abf.window[, 2], lty = m.lty, lwd = m.lwd, col = m.col)

}

chromosome.view <- function(baf.windows, ratio.windows, mut.tab = NULL, segments = NULL,
                            min.N.baf = 1, min.N.ratio = 1e4, main = "", vlines = FALSE, CNr = 2,
                            cellularity = NULL, dna.index = NULL, avg.depth.ratio = NULL) {
   if (is.null(segments)) {
      data.model <- NULL
      } else {
         if ("CNt" %in% colnames(segments)) {
            if (length(c(cellularity, dna.index, avg.depth.ratio)) != 3) {
               data.model <- NULL
            } else {
               data.model     <- list()
               CNt.max        <- max(segments$CNt, na.rm = TRUE) + 1
               CNt.min        <- 0
               data.model$baf <- theoretical.baf(CNr = CNr, CNt = CNt.max, cellularity = cellularity)
               if (CNr == 2) {
                  data.model$baf <- rbind(c(0,0,0.5,0), data.model$baf)
               } else {
                  data.model$baf <- rbind(c(0,0,1,0), data.model$baf)                  
               }   
               types          <- types.matrix(CNt.min = CNt.min, CNt.max = CNt.max, CNr = CNr)
               data.model$muf <- cbind(types, model.points(cellularity = cellularity, dna.index = dna.index,
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
      par(mar = c(0, 4, 0, 10), oma = c(5, 0, 4, 0), mfcol = c(3,1), xaxt='n', xpd = TRUE)
      mutation.colors <- c(
         'A>C' = rgb(red =   0, green = 178, blue = 238, alpha = 120, maxColorValue = 255),
         'T>G' = rgb(red =   0, green = 178, blue = 238, alpha = 120, maxColorValue = 255),
         'A>G' = rgb(red = 255, green =  64, blue =  64, alpha = 120, maxColorValue = 255),
         'T>C' = rgb(red = 255, green =  64, blue =  64, alpha = 120, maxColorValue = 255),
         'A>T' = rgb(red =  34, green = 139, blue =  34, alpha = 120, maxColorValue = 255),
         'T>A' = rgb(red =  34, green = 139, blue =  34, alpha = 120, maxColorValue = 255),
         'C>A' = rgb(red = 139, green =  90, blue =   0, alpha = 120, maxColorValue = 255),
         'G>T' = rgb(red = 139, green =  90, blue =   0, alpha = 120, maxColorValue = 255),
         'C>G' = rgb(red = 127, green =   0, blue = 255, alpha = 120, maxColorValue = 255),
         'G>C' = rgb(red = 127, green =   0, blue = 255, alpha = 120, maxColorValue = 255),
         'C>T' = rgb(red = 255, green = 215, blue =   0, alpha = 120, maxColorValue = 255),
         'G>A' = rgb(red = 255, green = 215, blue =   0, alpha = 120, maxColorValue = 255)
      )  
      plot(x = mut.tab$n.base, y = mut.tab$F, 
           ylab = "Mutant allele frequency", las = 1, pch = 19,
           col = c(mutation.colors, 'NA' = NA)[as.character(mut.tab$mutation)],
           ylim = c(min(mut.tab$F, na.rm = TRUE), 1), xlim = xlim)
      unique.colors <- unique(mutation.colors)
      labels <- sapply(unique.colors, function(a) paste(names(mutation.colors)[mutation.colors == a], collapse = ", "))
      #legend("topleft", legend = labels, fill = unique.colors, border = NA, bty = "n")
      legend(y = "center", x  = "right", legend = labels,
             inset = c(-25 * strwidth("a", units = 'figure'), 0),
             fill = unique.colors, border = NA, bty = "n")
      if (!is.null(segments)){
         if (vlines) {
            abline(v = segments$end.pos, lwd = 0.9, lty = 2)
         }   
         if (!is.null(data.model)) {
            for (i in 1:nrow(segments)) {
                segments(x0 = segments$start.pos[i], x1 = segments$end.pos[i], 
                         y0 = unique(data.model$muf$mufreqs[data.model$muf$CNt == segments$CNt[i]]), lwd = 0.4, lty = "24")
            }
         }
      }

   }
   plotWindows(baf.windows, ylab = "B allele frequency", 
               xlim = xlim, ylim = c(0, 0.5), las = 1,
               n.min = min.N.baf)
   if (!is.null(segments)){
      if (vlines) {
         abline(v = segments$end.pos, lwd = 0.9, lty = 2)
      }
      segments(x0 = segments$start.pos, y0 = segments$Bf, x1=segments$end.pos, y1 = segments$Bf, col = "red", lwd = 3)
      if (!is.null(data.model)) {
         for (i in 1:nrow(segments)) {
            segments(x0 = segments$start.pos[i], x1 = segments$end.pos[i], 
                     y0 = unique(data.model$baf$BAF[data.model$baf$CNt == segments$CNt[i]]), lwd = 0.4, lty = "24")
         }
      }
   }
   plotWindows(ratio.windows, ylab = "Depth ratio", 
               las = 1, n.min = min.N.ratio, ylim = c(0, 2.5))
   if (!is.null(segments)){
      if (vlines) {
         abline(v = segments$end.pos, lwd = 0.9, lty = 2)
      }   
      segments(x0 = segments$start.pos, y0 = segments$depth.ratio, x1=segments$end.pos, y1 = segments$depth.ratio, col = "red", lwd = 3)
      if (!is.null(data.model)) {
         ratios.theoric <- unique(data.model$muf[,c('CNt', 'depth.ratio')])

         segments(x0 = rep(min(segments$start.pos, na.rm =TRUE), times = nrow(ratios.theoric)),
                  x1 = rep(max(segments$end.pos, na.rm = TRUE), times = nrow(ratios.theoric)), 
                  y0 = ratios.theoric$depth.ratio, lwd = 0.4, lty = "24")
         text(x = rep(min(segments$start.pos, na.rm =TRUE), times = nrow(ratios.theoric)),
              y = ratios.theoric$depth.ratio, labels = ratios.theoric$CNt, pos = 2, offset = 0.5, cex = 0.8)
         
      }
   }
   mtext(main, 3, outer = TRUE, cex = par("cex.main"), line = 2)
   mtext("Position (Mb)", 1, outer = TRUE, cex = par("cex.main"), line = 3)

   mtext(at = seq(xlim[1], xlim[2], by = 1e7), text = round(seq(xlim[1]/1e6, xlim[2]/1e6, by = 10), 0), side = 1, cex = 0.6)
}

genome.view <- function(baf.windows, ratio.windows, segments = NULL, main = "", 
                            min.N.baf = 1, min.N.ratio = 1e4, CNr = rep(2, length(ratio.windows)),
                            cellularity = NULL, dna.index = NULL, avg.depth.ratio = NULL) {
   chr.metrics <- list()
   for (i in 1:length(ratio.windows)) {
      chr.metrics[[i]] <- range(ratio.windows[[i]]$mean, na.rm = TRUE)
   }
   chr.metrics <- do.call(rbind, chr.metrics)
   x0 <- chr.metrics[1,1]
   
}
