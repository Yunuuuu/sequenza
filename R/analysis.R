read.abfreq <- function (file, nrows = -1, fast = TRUE, gz = TRUE, 
    colClasses = c('factor', 'integer', 'factor', 'integer', 
      'integer', 'numeric', 'numeric', 'numeric', 'factor', 
      'numeric', 'numeric', "factor", "factor"), chr.name = NULL, ...) {
  if(fast && nrows == -1) {
    if(gz) {
       if (!is.null(chr.name)) {
          wc <- system(paste(paste('zgrep -c "^', chr.name, '\t"', sep = ''), file, sep = ' '), intern = TRUE)
       } else {
          wc <- system(paste('gunzip -c', file, '| wc'), intern = TRUE)
       }
    } else {
       if (!is.null(chr.name)) {
          wc <- system(paste(paste('grep -c "^', chr.name, '\t"', sep = ''), file, sep = ' '), intern = TRUE)
       } else {
          wc <- system(paste('wc', file), intern = TRUE)
       }
    }
    if (is.null(chr.name)) {
       wc <- sub("^ +", "", wc)
       wc <- strsplit(wc, ' ')[[1]][1]
    }
    nrows <- max(as.integer(wc), 1)
    message('Reading ', nrows, ' lines...')
  }
   if (!is.null(chr.name)) {
      if (gz) {
            grep.part <- paste("zgrep '^", chr.name, "\t'", sep = "")
         } else {
            grep.part <- paste("grep '^", chr.name, "\t'", sep = "")
         }
      abf.data   <- read.delim(pipe(paste(grep.part, file, sep = " ")), nrows = nrows, colClasses = colClasses, ...)
      head       <- colnames(read.table(file, header = TRUE, nrow = 1 ))
      colnames(abf.data) <- head
      abf.data
   } else { 
      read.delim(file, nrows = nrows, colClasses = colClasses, ...)
   }
}

read.acgt <- function (file, nrows = -1, fast = TRUE, gz = TRUE, 
    colClasses = c('factor', 'integer', 'factor', 'integer', 
      'integer', 'integer', 'integer', 'integer'), ...) {
  if(fast && nrows == -1) {
    if(gz) {
      wc <- system(paste('gunzip -c', file, '| wc'), intern = TRUE)
    } else {
      wc <- system(paste('wc', file), intern = TRUE)
    }
    wc <- sub("^ +", "", wc)
    wc <- strsplit(wc, ' ')[[1]][1]
    nrows <- max(as.integer(wc), 1)
    message('Reading ', nrows, ' lines...')
  } 
  read.delim(file, nrows = nrows, colClasses = colClasses, ...)
}

gc.norm <- function (ratio, gc) {
   dr.by.gc <- split(ratio, gc)
   raw <- t(sapply(dr.by.gc, quantile, probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
   dr.by.gc.median <- sapply(dr.by.gc, median)
   adj <- sweep(raw, 1, dr.by.gc.median, '/')
   adjratio = ratio / unsplit(dr.by.gc.median, gc)
   list(raw = raw, adj = adj, gc.values = as.numeric(names(dr.by.gc)), ratio = adjratio)
}

gc.sample.stats <- function (filename, gz = TRUE) {
   cat.command = paste("cat", filename)
   if (gz) {
      cat.command = paste("gunzip -c", filename)
   }
   gc.data   <- read.table(pipe(paste(cat.command, "| awk '{print $6,$10}'")), header = TRUE)
   dr.by.gc <- split(gc.data$depth.ratio, gc.data$GC.percent)
   raw <- t(sapply(dr.by.gc, quantile, probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
   dr.by.gc.median <- sapply(dr.by.gc, median)
   adj <- sweep(raw, 1, dr.by.gc.median, '/')
   list(raw = raw, adj = adj, gc.values = as.numeric(names(dr.by.gc)))
}

windowValues <- function(x, positions, chromosomes, window = 1e6, overlap = 0, verbose = TRUE,
                          weight = rep.int( x = 1, times = length(x)), start.coord = 1) {
   results          <- list()
   results$windows  <- list()
   weight   <- sqrt(weight)
   xw       <- x * weight
   overlap  <- as.integer(overlap)
   window.o <- window - round(window * (overlap / (overlap + 1)), 0)
   mat.w    <- data.frame(chr = chromosomes, pos = positions, x = x, xw = x * weight,
                          w = weight, stringsAsFactors = FALSE)
   mat.w    <- split(mat.w, mat.w$chr)
   do.windows <- function(x.i, w.i, xw.i, breaks, overlap){
      coords    <- data.frame(start = breaks[-length(breaks)], end = breaks[-1])
      quartiles <- do.call(rbind, 
                           lapply(X = x.i, FUN = function(x) quantile(x, probs = c(0.25, 0.75),
                                                                      na.rm = TRUE)))
      sum.w     <- sapply(X = w.i, FUN = function(x) sum(x, na.rm = TRUE))
      sum.xw    <- sapply(X = xw.i, FUN = function(x) sum(x, na.rm = TRUE))
      size      <- sapply(X = x.i, FUN = length)
      data.frame(coords, mean = sum.xw/sum.w, q0 = quartiles[,1],
                 q1 = quartiles[,2], N = size, row.names = 1:length(size))            
   }
   for (i in 1:length(mat.w)) {
      results$windows[[i]] <- list()
      range.pos            <- range(mat.w[[i]]$pos, na.rm = TRUE)
      if (!is.null(start.coord)) {
         range.pos[1] <- as.integer(start.coord)
      }
      beam.coords          <- seq(range.pos[1], range.pos[2], by = window.o)
      if (max(beam.coords) != range.pos[2] ) {
         beam.coords <- c(beam.coords, range.pos[2])
      }
      f.windows <- cut(x = mat.w[[i]]$pos, breaks = beam.coords)
      xw      <- split(x = mat.w[[i]]$xw, f = f.windows)
      w       <- split(x = mat.w[[i]]$w, f = f.windows)
      x       <- split(x = mat.w[[i]]$x, f = f.windows)
      if (overlap > 0 ) {
         if (verbose) {
            cat(paste("chromosome:", names(mat.w)[i], "from:", range.pos[1],"to:", range.pos[2], "window:", window, "overlapping windows:", overlap, "\n",sep=" "))
         }
      } else {
         if (verbose) {
            cat(paste("chromosome:", names(mat.w)[i], "from:", range.pos[1],"to:", range.pos[2], "window:", window,"\n",sep=" "))
         }
      }
      
      results$windows[[i]] <- do.windows(x.i = x, w.i = w, xw.i = xw,
                                         breaks = beam.coords, overlap = overlap)
   }
   names(results$windows) <- names(mat.w)
   results
}

get.ci <- function(mat, interval = 0.95) {
   incr.col <- function(x) {
      for (i in 1:length(x)) {
         if (i == 1) {
            x[i] <- x[i]
         } else {
            x[i] <- x[i] + x[i - 1]
         }
      }
      x
   }
   results  <- list()
   values   <- sort(unique(mat[, 1]))
   mat      <- t(mapply(x = values, FUN = function(x ,mat) cbind(x, max(mat[mat[, 1] == x, 2])), MoreArgs = list(mat = mat)))
   L.max    <- which.max(mat[, 2]) 
   sum.all  <- log1p(sum(exp(mat[-L.max, 2]-mat[L.max, 2]))) + mat[L.max, 2] 
   exp.vals <- 2^(mat[, 2] - sum.all)
   exp.vals[is.infinite(exp.vals)] <- 0
   values   <- cbind(mat[, 1], expL = exp.vals)
   val.95   <- quantile(values[, 2], prob = interval, na.rm = TRUE)
   values.s <- values[values[, 2] >= val.95, ]
   if (is.null(dim(values.s))) {
      up.v  <- values.s[1]
      low.v <- values.s[1]
      max.l <- values.s[1]
   } else {
      up.v     <- max(values.s[, 1])
      low.v    <- min(values.s[, 1]) 
      max.l    <- values[which.max(values[, 2]), 1]
   }
   results$values  <- values
   results$confint <- c(low.v, up.v)
   results$max.l   <- max.l
   results
}

# merge.baf.ratio <- function(baf.segments, ratio.segments) {
#    baf.table    <- lapply(1:length(baf.segments), FUN = function(x) cbind(chromosome = names(baf.segments)[x], 
#                                                                         as.data.frame(do.call(rbind, baf.segments[[x]]))))
#    ratio.table  <- lapply(1:length(ratio.segments), FUN = function(x) cbind(chromosome = names(ratio.segments)[x],
#                                                                           as.data.frame(do.call(rbind, ratio.segments[[x]]))))
#    baf.table    <- do.call(rbind, baf.table)
#    ratio.table  <- do.call(rbind, ratio.table)
#    baf.table$chromosome <- as.character(baf.table$chromosome)
#    ratio.table$chromosome <- as.character(ratio.table$chromosome)
#    baf.index    <- sapply(1:nrow(baf.table), FUN = function(x) paste(baf.table[x, 1:2], collapse ="_"))
#    ratio.index  <- sapply(1:nrow(ratio.table), FUN = function(x) paste(ratio.table[x, 1:2], collapse = "_"))
#    baf.table    <- cbind(index = baf.index, baf.table[, c(5,7)])
#    ratio.table  <- cbind(index = ratio.index, ratio.table[, c(1,2,3,5,7)])
#    colnames(ratio.table) <- c("index", "chromosome", "start", "end", "ratio", "N.depth")
#    colnames(baf.table) <- c("index", "Bf", "N.BAF")
#    merged.table <- merge(ratio.table, baf.table, all = TRUE, sort = FALSE)
#    merged.table <- merged.table[, -1]

#    merged.table
# }


mut.fractions <- function(AB.sample, Af) {
  F = 1 - Af
   base.mut <- lapply(X = AB.sample, FUN = function(x) unlist(strsplit(as.character(x), split = '[:]')))
   frequencify <- function (x) {
      base.name <- substr(unlist(x), 1, 1)
      base.val  <- as.numeric(substr(unlist(x), 2, nchar(x)))
      setNames(base.val, base.name)
   }
   base.freqs <- lapply(X = base.mut, FUN = frequencify)   
   n.base.mut <- do.call(c, lapply(X = base.mut, FUN = length))
   max.fq <- function (x) {
      freq.rel <- base.freqs[[x]] / F[x]
      f.max    <- which.max(freq.rel)
      c(freq.rel[f.max], names(base.freqs[[x]])[f.max], base.freqs[[x]][f.max])
   }
   max.freqs  <- do.call(rbind, lapply(1:length(F), max.fq))
   data.frame(base.count = as.integer(n.base.mut), maj.base.freq = as.numeric(max.freqs[, 1]),
              base = as.character(max.freqs[,2]), freq = as.numeric(max.freqs[,3]))
}

mutation.table <- function(abf.tab, mufreq.treshold = 0.15, min.reads = 40, max.mut.types = 3,
                           min.type.freq = 0.9, segments = NULL) {
   hom.filt    <- abf.tab$ref.zygosity == 'hom'
   abf.tab     <- abf.tab[hom.filt, ]
   reads.filt  <- abf.tab$good.s.reads >= min.reads
   abf.tab     <- abf.tab[reads.filt, ]
   mufreq.filt <- abf.tab$Af <= (1 - mufreq.treshold)
   abf.tab     <- abf.tab[mufreq.filt, ]
   if (!is.null(segments)) {
      for (i in 1:nrow(segments)) {
         pos.filt <- abf.tab$chromosome == segments$chrom[i] & abf.tab$n.base >= segments$start.pos[i] & abf.tab$n.base <= segments$end.pos[i]
         abf.tab$adjusted.ratio[pos.filt] <- segments$depth.ratio[i]
      }
   }
   mu.fracts   <- mut.fractions(AB.sample = abf.tab$AB.sample, Af = abf.tab$Af)
   mufreq.filt <- mu.fracts$freq >= mufreq.treshold
   type.filt   <- mu.fracts$base.count <= max.mut.types
   prop.filt   <- mu.fracts$maj.base.freq <= min.type.freq
   mut.type    <- paste(abf.tab$AB.germline, mu.fracts$base, sep = '>')   
   abf.tab     <- abf.tab[,c('chromosome', 'n.base', 'GC.percent', 'good.s.reads', 'adjusted.ratio')]
   abf.tab     <- cbind(abf.tab, F = mu.fracts$freq, mutation = mut.type)
   abf.tab[mufreq.filt, ]
}

find.breaks <- function(abf.baf, gamma = 80, kmin = 10, baf.thres = c(0, 0.5), verbose = FALSE, ...) {
   require(copynumber)
   chromosome <- gsub(x = abf.baf$chromosome, pattern = "chr", replacement = "")
   logR = data.frame(chrom = chromosome, 
                     pos = abf.baf$n.base,
                     s1 = log2(abf.baf$adjusted.ratio))
   BAF = data.frame(chrom = chromosome, 
                    pos = abf.baf$n.base,
                    s1 = abf.baf$Bf)
   logR.wins <- winsorize(logR, verbose = verbose)
   allele.seg <- aspcf(logR = logR.wins, BAF = BAF, baf.thres = baf.thres,
                       verbose = verbose, gamma = gamma, kmin = kmin, ...)
    if (length(grep("chr", abf.baf$chromosome)) > 0) { 
        allele.seg$chrom <- paste("chr", allele.seg$chrom, sep = "")
    }
    allele.seg[, c("chrom", "start.pos", "end.pos")]
}

segment.breaks <- function(abf.tab, breaks) {
   w.r     <- sqrt(abf.tab$depth.sample)
   rw      <- abf.tab$adjusted.ratio * w.r
   w.b     <- sqrt(abf.tab$good.s.reads)
   bw      <- abf.tab$Bf * w.b
   abf.tab <- cbind(abf.tab[, c("chromosome", "n.base", "ref.zygosity")],
                    rw = rw, w.r = w.r, bw = bw, w.b = w.b)
   chromosomes <- unique(abf.tab$chromosome)
   segments <- list()
   for (i in 1:length(chromosomes)) {
      abf.i       <- abf.tab[abf.tab$chromosome == chromosomes[i], ]
      abf.b.i     <- abf.i[abf.i$ref.zygosity == 'het', ]
      breaks.i    <- breaks[breaks$chrom == chromosomes[i], ]
      nb          <- nrow(breaks.i)
      breaks.vect <- do.call(cbind, split.data.frame(breaks.i[,c("start.pos", "end.pos")], f = 1:nb))
      fact.r.i    <- cut(abf.i$n.base, breaks.vect)
      fact.b.i    <- cut(abf.b.i$n.base, breaks.vect)
      seg.i.s.r   <- sapply(X = split(abf.i$w.r, f = fact.r.i), FUN = length)
      seg.i.s.b   <- sapply(X = split(abf.b.i$w.b, f = fact.b.i), FUN = length)      
      seg.i.rw    <- sapply(X = split(abf.i$rw, f = fact.r.i), FUN = function(x) sum(x, na.rm = TRUE))
      seg.i.w.r   <- sapply(X = split(abf.i$w.r, f = fact.r.i), FUN = function(x) sum(x, na.rm = TRUE))
      seg.i.bw    <- sapply(X = split(abf.b.i$bw, f = fact.b.i), FUN = function(x) sum(x, na.rm = TRUE))
      seg.i.w.b   <- sapply(X = split(abf.b.i$w.b, f = fact.b.i), FUN = function(x) sum(x, na.rm = TRUE))
      segments.i <- data.frame(chromosome  = chromosomes[i], start.pos = as.numeric(breaks.vect[-length(breaks.vect)]),
                               end.pos = as.numeric(breaks.vect[-1]), Bf = seg.i.bw/seg.i.w.b, N.BAF = seg.i.s.b,
                               depth.ratio = seg.i.rw/seg.i.w.r, N.ratio = seg.i.s.r, stringsAsFactors = FALSE)
      segments[[i]] <- segments.i[seq(from = 1, to = nrow(segments.i), by = 2),]
   }
   segments <- do.call(rbind, segments)
   row.names(segments) <- 1:nrow(segments)
   segments
}

# segment.chromosome <- function(x, breaks) {
#    segments <- list()
#    chromosome = unique(x$chromosome)
#    for (i in 1:length(breaks)) {
#       if (i ==1) {
#          start.i = 1
#       } else {
#          start.i = breaks[i-1]
#       }
#       end.i = breaks[i]
#       pos.filter.i <- x$n.base >= start.i & x$n.base <= end.i
#       het.filter.i <- pos.filter.i & x$ref.zygosity == 'het'
#       start.j <- min(x$n.base[pos.filter.i])
#       # Bf.i    <- median(rep(x$Bf[het.filter.i], 
#       #                 times = round( x$depth.sample[het.filter.i] * x$sample.reads.above.quality[het.filter.i], 0)))
#       # ratio.i <- median(rep(x$adjusted.ratio[pos.filter.i], 
#       #                 times = x$depth.sample[pos.filter.i]))
#        Bf.i    <- weighted.mean(x = x$Bf[het.filter.i], w = sqrt(x$depth.sample[het.filter.i] * x$sample.reads.above.quality[het.filter.i]))
#        ratio.i <- weighted.mean(x = x$adjusted.ratio[pos.filter.i], w = sqrt(x$depth.sample[pos.filter.i]))
#       segments[[i]] <- data.frame(chromosome  = chromosome,
#                                   start       = start.j,
#                                   end         = end.i,
#                                   Bf          = Bf.i,
#                                   depth.ratio = ratio.i,
#                                   N.BAF       = length(het.filter.i[het.filter.i == TRUE]),
#                                   N.ratio     = length(pos.filter.i[pos.filter.i == TRUE]))
#    }
#    do.call(rbind, segments)
# }
