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
   gc.values   <- sort(unique(gc))
   gc.stats  <- list()
   pb        <- txtProgressBar(min = 1, max = length(gc.values), style = 3)
   for (ii in 1:length(gc.values)) {
      selected        <- gc == gc.values[ii]
      ratio.values    <-  ratio[selected]
      gc.stats$raw[[ii]] <- quantile(x = ratio.values , probs = c(0.25, 0.5, 0.75))
      ratio.values    <- ratio.values / median(ratio.values)
      gc.stats$adj[[ii]] <- quantile(x = ratio.values , probs = c(0.25, 0.5, 0.75))
      ratio[selected] <- ratio.values
      setTxtProgressBar(pb, ii)
   }
   cat("\n")
   gc.stats$raw <- do.call(rbind, gc.stats$raw)
   gc.stats$adj <- do.call(rbind, gc.stats$adj)   
   gc.stats$gc.values <- gc.values
   gc.stats$ratio <- ratio
   gc.stats
}

gc.sample.stats <- function (filename, gz = TRUE) {
   cat.command = paste("cat", filename, sep = " ")
   if (gz == TRUE) {
      cat.command = paste("gzcat", filename, sep = " ")
   }
   gc.data   <- read.table(pipe(paste(cat.command, "| awk '{print $6,$10}'", sep = " ")), header = TRUE)
   gc.values   <- sort(unique(gc.data$GC.percent))
   gc.stats  <- list()
   pb        <- txtProgressBar(min = 1, max = length(gc.values), style = 3)
   for (ii in 1:length(gc.values)) {
      selected        <- gc.data$GC.percent == gc.values[ii]
      ratio.values    <- gc.data$depth.ratio[selected]
      gc.stats$raw[[ii]] <- quantile(x = ratio.values , probs = c(0.25, 0.5, 0.75))
      ratio.values    <- ratio.values / median(ratio.values)
      gc.stats$adj[[ii]] <- quantile(x = ratio.values , probs = c(0.25, 0.5, 0.75))
      setTxtProgressBar(pb, ii)
   }
   cat("\n")
   gc.stats$raw <- do.call(rbind, gc.stats$raw)
   gc.stats$adj <- do.call(rbind, gc.stats$adj)
   gc.stats$gc.values <- gc.values
   gc.stats
}

windowValues <- function(x, positions, chromosomes, window = 1e6, overlap = 0,
                          weight = rep.int( x = 1, times = length(x)), start.coord = NULL) {
   chr.list         <- unique(chromosomes)
   results          <- list()
   results$windows  <- list()
   weight <- sqrt(weight)
   # x.adj            <- x
   for (i in 1:length(chr.list)) {
      results$windows[[i]] <- list()
      is.chr.i              <- chromosomes == chr.list[i]
      range.pos             <- range(positions[is.chr.i], na.rm = TRUE)
      if (!is.null(start.coord)) range.pos[1] <- as.integer(start.coord)
      beam.coords           <- seq(range.pos[1], range.pos[2], by = window)
      if (max(beam.coords) != range.pos[2] ) {
         beam.coords <- c(beam.coords, range.pos[2])
      }
      pos.i <- positions[is.chr.i]
      rat.i <- x[is.chr.i]
      #rat.o <- rat.i
      wgt.i <- weight[is.chr.i]
      if (overlap > 0 ) {
         segs        <- do.call(rbind, lapply(X = 1:(length(beam.coords) -1), FUN = function(x) c(beam.coords[x], beam.coords[x + 1])))
         segs.middle <- apply(X = segs, 2, FUN = function(x) x - overlap)
         segs.middle <- rbind(segs.middle[-1,], c(max(segs.middle, na.rm = TRUE), max(segs, na.rm = TRUE)))
         beam.coords <- sort(unique(c(segs, segs.middle)))
         cat(paste("chromosome:", chr.list[i], "from:", range.pos[1],"to:", range.pos[2], "window:", window, "overlaps", overlap, "\n",sep=" "))
      } else {
         cat(paste("chromosome:", chr.list[i], "from:", range.pos[1],"to:", range.pos[2], "window:", window,"\n",sep=" "))
      }
      if (length(beam.coords) > 2) {
         pb      <- txtProgressBar(min = 1, max = length(beam.coords) - 1, style = 3)
         progressb <- TRUE
      } else {
         progressb <- FALSE
      }
      for (n in 1:(length(beam.coords) - 1)) {
         if (overlap <= 0) {
            selected <- pos.i >= beam.coords[n] & pos.i < beam.coords[n + 1] 
            # rat.i.quartiles <- quantile(x = rep.int(x = rat.i[selected], times = wgt.i[selected]),
            #                             probs = c(0.25, 0.5, 0.75))
            rat.i.quartiles <- quantile(x = rat.i[selected], probs = c(0.25, 0.75))
            rat.i.mean      <- weighted.mean(x = rat.i[selected], w = wgt.i[selected])            
            results$windows[[i]][[n]] <- c(start = beam.coords[n],
                                           end = beam.coords[n + 1], mean = rat.i.mean,
                                           q = rat.i.quartiles[1], q = rat.i.quartiles[2],
                                           N = length(selected[selected == TRUE]))
         } else {
            if (n < (length(beam.coords) - 2)) {
               next.coord <- n + 2
            } else {
               next.coord <- n + 1
            }
            selected <- pos.i >= beam.coords[n] & pos.i < beam.coords[next.coord]
            # rat.i.quartiles <- quantile(x = rep.int(x = rat.i[selected], times = wgt.i[selected]),
            #                             probs = c(0.25, 0.5, 0.75))
            rat.i.quartiles <- quantile(x = rat.i[selected], probs = c(0.25, 0.75))
            rat.i.mean      <- weighted.mean(x = rat.i[selected], w = wgt.i[selected])            
            results$windows[[i]][[n]] <- c(start = beam.coords[n],
                                            end = beam.coords[next.coord], mean = rat.i.mean,
                                            q = rat.i.quartiles[1], q = rat.i.quartiles[2],
                                            N = length(selected[selected == TRUE]))
         }
         # rat.o[selected] <- as.numeric(rat.i.quartiles[2])
         if (progressb) {
            setTxtProgressBar(pb, n)
         }
      }
      results$windows[[i]] <- as.data.frame(do.call(rbind, results$windows[[i]]))
      # x.adj[is.chr.i]  <- rat.o
      cat("\n")
   }
   names(results$windows) <- chr.list
   # results$points <- x.adj
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
   L.max    <- which.max(mat[,2]) 
   sum.all  <- log1p(sum(exp(mat[-L.max, 2]-mat[L.max, 2]))) + mat[L.max, 2] 
   exp.vals <- 2^(mat[, 2] - sum.all)
   exp.vals[is.infinite(exp.vals)] <- 0
   values   <- cbind(mat[,1], expL = exp.vals)
   val.95   <- sum(values[, 2]) * interval
   sorted.v <- values[order(values[, 2], decreasing = TRUE),]
   sorted.v <- cbind(sorted.v,incr.col(sorted.v[, 2]))
   sorted.v <- sorted.v[sorted.v[,3] <= val.95, ]
   up.v     <- max(sorted.v[,1])
   low.v    <- min(sorted.v[,1]) 
   max.l    <- sorted.v[1,1]
   results$values  <- values
   results$confint <- c(low.v,up.v)
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
   data.frame(base.count = as.integer(n.base.mut), maj.base.freq = as.numeric(max.freqs[,1]),
              base = as.character(max.freqs[,2]), freq = as.numeric(max.freqs[,3]))
}

mutation.table <- function(abf.tab, mufreq.treshold = 0.15, min.reads = 40, max.mut.types = 3,
                           min.type.freq = 0.9, segments = NULL) {
   hom.filt    <- abf.tab$ref.zigosity == 'hom'
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
   if (length(grep("chr", abf.baf$chromosome)) == 0) {
      allele.seg[,c('chrom', 'start.pos', 'end.pos')]
   } else {
      allele.seg$chrom <- paste("chr", allele.seg$chrom, sep = '')
      allele.seg[,c('chrom', 'start.pos', 'end.pos')]
   }
}

segment.breaks <- function(abf.tab, breaks) {
   segments <- list()
   for (i in 1:nrow(breaks)) {
      # pos.filt <- abf.tab$chromosome == breaks$chrom[i] & abf.tab$n.base >= breaks$start.pos[i] & abf.tab$n.base <= breaks$end.pos[i]
      data.i  <- abf.tab[abf.tab$chromosome == breaks$chrom[i] & abf.tab$n.base >= breaks$start.pos[i] & abf.tab$n.base <= breaks$end.pos[i], ]
      het.i   <- data.i[data.i$ref.zigosity == 'het',]
      Bf.i    <- weighted.mean(x = het.i$Bf, w = sqrt(het.i$good.s.reads))
      ratio.i <- weighted.mean(x = data.i$adjusted.ratio, w = sqrt(data.i$depth.sample))
      segments[[i]] <- data.frame(chromosome  = breaks$chrom[i],
                                  start.pos   = breaks$start.pos[i],
                                  end.pos     = breaks$end.pos[i],
                                  Bf          = Bf.i,
                                  depth.ratio = ratio.i,
                                  N.BAF       = nrow(het.i),
                                  N.ratio     = nrow(data.i))
   }
   do.call(rbind, segments)
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
#       het.filter.i <- pos.filter.i & x$ref.zigosity == 'het'
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
