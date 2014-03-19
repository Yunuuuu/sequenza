read.seqz <- function (file, nrows = -1, fast = FALSE, gz = TRUE, header = TRUE,
    colClasses = c('character', 'integer', 'character', 'integer',
      'integer', 'numeric', 'numeric', 'numeric', 'character',
      'numeric', 'numeric', "character", "character"), chr.name = NULL, n.lines = NULL, ...) {
   if (!is.null(n.lines) & is.null(chr.name)) fast <-  FALSE
   if(fast && nrows == -1) {
    if(gz) {
       if (!is.null(chr.name)) {
          wc <- system(paste('gzip -d -c ',file,' | grep -c "^', chr.name, '\t"', sep = ''), intern = TRUE)
       } else {
          wc <- system(paste('gzip -d -c', file, '| wc'), intern = TRUE)
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
         grep.part <- paste("gzip -d -c ", file," | grep '^", chr.name, "\t'", sep = "")
      } else {
         grep.part <- paste("grep '^", chr.name, "\t'", sep = "")
      }
      abf.data   <- read.delim(pipe(grep.part), nrows = nrows, colClasses = colClasses, header = FALSE, ...)
      if (header == TRUE) {
         head       <- colnames(read.table(file, header = TRUE, nrows = 1 ))
         colnames(abf.data) <- head
      }
      abf.data
   } else {
      if (!is.null(n.lines)){
         if (!is.numeric(n.lines) | length(n.lines) != 2) stop("n.lines must be a vector of 2 integers")
         n.lines <- round(sort(n.lines), 0)
         if (header == TRUE) {
            n.lines <- n.lines + 1
         }
         if(gz) {
            abf.data <- read.delim(pipe(paste("gzip -d -c", file,"| sed -n '", paste(n.lines[1], n.lines[2], sep = ","),"p'")),
                       colClasses = colClasses, nrows = 1 + n.lines[2] - n.lines[1], header = FALSE,...)
         }  else{
            abf.data <- read.delim(pipe(paste("sed -n '", paste(n.lines[1], n.lines[2], sep = ","),"p'", file)),
                       colClasses = colClasses, nrows = 1 + n.lines[2] - n.lines[1], header = FALSE, ...)
         }
         if (header == TRUE) {
            head  <- colnames(read.table(file, header = TRUE, nrows = 1 ))
            colnames(abf.data) <- head
         }
         abf.data
      } else {
         read.delim(file, nrows = nrows, colClasses = colClasses, header = header, ...)
      }
   }
}

read.acgt <- function (file, colClasses = c('character', 'integer', 'character', 'integer',
                                            'integer', 'integer', 'integer', 'integer'), ...) {
   read.seqz(file = file , colClasses = colClasses, ...)
}

gc.norm <- function (x, gc) {
   dr.by.gc <- split(x, gc)
   raw <- t(sapply(dr.by.gc, quantile, probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
   dr.by.gc.median <- sapply(dr.by.gc, median, na.rm = TRUE)
   dr.by.gc.mean <- sapply(dr.by.gc, mean, na.rm = TRUE)
   adj <- sweep(raw, 1, dr.by.gc.median, '/')
   list(raw = raw, adj = adj, gc.values = as.numeric(names(dr.by.gc)),
        raw.mean = dr.by.gc.mean, raw.median = dr.by.gc.median)
}

gc.sample.stats <- function (file, gz = TRUE) {
   colClasses = c('character', 'numeric', 'numeric')
   if (gz) {
      abf.data <- read.delim(pipe(paste('gzip -d -c', file, '| cut -f 1,6,10')), colClasses = colClasses)
   } else {
      abf.data <- read.delim(pipe(paste('cut -f 1,6,10', file)), colClasses = colClasses)
   }
   gc.stats <- gc.norm(x = abf.data$depth.ratio,
                       gc = abf.data$GC.percent)
   chr.ord  <- unique(abf.data$chromosome)
   chr.dim  <- lapply(X = split(abf.data$chromosome, abf.data$chromosome), FUN = length)
   chr.dim  <- data.frame(chr = chr.ord, n.lines = do.call(rbind,chr.dim[chr.ord]))
   chr.dim$start <- cumsum(c(1, chr.dim$n.lines[-length(chr.dim$n.lines)]))
   chr.dim$end   <- chr.dim$start + chr.dim$n.lines - 1
   gc.stats$file.metrics <- chr.dim
   gc.stats
}

windowValues <- function(x, positions, chromosomes, window = 1e6, overlap = 0,
                         weight = rep.int( x = 1, times = length(x)), start.coord = 1) {
  weight   <- sqrt(weight)
  overlap  <- as.integer(overlap)
  window.offset <- window - round(window * (overlap / (overlap + 1)))
  chr.ordered <- unique(chromosomes)
  data.splitByChr    <- split(data.frame(pos = positions, x = x, weight = weight), 
                              f = factor(chromosomes, levels = chr.ordered))
  lapply(data.splitByChr, function(data.oneChr) {
    range.pos <- range(data.oneChr$pos, na.rm = TRUE)
    if (!is.null(start.coord)) {
      range.pos[1] <- as.integer(start.coord)
    }
    beam.coords <- seq(range.pos[1], range.pos[2], by = window.offset)
    if (max(beam.coords) != range.pos[2] ) {
      beam.coords <- c(beam.coords, range.pos[2])
    }
    nWindows <- length(beam.coords) - overlap - 1
    pos.cut <- cut(data.oneChr$pos, breaks = beam.coords)
    x.split <- split(data.oneChr$x, f = pos.cut)
    weight.split <- split(data.oneChr$weight, f = pos.cut)
    window.starts <- beam.coords[1:nWindows]
    window.ends <- beam.coords[(1:nWindows) + 1 + overlap]
    idx.list <- lapply(1:nWindows, function(ii) ii + (0:overlap))
    x.window <- lapply(idx.list, function(idx) unlist(x.split[idx], use.names = FALSE))
    weight.window <- lapply(idx.list, function(idx) unlist(weight.split[idx], use.names = FALSE))
    window.means <- mapply(weighted.mean, x = x.window, w = weight.window)
    window.quantiles <- sapply(x.window, quantile, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
    window.counts <- sapply(x.window, length)
    data.frame(start = window.starts, end = window.ends, mean = window.means, 
               q0 = window.quantiles[1,], q1 = window.quantiles[2,], N = window.counts)
  })
}

get.ci <- function(cp.table, interval = 0.95) {
  znormsort <- sort(cp.table$z, decreasing = TRUE)
  znormcumLik <- cumsum(znormsort)
  n <- sapply(interval, function(x) sum(znormcumLik < x) + 1)
  LikThresh <- znormsort[n]
  values.x <- data.frame(x = cp.table$x, y = apply(cp.table$z, 1, max))
  values.y <- data.frame(x = apply(cp.table$z, 2, max), y = cp.table$y)
  up.x  <- max(values.x$x[values.x$y >= LikThresh])
  low.x <- min(values.x$x[values.x$y >= LikThresh])
  max.x <- values.x$x[which.max(values.x$y)]
  up.y  <- max(values.y$y[values.y$x >= LikThresh])
  low.y <- min(values.y$y[values.y$x >= LikThresh])
  max.y <- values.y$y[which.max(values.y$x)]
  results <- list()
  results$values.x <- values.x
  results$confint.x <- c(low.x, up.x)
  results$max.x <- max.x
  results$values.y <- values.y
  results$confint.y <- c(low.y, up.y)
  results$max.y <- max.y
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


mut.fractions <- function(AB.sample, Af, sample.strand) {
  F = 1 - Af
   base.mut <- lapply(X = AB.sample, FUN = function(x) unlist(strsplit(as.character(x), split = '[:]')))
   base.fw  <- lapply(X = sample.strand, FUN = function(x) unlist(strsplit(as.character(x), split = '[:]')))
   frequencify <- function (x) {
      base.name <- substr(unlist(x), 1, 1)
      base.val  <- as.numeric(substr(unlist(x), 2, nchar(x)))
      setNames(base.val, base.name)
   }
   base.freqs <- lapply(X = base.mut, FUN = frequencify)
   fw.freqs   <- lapply(X = base.fw, FUN = frequencify)
   n.base.mut <- do.call(c, lapply(X = base.mut, FUN = length))
   max.fq <- function (x) {
      freq.rel <- base.freqs[[x]] / F[x]
      f.max    <- which.max(freq.rel)
      c(freq.rel[f.max], names(base.freqs[[x]])[f.max], base.freqs[[x]][f.max], fw.freqs[[x]][f.max])
   }
   max.freqs  <- do.call(rbind, lapply(1:length(F), max.fq))
   data.frame(base.count = as.integer(n.base.mut), maj.base.freq = as.numeric(max.freqs[, 1]),
              base = as.character(max.freqs[,2]), freq = as.numeric(max.freqs[,3]),
              fw.freq = as.numeric(max.freqs[,4]))
}

mutation.table <- function(abf.tab, mufreq.treshold = 0.15, min.reads = 40, min.reads.normal = 10,
                           max.mut.types = 3, min.type.freq = 0.9, min.fw.freq = 0,
                           segments = NULL) {
   chroms      <- unique(abf.tab$chromosome)
   hom.filt    <- abf.tab$ref.zygosity == 'hom'
   abf.tab     <- abf.tab[hom.filt, ]
   reads.filt  <- abf.tab$good.s.reads >= min.reads & abf.tab$depth.normal >= min.reads.normal
   abf.tab     <- abf.tab[reads.filt, ]
   mufreq.filt <- abf.tab$Af <= (1 - mufreq.treshold)
   abf.tab     <- abf.tab[mufreq.filt, ]
   if (!is.null(segments)) {
      for (i in 1:nrow(segments)) {
         pos.filt <- abf.tab$chromosome == segments$chrom[i] & abf.tab$n.base >= segments$start.pos[i] & abf.tab$n.base <= segments$end.pos[i]
         abf.tab$adjusted.ratio[pos.filt] <- segments$depth.ratio[i]
      }
   }
   abf.dummy   <- data.frame(chromosome = chroms, n.base = 1, GC.percent = NA, good.s.reads = NA,
                             adjusted.ratio = NA, F = 0, mutation = 'NA', stringsAsFactors= FALSE)
   if (nrow(abf.tab) >= 1) {
      mu.fracts   <- mut.fractions(AB.sample = abf.tab$AB.sample, Af = abf.tab$Af,
                                   sample.strand = abf.tab$sample.strand)
      mufreq.filt <- mu.fracts$freq >= mufreq.treshold
      type.filt   <- mu.fracts$base.count <= max.mut.types
      prop.filt   <- mu.fracts$maj.base.freq >= min.type.freq
      if (!is.na(min.fw.freq)) {
         fw.2 = 1 - min.fw.freq
         fw.2 <- sort(c(fw.2, min.fw.freq))
         fw.filt     <- mu.fracts$fw.freq > fw.2[1] & mu.fracts$fw.freq < fw.2[2]
         mufreq.filt <- mufreq.filt & type.filt  & prop.filt & fw.filt
      } else {
         mufreq.filt <- mufreq.filt & type.filt  & prop.filt
      }
      mut.type    <- paste(abf.tab$AB.germline, mu.fracts$base, sep = '>')
      abf.tab     <- abf.tab[,c('chromosome', 'n.base', 'GC.percent', 'good.s.reads', 'adjusted.ratio')]
      abf.tab     <- cbind(abf.tab, F = mu.fracts$freq, mutation = mut.type)
      rbind(abf.tab[mufreq.filt, ], abf.dummy)
   } else {
      abf.dummy
   }
}

find.breaks <- function(abf.baf, gamma = 80, kmin = 10, baf.thres = c(0, 0.5), verbose = FALSE, ...) {
   chromosome <- gsub(x = abf.baf$chromosome, pattern = "chr", replacement = "")
   logR = data.frame(chrom = chromosome,
                     pos = abf.baf$n.base,
                     s1 = log2(abf.baf$adjusted.ratio))
   BAF = data.frame(chrom = chromosome,
                    pos = abf.baf$n.base,
                    s1 = abf.baf$Bf)
   logR.wins <- copynumber::winsorize(logR, verbose = verbose)
   allele.seg <- copynumber::aspcf(logR = logR.wins, BAF = BAF, baf.thres = baf.thres,
                       verbose = verbose, gamma = gamma, kmin = kmin, ...)
    if (length(grep("chr", abf.baf$chromosome)) > 0) {
        allele.seg$chrom <- paste("chr", allele.seg$chrom, sep = "")
    }
    allele.seg[allele.seg$end.pos - allele.seg$start.pos != 0,
               c("chrom", "start.pos", "end.pos")]
}

segment.breaks <- function(abf.tab, breaks, weighted.mean = TRUE) {
   if (weighted.mean == TRUE){
      w.r     <- sqrt(abf.tab$depth.sample)
      rw      <- abf.tab$adjusted.ratio * w.r
      w.b     <- sqrt(abf.tab$good.s.reads)
      bw      <- abf.tab$Bf * w.b
      abf.tab <- cbind(abf.tab[, c("chromosome", "n.base", "ref.zygosity")],
                    rw = rw, w.r = w.r, bw = bw, w.b = w.b)
   }
   chr.order <- unique(abf.tab$chromosome)
   abf.tab <- split(abf.tab, f = abf.tab$chromosome)
   segments <- list()
   for (i in 1:length(abf.tab)) {
      abf.b.i     <- abf.tab[[i]][abf.tab[[i]]$ref.zygosity == 'het', ]
      breaks.i    <- breaks[breaks$chrom == names(abf.tab)[i], ]
      nb          <- nrow(breaks.i)
      breaks.vect <- do.call(cbind, split.data.frame(breaks.i[,c("start.pos", "end.pos")], f = 1:nb))
      fact.r.i    <- cut(abf.tab[[i]]$n.base, breaks.vect)
      fact.b.i    <- cut(abf.b.i$n.base, breaks.vect)
      seg.i.s.r   <- sapply(X = split(abf.tab[[i]]$chromosome, f = fact.r.i), FUN = length)
      seg.i.s.b   <- sapply(X = split(abf.b.i$chromosome, f = fact.b.i), FUN = length)
      if (weighted.mean == TRUE){
         seg.i.rw    <- sapply(X = split(abf.tab[[i]]$rw, f = fact.r.i), FUN = function(a) sum(a, na.rm = TRUE))
         seg.i.w.r   <- sapply(X = split(abf.tab[[i]]$w.r, f = fact.r.i), FUN = function(a) sum(a, na.rm = TRUE))
         seg.i.bw    <- sapply(X = split(abf.b.i$bw, f = fact.b.i), FUN = function(a) sum(a, na.rm = TRUE))
         seg.i.w.b   <- sapply(X = split(abf.b.i$w.b, f = fact.b.i), FUN = function(a) sum(a, na.rm = TRUE))
         segments.i <- data.frame(chromosome  = names(abf.tab)[i], start.pos = as.numeric(breaks.vect[-length(breaks.vect)]),
                               end.pos = as.numeric(breaks.vect[-1]), Bf = seg.i.bw/seg.i.w.b, N.BAF = seg.i.s.b,
                               depth.ratio = seg.i.rw/seg.i.w.r, N.ratio = seg.i.s.r, stringsAsFactors = FALSE)
      } else {
        seg.i.r   <- sapply(X = split(abf.tab[[i]]$adjusted.ratio, f = fact.r.i), FUN = mean)
        seg.i.b   <- sapply(X = split(abf.b.i$Bf, f = fact.b.i), FUN = mean)
        segments.i <- data.frame(chromosome  = names(abf.tab)[i], start.pos = as.numeric(breaks.vect[-length(breaks.vect)]),
                                 end.pos = as.numeric(breaks.vect[-1]), Bf = seg.i.b, N.BAF = seg.i.s.b,
                                 depth.ratio = seg.i.r, N.ratio = seg.i.s.r, stringsAsFactors = FALSE)
      }
      segments[[i]] <- segments.i[seq(from = 1, to = nrow(segments.i), by = 2),]
   }
   segments <- do.call(rbind, segments[as.factor(chr.order)])
   row.names(segments) <- 1:nrow(segments)
   segments
}
