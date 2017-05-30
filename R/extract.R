sequenza.extract <- function(file, gz = TRUE, window = 1e6, overlap = 1,
    gamma = 80, kmin = 10, gamma.pcf = 140, kmin.pcf = 40,
    mufreq.treshold = 0.10, min.reads = 40, min.reads.normal = 10,
    min.reads.baf = 1, max.mut.types = 1, min.type.freq = 0.9,
    min.fw.freq = 0, verbose = TRUE, chromosome.list = NULL,
    breaks = NULL, breaks.method = "het", assembly = "hg19",
    weighted.mean = TRUE, normalization.method = "mean",
    ignore.normal = FALSE, parallel = 2L, gc.stats = NULL){

    if (is.null(gc.stats)) {
        gc.stats <- gc.sample.stats(file, gzip = gz, verbose = verbose,
            parallel = parallel)
    }
    chr.vect <- as.character(gc.stats$file.metrics$chr)
    if (normalization.method != "mean") {
        gc.normal.vect  <- median_gc(gc.stats$normal)
        gc.tumor.vect  <- median_gc(gc.stats$tumor)
    } else {
        gc.normal.vect  <- mean_gc(gc.stats$normal)
        gc.tumor.vect  <- mean_gc(gc.stats$tumor)
    }
    windows.baf   <- list()
    windows.ratio <- list()
    windows.normal <- list()
    windows.tumor <- list()
    mutation.list <- list()
    segments.list <- list()
    if (is.null(dim(breaks))) {
        breaks <- NULL
    }
    if (is.null(chromosome.list)) {
        chromosome.list <- chr.vect
    } else {
        chromosome.list <- chromosome.list[chromosome.list %in% chr.vect]
    }
    for (chr in chromosome.list){
        if (verbose){
            message("Processing ", chr, ": ", appendLF = FALSE)
        }
        tbi <- file.exists(paste(file, "tbi", sep = "."))
        if (tbi) {
            seqz.data   <- read.seqz(file, gzip = gz, chr_name = chr)
        } else {
            file.lines <- gc.stats$file.metrics[which(chr.vect == chr), ]
            seqz.data   <- read.seqz(file, gzip = gz,
                n_lines = c(file.lines$start, file.lines$end))
        }
        if (ignore.normal) {
            seqz.data$adjusted.ratio <- round(
                (seqz.data$depth.tumor /
                    gc.tumor.vect[as.character(seqz.data$GC.percent)])
                , 3)
        } else {
            seqz.data$adjusted.ratio <- round(
                (seqz.data$depth.tumor /
                    gc.tumor.vect[as.character(seqz.data$GC.percent)]) /
                (seqz.data$depth.normal /
                    gc.normal.vect[as.character(seqz.data$GC.percent)])
                , 3)
        }
        seqz.r.win <- windowValues(x = seqz.data$adjusted.ratio,
            positions = seqz.data$position,
            chromosomes = seqz.data$chromosome,
            window = window, overlap = overlap,
            weight = seqz.data$depth.normal)
        seqz.n.win <- windowValues(x = seqz.data$depth.normal,
            positions = seqz.data$position,
            chromosomes = seqz.data$chromosome,
            window = window, overlap = overlap)
        seqz.t.win <- windowValues(x = seqz.data$depth.tumor,
            positions = seqz.data$position,
            chromosomes = seqz.data$chromosome,
            window = window, overlap = overlap)

        seqz.hom <- seqz.data$zygosity.normal == "hom"
        seqz.het <- seqz.data[!seqz.hom, ]
        het.filt <- seqz.het$good.reads >= min.reads.baf
        het_ok <- nrow(seqz.het) > 0
        if (is.null(breaks)) {
            breaks_chr <- NULL
        } else {
            breaks_chr <- breaks[breaks$chrom == chr, ]
        }
        if (het_ok) {
            seqz.b.win <- windowBf(Af = seqz.het$Af, Bf = seqz.het$Bf,
                good.reads = seqz.het$good.reads,
                chromosomes = seqz.het$chromosome,
                positions = seqz.het$position, conf = 0.95,
                window = window, overlap = overlap)
        } else {
            seqz.b.win <- list()
            seqz.b.win[[1]] <- data.frame(start = min(seqz.data$position,
                na.rm = TRUE), end = max(seqz.data$position, na.rm = TRUE),
                mean = 0, q0 = 0,  q1 = 0, N = 1)
        }
        if (het_ok) {
            breaks_chr <- extract_breaks(data = seqz.data, data_het = seqz.het,
                ratio = seqz.r.win, baf = seqz.b.win,
                gamma = gamma, kmin = kmin, breaks = breaks_chr,
                gamma.pcf = gamma.pcf, kmin.pcf = kmin.pcf,
                assembly = assembly, chromosome = chr,
                method = breaks.method)
        } else {
            if (breaks.method == "full") {
                breaks_chr <- extract_breaks(data = seqz.data,
                    data_het = seqz.het, ratio = seqz.r.win, baf = seqz.b.win,
                    gamma = gamma, kmin = kmin,
                    gamma.pcf = gamma.pcf, kmin.pcf = kmin.pcf,
                    assembly = assembly, chromosome = chr,
                    method = breaks.method)
            }
        }
        if (is.null(breaks_chr) || nrow(breaks_chr) == 0 ||
            length(breaks_chr) == 0) {
            breaks_chr <- data.frame(chrom = chr,
                start.pos = min(seqz.data$position, na.rm = TRUE),
                end.pos = max(seqz.data$position, na.rm = TRUE))
        }
        seg.s1 <- segment.breaks(seqz.tab = seqz.data, breaks = breaks_chr,
            min.reads.baf = min.reads.baf, weighted.mean = weighted.mean)

        mut.tab   <- mutation.table(seqz.data,
            mufreq.treshold = mufreq.treshold,
            min.reads = min.reads, min.reads.normal = min.reads.normal,
            max.mut.types = max.mut.types, min.type.freq = min.type.freq,
            min.fw.freq = min.fw.freq, segments = seg.s1)

        windows.ratio[[which(chromosome.list == chr)]] <- seqz.r.win[[1]]
        windows.normal[[which(chromosome.list == chr)]] <- seqz.n.win[[1]]
        windows.tumor[[which(chromosome.list == chr)]] <- seqz.t.win[[1]]
        windows.baf[[which(chromosome.list == chr)]]   <- seqz.b.win[[1]]
        segments.list[[which(chromosome.list == chr)]] <- seg.s1
        mutation.list[[which(chromosome.list == chr)]] <- mut.tab
        if (verbose){
            message(nrow(mut.tab), ' variant calls; ',
                nrow(seqz.het), ' heterozygous positions; ',
                sum(seqz.hom), ' homozygous positions.')
        }
    }
    names(windows.baf)   <- chromosome.list
    names(windows.ratio) <- chromosome.list
    names(windows.normal) <- chromosome.list
    names(windows.tumor) <- chromosome.list
    names(mutation.list) <- chromosome.list
    names(segments.list) <- chromosome.list
    return(list(BAF = windows.baf, ratio = windows.ratio,
        normal = windows.normal, tumor = windows.tumor,
        mutations = mutation.list, segments = segments.list,
        chromosomes = chromosome.list, gc = gc.stats,
                avg.depth = 1))
}
