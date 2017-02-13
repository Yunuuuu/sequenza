read.seqz <- function(file, n_lines = NULL, gzip = TRUE,
    col_types = "ciciidddcddccc", chr_name = NULL,
    buffer = 33554432, parallel = 2L,
    col_names = c("chromosome", "position", "base.ref", "depth.normal",
                  "depth.tumor", "depth.ratio", "Af", "Bf", "zygosity.normal",
                  "GC.percent", "good.reads", "AB.normal", "AB.tumor",
                  "tumor.strand"), ...) {

    if (is.null(n_lines)) {
        skip <- 1
        n_max <- Inf
    } else {
        n_lines <- round(sort(n_lines), 0)
        skip <- n_lines[1]
        if (skip > 1) {
            skip <- skip - 1
        }
        n_max <- n_lines[2] - skip
    }
    if (!is.null(chr_name)) {
        tbi <- file.exists(paste(file, "tbi", sep = "."))
        if (tbi) {
            read.seqz.tbi(file, chr_name, col_names, col_types)
        } else {
            read.seqz.chr(file, chr_name = chr_name, col_types = col_types,
                col_names = col_names, gzip = gzip,
                buffer = buffer, parallel = parallel)
        }
    } else {
        readr::read_tsv(file = file, col_types = col_types, skip = skip,
            n_max = n_max, col_names = col_names, ...)
    }
}

gc.sample.stats <- function (file, col_types = "c----d---d----", ...) {
    seqz_data <- readr::read_tsv(file, col_types = col_types, ...)
    #gc_stats <- gc_norm(x  = seqz_data[, 2],
    #                   gc = seqz_data[, 3])
    chr_ord  <- unique(seqz_data$chromosome)
    chr_dim  <- lapply(X = split(seqz_data$chromosome,
        seqz_data$chromosome), FUN = length)
    chr_dim  <- data.frame(chr = chr_ord,
        n_lines = do.call(rbind, chr_dim[chr_ord]))
    chr_dim$start <- cumsum(c(1, chr_dim$n_lines[-length(chr_dim$n_lines)]))
    chr_dim$end   <- chr_dim$start + chr_dim$n_lines - 1
    #gc_stats$file_metrics <- chr_dim
    #gc_stats
    chr_dim
}

read.seqz.chr <- function(file, chr_name, col_names, col_types,
    gzip, buffer, parallel) {
    if (gzip == TRUE) {
        con <- gzfile(file, "rb")
    } else {
        con <- file(file, "rb")
    }
    suppressWarnings(skip_line <- readLines(con, n = 1))
    remove(skip_line)
    parse_chunck <- function(x, chr_name, col_names, col_types) {
        x <- readr::read_tsv(file = paste(mstrsplit(x), collapse = "\n"),
            col_types = col_types, skip = 0, n_max = Inf,
            col_names = col_names, progress = FALSE)
        x[x$chromosome == chr_name, ]
    }
    res <- chunk.apply(input = con, FUN = parse_chunck, chr_name = chr_name,
        col_names = col_names, col_types = col_types, CH.MAX.SIZE = buffer,
        parallel = parallel)
    close(con)
    res
}

read.seqz.tbi <- function(file, chr_name, col_names, col_types) {
    res <- seqminer::tabix.read(file, chr_name)
    res <- readr::read_tsv(file = paste(mstrsplit(res), collapse = "\n"),
        col_types = col_types, skip = 0, n_max = Inf,
        col_names = col_names, progress = FALSE)
}
