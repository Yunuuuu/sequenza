theoretical_depth_ratio <- function(CNt, CNn = 2, cellularity, ploidy,
    normal_ploidy = 2, avg_depth_ratio = 1) {
    cellu_copy_term <- (1 - cellularity) + (CNt / CNn * cellularity)
    ploidy_cellu_term <- (ploidy / normal_ploidy * cellularity)
        + 1 - cellularity
    avg_depth_ratio * cellu_copy_term / ploidy_cellu_term
}

theoretical_baf <- function(cellularity, CNt, B, CNn = 2) {
    baf <- ( (B * cellularity) + ( 1 - cellularity) ) /
        ( (CNt * cellularity) + CNn * ( 1 - cellularity) )
    baf[CNn <= 1] <- NA
    baf
}

theoretical_mufreq <- function(Mt, CNt, CNn = 2, cellularity) {
    normal_alleles <- (CNt - Mt) * cellularity + CNn * (1 - cellularity)
    all_alleles <- (CNt * cellularity) + CNn * (1 - cellularity)
    1 - (normal_alleles / all_alleles)
}

baf_types_matrix <- function(CNt_min, CNt_max, CNn = 2) {
    cn_ratio_vect <- seq(from = CNt_min / CNn, to = CNt_max / CNn,
        by = 1 / CNn)
    CNt <- cn_ratio_vect * CNn
    if (CNn < 2) {
        b_comb <- lapply(CNt, FUN = function(x) 0)
    } else {
        b_comb <- lapply(CNt, FUN = function(x) {
            seq(from = 0, to = trunc(x / 2))
        })
    }
    times_b <- sapply(b_comb, length)
    CNt <- rep(CNt, times = times_b)
    B <- unlist(b_comb)
    cbind(CNn = CNn, CNt = CNt, B = B, A = CNt - B)
}

mufreq_types_matrix <- function(CNt_min, CNt_max, CNn = 2) {
    cn_ratio_vect <- seq(from = CNt_min / CNn,
        to = CNt_max / CNn, by = 1 / CNn)
    CNt <- cn_ratio_vect * CNn
    mut_comb <- lapply(CNt, FUN = function(x) seq(from = 0, to = x))
    times_muts <- sapply(mut_comb, length)
    cbind(CNn = CNn, CNt = rep(CNt, times = times_muts), Mt = unlist(mut_comb) )
}

baf_model_points <- function(cellularity, ploidy, baf_types, avg_depth_ratio) {
    depth_ratio <- theoretical_depth_ratio(cellularity = cellularity,
        ploidy = ploidy, CNn = baf_types[, "CNn"], CNt = baf_types[, "CNt"],
        avg_depth_ratio = avg_depth_ratio)
    baf <- theoretical_baf(cellularity = cellularity, CNn = baf_types[, "CNn"],
        CNt = baf_types[, "CNt"], B = baf_types[, "B"])
    cbind(BAF = baf, depth_ratio = depth_ratio)
}

mufreq_model_points <- function(cellularity, ploidy, mufreq_types,
    avg_depth_ratio) {
    mufreqs <- theoretical_mufreq(cellularity = cellularity,
        CNn = mufreq_types[, "CNn"], CNt = mufreq_types[, "CNt"],
        Mt = mufreq_types[, "Mt"])
    depth_ratio <- theoretical_depth_ratio(cellularity = cellularity,
        ploidy = ploidy, CNn = mufreq_types[, "CNn"],
        CNt = mufreq_types[, "CNt"], avg_depth_ratio = avg_depth_ratio)
    cbind(mufreqs, depth_ratio)
}
