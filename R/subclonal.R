baf.test <- function(Bf, cellularity, B, CNt, CNn = 2, sd.Bf = 0.1) {
   if (!(is.na(CNt) | is.na(B))) {
         if (B >= CNt) {
         B = CNt - B
      }
   } else {
      if (is.na(CNt)) {
         CNt = 1
      }
      if (is.na(B)) {
         B = 1
      }
   }
   model.baf  <- expected.baf(sd = sd.Bf, CNn = CNn, CNt = CNt, cellularity = cellularity)
   test.baf   <- model.baf[model.baf$B == B & model.baf$CNt == CNt, ]$BAF
   #min.offset <- 1e-323
   dt2(mean = Bf, sd = sd.Bf, x = test.baf, df = 5, log = TRUE)
}

baf.test.fit <- function(cellularity = seq(0.3, 1, by = 0.01), ...) {
   fit.cell <- function(cellularity) {
      baf.test(cellularity = cellularity, ...)
   }
   res <- sapply(cellularity, fit.cell)
   cbind(cellularity = cellularity, LPP = res)
}

ratio.test <- function(depth.ratio, cellularity, ploidy, CNt, CNn = 2, sd.ratio = 0.1, avg.depth.ratio = 1) {
   test.ratio  <- theoretical.depth.ratio(CNt = CNt, CNn = CNn, cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio)
   dt2(mean = depth.ratio, sd = sd.ratio, x = test.ratio, df = 5, log = TRUE)
}

ratio.test.fit <- function(cellularity = seq(0.3, 1, by = 0.01), ...) {
   fit.cell <- function(cellularity) {
      ratio.test(cellularity = cellularity, ...)
   }
   res <- sapply(cellularity, fit.cell)
   cbind(cellularity = cellularity, LPP = res)
}

mufreq.ccf <- function(cellularity, mufreq, Mt, CNt, CNn = 2, N, ci = 0.95) {
   calc.ccf <- function (x) {
      rel.freq    <- x/cellularity
      tumor.comp  <- cellularity * CNt
      normal.comp <- CNn * (1 - cellularity)
      mutation.multiplicity <- rel.freq * (tumor.comp + normal.comp)
      mutation.multiplicity/Mt

   }
   mufreq <- binom.test(x = round(mufreq * N, 0), n = N, conf.level = 0.95)
   ccf.left  <- calc.ccf(mufreq$conf.int[1])
   ccf.right <- calc.ccf(mufreq$conf.int[2])
   ccf       <- calc.ccf(mufreq$estimate)
   cbind(CCF.left = ccf.left, CCF = ccf, CCF.right = ccf.right)
}

ratio.ccf <- function(depth.ratio, CNt, CNn = 2, cellularity, ploidy,
                      normal.ploidy = 2, avg.depth.ratio = 1, sd, N, ci = 0.95) {
   dr <- theoretical.depth.ratio(CNt = 0:20, cellularity = cellularity, ploidy = ploidy,
                                 normal.ploidy = normal.ploidy, avg.depth.ratio = avg.depth.ratio,
                                 CNn = CNn)
   ratios <- data.frame(CNt=0:20, ratio = dr)
   norm_cn <- which.min(abs(ratios$ratio - 1))
   mid_cn = ratios$CNt[norm_cn]
   mid_dr = ratios$ratio[norm_cn]
   if (!is.na(CNt)) {
      if (CNt == mid_cn) {
         if(depth.ratio > mid_dr) {
            CNt = CNt + 1
         } else {
            CNt = CNt - 1
         }
      }
   }
   res <- ratio.test.fit (depth.ratio = depth.ratio, CNt = CNt, CNn = CNn,
                           ploidy = ploidy, avg.depth.ratio = avg.depth.ratio,
                           cellularity = seq(from = 0, to = 1, by = 0.01))
   calc.dr <- function(x) {
      theoretical.CNt(x, CNn = CNn, cellularity, ploidy, normal.ploidy = normal.ploidy,
                      avg.depth.ratio = avg.depth.ratio)
   }
   dr.f       <- calc.dr(depth.ratio)
   range.lpp  <- quantile(res[, 2], ci, na.rm = TRUE)
   if (is.na( range.lpp)) {
      cbind(CNt.float = dr.f, CCF.ratio.left = 1, CCF.ratio =  1, CCF.ratio.right = 1)
   } else {
      res        <- res[res[, 2] >= range.lpp, ]
      ccf.left   <- res[which.min(res[, 1]), 1] / cellularity
      ccf        <- res[which.max(res[, 2]), 1] / cellularity
      ccf.right  <- res[which.max(res[, 1]), 1] / cellularity
      cbind(CNt.float = dr.f, CCF.ratio.left = ccf.left, CCF.ratio =  ccf, CCF.ratio.right = ccf.right)
   }
}

baf.ccf <- function(cellularity, Bf, B, CNt, CNn = 2, sd, N, ci = 0.95, offset = 1e-10) {
   if (!(is.na(CNt) | is.na(B))) {
      if ((CNt - B) == B) {
         if (B > 0) {
            B = B - 1
         } else {
            CNt = CNt + 1
         }
      }
      res <- baf.test.fit(Bf = Bf, CNt = CNt, CNn = CNn, sd.Bf = sd, B = B,
                       cellularity = seq(from = 0, to = 1, by = 0.01))
      range.lpp  <- quantile(res[, 2], ci, na.rm = TRUE)
   } else {
      range.lpp <- NA
   }
   if (is.na(range.lpp)) {
      cbind(CCF.baf.left = 1, CCF.baf = 1, CCF.baf.right = 1)
   } else {
      res        <- res[res[, 2] >= range.lpp, ]
      ccf.left   <- res[which.min(res[, 1]), 1] / cellularity
      ccf        <- res[which.max(res[, 2]), 1] / cellularity
      ccf.right  <- res[which.max(res[, 1]), 1] / cellularity
      cbind(CCF.baf.left = ccf.left, CCF.baf = ccf, CCF.baf.right = ccf.right)
   }
}

get_clust_info <- function(segs, min.size = 1e3) {

   seg.size <- segs$end.pos - segs$start.pos
   segs     <- segs[seg.size >= min.size, ]
   seg.size <- seg.size[seg.size >= min.size]
   tot.size <- sum(as.numeric(seg.size))
   split_cl <- split(x = data.frame(baf = segs$CCF.baf, ratio = segs$CCF.ratio, size = seg.size),
                     f = segs$cluster)
   size     <- lapply(split(x = seg.size/tot.size, f = segs$cluster),
                      FUN = function(x) sum(x, na.rm = TRUE))
   ccf.b    <-  lapply(split_cl,
                       FUN = function(x) weighted.mean(x$baf, w = sqrt(x$size), na.rm = TRUE))
   ccf.r    <-  lapply(split_cl,
                       FUN = function(x) weighted.mean(x$ratio, w = sqrt(x$size), na.rm = TRUE))
   list(size = size, CCF.baf = ccf.b, CCF.ratio = ccf.r)
}

chromosomes.info <- function() {
   #remove chr9 weird area chr9:47367679-50367679 -> chr9:37890000-71410000
   data.frame(
      chrom = c(1:22, "X", "Y"),
      start.cent = c(121535434, 92326171, 90504854, 49660117, 46405641,
                     58830166, 58054331, 43838887, 37890000, 39254935,
                     51644205, 34856694, 16000000, 16000000, 17000000,
                     35335801, 22263006, 15460898, 24681782, 26369569,
                     11288129, 13000000, 58632012, 10104553),
      end.cent = c(124535434, 95326171, 93504854, 52660117, 49405641,
                   61830166, 61054331, 46838887, 71410000, 42254935,
                   54644205, 37856694, 19000000, 19000000, 20000000,
                   38335801, 25263006, 18460898, 27681782, 29369569,
                   14288129, 16000000, 61632012, 13104553),
      end.pos = c(249250621, 243199373, 198022430, 191154276, 180915260,
                  171115067, 159138663, 146364022, 141213431, 135534747,
                  135006516, 133851895, 115169878, 107349540, 102531392,
                  90354753, 81195210, 78077248, 63025520, 59128983, 51304566,
                  48129895, 155270560, 59373566)
   )
}

compare_segs_cluster <- function(seg_1, seg_2, c_info, max_gap = 1e3) {
   both_before_cent <- seg_1$end.pos <  c_info$start.cent & seg_1$end.pos < c_info$start.cent
   both_after_cent <- seg_1$start.pos > c_info$end.cent
   if (both_before_cent | both_after_cent) {
      if (seg_1$end.pos <= (seg_2$start.pos + max_gap) ) {
         if (seg_1$CNt == seg_2$CNt & seg_1$A == seg_2$A & seg_1$cluster == seg_2$cluster) {
            seg_tot = seg_1
            seg_tot$end.pos = seg_2$end.pos
            return(seg_tot)
         } else {
            return(seg_1)
         }
      } else {
         return(seg_1)
      }
   } else {
      return(seg_1)
   }
}

merge_segs_cluster <- function(segs, max_gap = 1e3) {
   chromosomes_info <- chromosomes.info()
   if (length(grep("chr", segs$chromosome)) > 0) {
      chromosomes_info$chrom <- paste("chr", chromosomes_info$chrom, sep = "")
   }
   chromosomes <- unique(segs$chromosome)
   segs.split  <- split(segs, segs$chromosome)
   chr_list <- list()
   for (chr in chromosomes) {
      c_info   <- chromosomes_info[chromosomes_info$chrom == chr, ]
      segs_chr <- segs.split[[chr]]
      merged   <- list()
      segs_t   <- NULL
      if (nrow(segs_chr) == 1) {
         i = 0
      } else {
         for (i in 1:(nrow(segs_chr) - 1)) {
            if (is.null(segs_t)) {
               segs_t <- compare_segs_cluster(segs_chr[i, ], segs_chr[i + 1, ], c_info)
            } else {
               segs_t <- compare_segs_cluster(segs_t, segs_chr[i, ], c_info)
            }
            if (segs_t$end.pos == segs_chr[i, 'end.pos'] &
                segs_t$start.pos == segs_chr[i, 'start.pos']) {
               merged <- c(merged, list(segs_t))
               segs_t <- NULL
            } else if(segs_t$end.pos == segs_chr[i, 'end.pos'] &
                      segs_t$start.pos != segs_chr[i, 'start.pos']) {
               merged <- c(merged, list(segs_t))
            } else if (segs_t$end.pos != segs_chr[i, 'end.pos'] &
                       segs_t$start.pos != segs_chr[i, 'start.pos'])
               merged <- c(merged, list(segs_t))
            segs_t <- NULL
         }
      }
      if (!is.null(segs_t)) {
         merged <- c(merged, list(segs_t))
      } else {
         merged <- c(merged, list(segs_chr[i + 1, ]))
      }
      chr_list <- c(chr_list, list(do.call(rbind, merged)))
   }
   return(do.call(rbind, chr_list))
}

sequenza.subclonal <- function(sequenza.extract, cp.table = NULL, sample.id, out.dir = getwd(),
                               cellularity = NULL, ploidy = NULL, female = TRUE, CNt.max = 20, subclonal_filter = 0.2,
                               ratio.priority = FALSE, XY = c(X = "X", Y = "Y"), chromosome.list = 1:24){
   if(!file.exists(out.dir)) {
      dir.ok <- dir.create(path = out.dir, recursive = TRUE)
      if(!dir.ok) stop('Directory does not exist and cannot be created: ', out.dir)
   }
   makeFilename <- function(x) file.path(out.dir, paste(sample.id, x, sep = '_'))

   muts.file   <- makeFilename("mutations_ccf.txt")
   segs.file   <- makeFilename("segments_floats.txt")
   segs.d.plot <- makeFilename("dirichlet_segment.pdf")
   clust.file   <- makeFilename("segments_clusters.txt")
   clust.plot   <- makeFilename("segments_clusters.pdf")
   clust.genome <- makeFilename("genome_clusters.pdf")

   seg.tab     <- do.call(rbind, sequenza.extract$segments[chromosome.list])
   #seg.tab     <- na.exclude(seg.tab)

   seg.len     <- (seg.tab$end.pos - seg.tab$start.pos)/1e6

   avg.depth.ratio <- 1

   if (is.null(cp.table) && (is.null(cellularity) || is.null(ploidy))){
      stop("Either the cp.table or both cellularity and ploidy argument are required.")
   }
   if (!is.null(cp.table)){
      cint <- get.ci(cp.table)
      if (!is.null(cellularity) || !is.null(ploidy)) {
         if (is.null(cellularity)) cellularity <- cint$max.cellularity
         if (is.null(ploidy)) ploidy <- cint$max.ploidy
      } else {
         cellularity <- cint$max.cellularity
         ploidy <- cint$max.ploidy
      }
   }
   mut.tab     <- na.exclude(do.call(rbind, sequenza.extract$mutations[chromosome.list]))
   get.ccf.dr <- function(x, CNn) {
      ratio.ccf(cellularity = cellularity,
                ploidy = ploidy,
                avg.depth.ratio = avg.depth.ratio,
                depth.ratio = as.numeric(x['depth.ratio']),
                CNt = as.numeric(x['CNt']),
                sd = as.numeric(x['sd.ratio']),
                N = as.numeric(x['N.ratio']),
                CNn = CNn)
   }
   get.ccf.baf <- function(x, CNn) {
      baf.ccf(cellularity = cellularity,
              Bf = as.numeric(x['Bf']),
              CNt = as.numeric(x['CNt']),
              B = as.numeric(x['B']),
              sd = as.numeric(x['sd.BAF']),
              N = as.numeric(x['N.BAF']),
              CNn = CNn)
   }
   secondary_cn_subclonal <- function(segments, cellularity, ploidy, avg.depth.ratio) {
      subclonal <- segments$status == 'subclonal'
      types  <- types.matrix(CNt.min=0, CNt.max=max(segments$CNt, na.rm = TRUE), CNn=2)
      points <- model.points(cellularity, ploidy, avg.depth.ratio=avg.depth.ratio, types = types)
      points <- cbind(types, points)
      points <- unique(cbind(points$CNt, points$depth.ratio))
      #mid_cnt <-  points[which.min(abs(points[, 2] - avg.depth.ratio)), 1]
      points <- setNames(nm= points[, 1], points[, 2])
      secondary_cnt <- segments$CNt
      more_cnt <- points[as.character(segments$CNt)] > segments$depth.ratio & subclonal & segments$CNt > 0 # & segments$CNt != (mid_cnt - 1)
      less_cnt <- points[as.character(segments$CNt)] < segments$depth.ratio & subclonal # & segments$CNt != (mid_cnt + 1)
      secondary_cnt[more_cnt] <- secondary_cnt[more_cnt] + 1
      secondary_cnt[less_cnt] <- secondary_cnt[less_cnt] - 1
      cbind(segments, subCNt = secondary_cnt)
   }
   if (female){
      segs.is.xy <- seg.tab$chromosome == XY["Y"]
      mut.is.xy  <- mut.tab$chromosome == XY["Y"]
   } else{
      segs.is.xy <- seg.tab$chromosome %in% XY
      mut.is.xy  <- mut.tab$chromosome %in% XY
      types.xy = types.matrix(CNt.min = 0, CNt.max = CNt.max, CNn = 1)
   }

   avg.sd.ratio  <- sum(seg.tab$sd.ratio * seg.tab$N.ratio, na.rm = TRUE)/sum(seg.tab$N.ratio, na.rm = TRUE)
   avg.sd.Bf     <- sum(seg.tab$sd.BAF * seg.tab$N.BAF, na.rm = TRUE)/sum(seg.tab$N.BAF, na.rm = TRUE)
   cn.alleles  <- baf.bayes(Bf = seg.tab$Bf[!segs.is.xy], CNt.max = CNt.max,
                            depth.ratio = seg.tab$depth.ratio[!segs.is.xy],
                            cellularity = cellularity, ploidy = ploidy,
                            avg.depth.ratio = avg.depth.ratio, sd.ratio = seg.tab$sd.ratio[!segs.is.xy],
                            weight.ratio = seg.len[!segs.is.xy], sd.Bf = seg.tab$sd.BAF[!segs.is.xy],
                            weight.Bf = 1, ratio.priority = ratio.priority, CNn = 2)

   seg.res    <- cbind(seg.tab[!segs.is.xy, ], cn.alleles)

   #no_ai      <- which(seg.res$A == seg.res$B & seg.res$CNt != 0)
   #dr <- theoretical.depth.ratio(CNt=0:CNt.max, cellularity = cellularity,
   #                              ploidy = ploidy)
   #ratios <- data.frame(CNt=0:CNt.max, ratio = dr)
   #norm_cn <- which.min(abs(ratios$ratio - 1))
   #mid_cn = ratios$CNt[norm_cn]
   #no_mean    <- which(seg.res$CNt == mid_cn)
   #if (mid_cn%%2 == 0) {
   #   no_cnv     <- sort(unique(c(no_ai, no_mean)))
   #} else {
   #   no_cnv     <- no_ai
   #}
   seg.res_tmp <- seg.res
   fix_CNt <- function(segs, mid_cn, avg.depth.ratio) {
      segs_gain <- which(segs$depth.ratio >= avg.depth.ratio & segs$CNt > 0)
      segs_loss <- which(segs$depth.ratio <= avg.depth.ratio & segs$CNt > 0)
      segs_mid <- which(segs$CNt == mid_cn)
      mid_gain <- sort(unique(intersect(segs_gain, segs_mid)))
      mid_loss <- sort(unique(intersect(segs_loss, segs_mid)))
      if (mid_cn%%2 == 0) {
         segs$CNt[mid_gain] <- segs$CNt[mid_gain] + 1
         segs$A[mid_gain] <- segs$A[mid_gain] + 1
         segs$CNt[mid_loss] <- segs$CNt[mid_loss] - 1
         segs$A[mid_gain] <- segs$A[mid_gain] - 1
      }

      zero_b <- which(segs$B == 0)
      more_b <- which(segs$B > 0)
      segs$B[sort(unique(intersect(segs_loss, more_b)))] <- segs$B[sort(unique(intersect(segs_loss, more_b)))] - 1
      segs$A[sort(unique(intersect(segs_loss, zero_b)))] <- segs$A[sort(unique(intersect(segs_loss, zero_b)))] - 1
      segs$A[segs_zero] <- segs$A[segs_zero] + 1
      segs$CNt[segs_zero] <- segs$CNt[segs_zero] + 1
      segs
   }
   #fixed_nocv <-  fix_CNt(seg.res_tmp[no_cnv, ], avg.depth.ratio = avg.depth.ratio)
   #seg.res_tmp$CNt[no_cnv] <- fixed_nocv$CNt
   #seg.res_tmp$A[no_cnv] <- fixed_nocv$A
   #seg.res_tmp$B[no_cnv] <- fixed_nocv$B
   seg.res_tmp$CNt[which(seg.res_tmp$CNt == 0)] <- 1
   ccf.dr     <- lapply(split(seg.res_tmp, seq(nrow(seg.res))), function(x) get.ccf.dr(x, CNn = 2))
   cat("done1\n")
   ccf.baf    <- lapply(split(seg.res_tmp, seq(nrow(seg.res))), function(x) get.ccf.baf(x, CNn = 2))
   cat("done2\n")
   cat(dim(seg.res), "segs\n")
   cat(dim(do.call(rbind, ccf.baf)), "baf\n")
   cat(dim(do.call(rbind, ccf.dr)), "dr\n")
   seg.res    <- cbind(seg.res, do.call(rbind, ccf.dr), do.call(rbind, ccf.baf))
   if (!female){
      if (sum(segs.is.xy) >= 1) {
         cn.alleles  <- baf.bayes(Bf = NA, CNt.max = CNt.max,
                                  depth.ratio = seg.tab$depth.ratio[segs.is.xy],
                                  cellularity = cellularity, ploidy = ploidy,
                                  avg.depth.ratio = avg.depth.ratio, sd.ratio = seg.tab$sd.ratio[segs.is.xy],
                                  weight.ratio = seg.len[segs.is.xy], sd.Bf = NA,
                                  weight.Bf = NA, ratio.priority = ratio.priority, CNn = 1)
         seg.xy     <- cbind(seg.tab[segs.is.xy, ], cn.alleles)
         ccf.dr     <- lapply(split(seg.xy, seq(nrow(seg.xy))), function(x) get.ccf.dr(x, CNn = 1))
         #ccf.baf    <- lapply(split(seg.xy, seq(nrow(seg.xy))), function(x) get.ccf.baf(x, CNn = 1))
         # set the CCF baf to 1 for homozygous genomes
         seg.xy     <- cbind(seg.xy, do.call(rbind, ccf.dr),  CCF.baf.left = 1, CCF.baf = 1, CCF.baf.right = 1)
         seg.res    <- rbind(seg.res, seg.xy)
      }
   }

   write.table(seg.res, file = segs.file,
               col.names = TRUE, row.names = FALSE, sep = "\t")
   if(nrow(mut.tab) > 0) {
      mut.alleles  <- mufreq.bayes(mufreq = mut.tab$F[!mut.is.xy], CNt.max = CNt.max,
                                   depth.ratio = mut.tab$adjusted.ratio[!mut.is.xy],
                                   cellularity = cellularity, ploidy = ploidy,
                                   avg.depth.ratio = avg.depth.ratio, CNn = 2)
      mut.res     <- cbind(mut.tab[!mut.is.xy, ], mut.alleles)
      if (!female){
         if (sum(mut.is.xy) >= 1) {
            mut.alleles  <- mufreq.bayes(mufreq = mut.tab$F[mut.is.xy], CNt.max = CNt.max,
                                         depth.ratio = mut.tab$adjusted.ratio[mut.is.xy],
                                         cellularity = cellularity, ploidy = ploidy,
                                         avg.depth.ratio = avg.depth.ratio, CNn = 1)
            mut.xy     <- cbind(mut.tab[mut.is.xy, ], mut.alleles)
            mut.res    <- rbind(mut.res, mut.xy)
         }
      }

      get.ccf.mut <- function(x) {
         Mta <- as.numeric(x['Mt'])
         if (Mta == 0) {Mta = 1}
         mufreq.ccf(cellularity = cellularity,
                    mufreq = as.numeric(x['F']),
                    CNt = as.numeric(x['CNt']),
                    Mt = Mta,
                    N = as.numeric(x['good.reads']),
                    CNn = as.numeric(x['CNn']))
      }
      ccfs <-  lapply(split(mut.res, seq(NROW(mut.res))), function(x) get.ccf.mut(x))
      mut.res <- cbind(mut.res, do.call(rbind, ccfs))
      write.table(mut.res, file = muts.file,
                  col.names = TRUE, row.names = FALSE, sep = "\t")
   }

   s2 <- matrix(c(10000, 0, 0, 1), ncol = 2)
   m2 <- c(180, 3)
   psiinv2 <- diag(c(1/10000, 1), 2)
   prior <- list(alpha = 1, nu1 = 4,
                 nu2 = 4, s2 = s2,
                 m2 = m2, psiinv2 = psiinv2,
                 tau1 = 0.01, tau2 = 0.01)
   nburn <- 5000
   nsave <- 8000
   nskip <- 3
   ndisplay <- 1000
   mcmc <- list(nburn = nburn,
                nsave = nsave,
                nskip = nskip,
                ndisplay = ndisplay)
   fit1 <- DPMdencens(left = cbind(seg.res$CCF.ratio.left, seg.res$CCF.baf.left),
                      right = cbind(seg.res$CCF.ratio.right, seg.res$CCF.baf.right),
                      ngrid = 100, prior = prior, mcmc = mcmc,
                      state = state, status = TRUE)

   dir.plot <- function (dp, colFn = colorRampPalette(c('white', 'red')), ci = 0.95, ...) {
      z <- matrix(rank(dp$fbiv[[1]]), nrow = nrow(dp$fbiv[[1]])) / length(dp$fbiv[[1]])
      map <- makecmap(c(ci, 1), colFn = colFn, include.lowest = FALSE)
      suppressWarnings(colorgram(x = dp$grid[, 1], y = dp$grid[, 2], z = z,
                map = map, outlier="white", key = NA, ...))
   }

   #seg.res  <- seg.res[seg.res$CNt > 0, ]
   seg.size <- (seg.res$end.pos - seg.res$start.pos) /1e6

   #seg.res.merged <- merge_segs_cluster(cbind(seg.res, cluster = fit1$state$ss))
   seg.res <- cbind(seg.res, cluster = fit1$state$ss)

   clust <- get_clust_info(seg.res)

   clust <- data.frame(cluster = names(clust$size),
                       size = do.call(c, clust$size),
                       CCF.ratio = do.call(c, clust$CCF.ratio),
                       CCF.baf = do.call(c, clust$CCF.baf),
                       stringsAsFactors = FALSE)
   #radius_clonal = 0.1
   #intercept_err = 0.75
   #clonal_1 = sqrt((1- clust$CCF.ratio)^2 + (1 - clust$CCF.baf)^2) <= radius_clonal
   #clonal_2 = sign((1 - intercept_err) * (clust$CCF.baf - 0) - (1 - 0) * (clust$CCF.ratio - intercept_err)) <= 0
   #clonal_3 = sign((1 - 0) * (clust$CCF.baf - intercept_err) - (1 - intercept_err) * (clust$CCF.ratio - 0)) >= 0


   #clonal = clonal_1 | clonal_2 | clonal_3
   #position = sign((Bx - Ax) * (Y - Ay) - (By - Ay) * (X - Ax))

   plot_rectangle <- function(intercept, ...) {
      segments(x0 = 0, y0 = intercept, x1 = 1 - intercept, y1 = 1, ...)
      segments(x0 = intercept, y0 = 0, x1 = 1, y1 = 1 - intercept, ...)
      segments(x0 = 0, y0 = intercept, x1 = intercept, y1 = 0, ...)
      segments(x0 = 1 - intercept, y0 = 1, x1 = 1, y1 = 1 - intercept, ...)

      segments(x0 = 0, y0 = 0, x1 = intercept, y1 = 1, ...)
      segments(x0 = 1, y0 = 1, x1 = 0, y1 = 1 - intercept, ...)
      segments(x0 = 0, y0 = 0, x1 = 1, y1 = intercept, ...)
      segments(x0 = 1, y0 = 1, x1 = 1 - intercept, y1 = 0, ...)
   }
   #clonal_1 = sign((1 - intercept_err) * (clust$CCF.baf - 0) - (1 - 0) * (clust$CCF.ratio - intercept_err)) <= 0
   intercept = subclonal_filter

   is_clonal <- function(x, y, intercept) {
      give_side <- function(x, y, a, b) {
         sign((b[1] - a[1]) * (y - a[2]) - (b[2] - a[2]) * (x - a[1]))
      }
      clonal_1 = give_side(x = x, y = y,
                           a = c(0, intercept), b =  c(1 - intercept, 1)) < 0
      clonal_2 = give_side(x = x, y = y,
                           a = c(intercept, 0), b =  c(1, 1 - intercept)) >0
      clonal_3 = give_side(x = x, y = y,
                           a = c(0, intercept), b =  c(intercept, 0)) > 0
      clonal_4 = give_side(x = x, y = y,
                           a = c(1 - intercept, 1), b =  c(1, 1 - intercept)) < 0
      clonal_5 = give_side(x = x, y = y,
                           a = c(0, 0), b =  c(intercept, 1)) < 0
      clonal_6 = give_side(x = x, y = y,
                           a = c(1, 1), b =  c(0, 1 - intercept)) > 0
      clonal_7 = give_side(x = x, y = y,
                           a = c(0, 0), b =  c(1, intercept)) > 0
      clonal_8 = give_side(x = x, y = y,
                           a = c(1, 1), b =  c(1 - intercept, 0)) < 0

      subclonal = clonal_1 & clonal_2 & clonal_3 & clonal_4 & clonal_5 & clonal_6 & clonal_7 & clonal_8
      !subclonal
   }

   clonal = is_clonal(x = clust$CCF.ratio, y = clust$CCF.baf, intercept = intercept)

   clust_status = setNames(c("subclonal", "clonal")[clonal + 1], clust$cluster)
   seg_status = clust_status[as.character(seg.res$cluster)]
   which.subclonal <- which(seg_status == "subclonal")
   seg_status[which.subclonal] <-  c("subclonal", "clonal")[
      is_clonal(x = seg.res$CCF.ratio.left[which.subclonal],
                y = seg.res$CCF.baf.left[which.subclonal],
                intercept = intercept) + 1 ]
   which.subclonal <- which(seg_status == "subclonal")
   seg_status[which.subclonal] <-  c("subclonal", "clonal")[
      is_clonal(x = seg.res$CCF.ratio.right[which.subclonal],
                y = seg.res$CCF.baf.right[which.subclonal],
                intercept = intercept) + 1 ]
   seg.res <- cbind(seg.res, status = seg_status)
   seg.res <- secondary_cn_subclonal(seg.res, cellularity, ploidy, avg.depth.ratio)
   write.table(seg.res, file = clust.file,
               col.names = TRUE, row.names = FALSE, sep = "\t")

   pdf(segs.d.plot, width = 6, height = 6)
      dir.plot(dp = fit1, n = 200, ci = 0.95, las = 1, xlab = "CCF Copy number", ylab = "CCF B allele",
               xlim = c(0,2), ylim = c(0,2), colFn = colorRampPalette(c('yellow', "orange", "red")))
      #contour(x = fit1$grid[, 1], y = fit1$grid[, 2], z = fit1$fbiv[[1]] ,
      #        levels = quantile(fit1$fbiv[[1]], c(.95)), drawlabels = FALSE, add = T,
      #        method = "edge", lty = 1, lwd = 1)
      abline(h = c(0, 0.5, 1, 1.5),
             v = c(0, 0.5, 1, 1.5),
             lty = 2, lwd = 0.8)
      points(x = clust$CCF.ratio, y = clust$CCF.baf,
             cex = log(clust$size * 100), pch = 21)
      text(clust$cluster, x = clust$CCF.ratio, y = clust$CCF.baf)
   dev.off()

   pdf(clust.plot, width = 6, height = 6)
      plot(x = 0, y = 0, ylim = c(0, 1), las = 1, bty = 'n', type = 'n',
           xlim = c(0, 1), xlab = "CCF Copy number", ylab = "CCF B allele")
      abline(a = 0, b = 1, lty = 2)
      abline(h = seq(0, 1, 0.25), lty = 2, lwd = 0.5)
      abline(v = seq(0, 1, 0.25), lty = 2, lwd = 0.5)
      #lines(x=c(intercept_err, 1), y = c(0, 1), lty = 2, lwd = 0.2)
      #lines(x=c(0, 1), y = c(intercept_err, 1), lty = 2, lwd = 0.2)
      plot_rectangle(intercept = intercept, lty = 2, lwd = 0.2 )
      points(x = clust$CCF.ratio, y = clust$CCF.baf,
             cex = log(clust$size * 100), pch = 21, col = NA,
             bg = c("lightblue", "lightsalmon")[clonal + 1])
      text(clust$cluster, x = clust$CCF.ratio, y = clust$CCF.baf,
           pos = 2, offset = clust$size)
      legend("topleft", legend = c("clonal", "subclonal"),
             fill = c("lightsalmon","lightblue"),
             bty = "n", border = FALSE)
      plot(seg.res$CCF.ratio, seg.res$CCF.baf, col = fit1$state$ss)
      plot(seg.res$CCF.ratio[seg.size >=3], seg.res$CCF.baf[seg.size >=3], col = fit1$state$ss[seg.size >=3],
           xlim = c(0, 1), ylim = c(0, 1), xlab = "CCF depth ratio", ylab = "CCF B allele frequency")
      plot(seg.res$Bf, seg.res$depth.ratio, col = fit1$state$ss, xlim = c(0, 0.5), ylim = c(0, 2.5))
      baf.ratio.model.fit(cellularity = cellularity, ploidy = ploidy, segs = seg.res, col = fit1$state$ss)
   dev.off()

   pdf(clust.genome, width = 15, height = 5)
      genome.view(seg.res, info.type = "colors",
                  col = c("lightblue", "lightsalmon")[(seg.res$status == "clonal") + 1])
      legend("topleft", legend = c("clonal", "subclonal"),
             fill = c("lightsalmon","lightblue"),
             bty = "n", border = FALSE)
   dev.off()


}


