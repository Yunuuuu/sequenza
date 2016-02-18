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
   calc.dr <- function(x) {
      theoretical.CNt(x, CNn = 2, cellularity, ploidy, normal.ploidy = normal.ploidy,
                      avg.depth.ratio = avg.depth.ratio)
   }
   ci <- qt(ci, df = N - 1) * sd / sqrt(N)
   dr.left  <- calc.dr(depth.ratio - ci)
   dr.right <- calc.dr(depth.ratio + ci)
   dr       <- calc.dr(depth.ratio)
   cbind(CNt.float = dr, CCF.ratio.left = dr.left/CNt, CCF.ratio = dr/CNt, CCF.ratio.right = dr.right/CNt)
}

baf.ccf <- function(cellularity, Bf, B, CNt, CNn = 2, sd, N, ci = 0.95) {
   eBf <- expected.baf(CNt = CNt, cellularity = cellularity, sd = sd)
   eBf <- eBf$BAF[eBf$B == B & eBf$CNt == CNt]
   calc.ccf <- function(x) {
      rel.freq.ebf <- eBf / cellularity
      rel.freq.bf  <- x / cellularity
      tumor.comp  <- cellularity * B
      normal.comp <- 1 - cellularity
      eb.multiplicity <- rel.freq.ebf * (tumor.comp + normal.comp)
      b.multiplicity  <- rel.freq.bf * (tumor.comp + normal.comp)
      b.multiplicity / eb.multiplicity
   }
   ci <- qt(ci, df = N - 1) * sd / sqrt(N)
   ccf.left  <- calc.ccf(Bf - ci)
   ccf.right <- calc.ccf(Bf + ci)
   ccf       <- calc.ccf(Bf)
   cbind(CCF.baf.left = ccf.left, CCF.baf = ccf, CCF.baf.right = ccf.right)
}

sequenza.subclonal <- function(sequenza.extract, cp.table = NULL, sample.id, out.dir = getwd(),
                               cellularity = NULL, ploidy = NULL, female = TRUE, CNt.max = 20,
                               ratio.priority = FALSE, XY = c(X = "X", Y = "Y"), chromosome.list = 1:24){
   if(!file.exists(out.dir)) {
      dir.ok <- dir.create(path = out.dir, recursive = TRUE)
      if(!dir.ok) stop('Directory does not exist and cannot be created: ', out.dir)
   }
   makeFilename <- function(x) file.path(out.dir, paste(sample.id, x, sep = '_'))

   muts.file   <- makeFilename("mutations_ccf.txt")
   segs.file   <- makeFilename("segments_floats.txt")
   segs.d.plot <- makeFilename("dirichlet_segment.pdf")

   seg.tab     <- do.call(rbind, sequenza.extract$segments[chromosome.list])
   seg.tab     <- na.exclude(seg.tab)
   seg.len     <- (seg.tab$end.pos - seg.tab$start.pos)/1e6
   #avg.depth.ratio <- mean(sequenza.extract$gc$adj[, 2])
   #avg.depth.ratio <- weighted.mean(x = seg.tab$depth.ratio, w = seg.len)
   #avg.depth.ratio <- center.ratio(seg.tab)
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
   ccf.dr     <- lapply(split(seg.res,seq(NROW(seg.res))), function(x) get.ccf.dr(x, CNn = 2))
   ccf.baf    <- lapply(split(seg.res,seq(NROW(seg.res))), function(x) get.ccf.baf(x, CNn = 2))
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
         ccf.dr     <- lapply(split(seg.xy, seq(NROW(seg.xy))), function(x) get.ccf.dr(x, CNn = 1))
         ccf.baf    <- lapply(split(seg.xy, seq(NROW(seg.xy))), function(x) get.ccf.baf(x, CNn = 1))
         seg.xy     <- cbind(seg.xy, do.call(rbind, ccf.dr), do.call(rbind, ccf.baf))
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
   nsave <- 5000
   nskip <- 3
   ndisplay <- 1000
   mcmc <- list(nburn = nburn,
                nsave = nsave,
                nskip = nskip,
                ndisplay = ndisplay)
   fit1 <- DPMdencens(left = cbind(seg.res$CCF.ratio.left[seg.res$CNt > 0], seg.res$CCF.baf.left[seg.res$CNt > 0]),
                      right = cbind(seg.res$CCF.ratio.right[seg.res$CNt > 0], seg.res$CCF.baf.right[seg.res$CNt > 0]),
                      ngrid = 100, prior = prior, mcmc = mcmc,
                      state = state, status = TRUE)

   dir.plot <- function (dp, colFn = colorRampPalette(c('white', 'red')), ci = 0.95, ...) {
      z <- matrix(rank(dp$fbiv[[1]]), nrow = nrow(dp$fbiv[[1]])) / length(dp$fbiv[[1]])
      map <- makecmap(c(ci, 1), colFn = colFn, include.lowest = FALSE)
      colorgram(x = dp$grid[, 1], y = dp$grid[, 2], z = z,
                map = map, outlier="white", key = NA, ...)
   }
   pdf(segs.d.plot, width = 4, height = 4)
      dir.plot(dp = fit1, n = 200, ci = 0.99, las = 1, xlab = "CCF ratio", ylab = "CCF Bf", xlim = c(0,2), ylim = c(0,2))
      contour(x = fit1$grid[, 1], y = fit1$grid[, 2], z = fit1$fbiv[[1]] ,
              levels = quantile(fit1$fbiv[[1]], c(.999)), drawlabels = FALSE, add = T,
              method = "edge", lty = 1, lwd = 1)
      abline(h = c(0,0.5,1), v = c(0,0.5,1), lty = 2, lwd = 0.8)
   dev.off()
}


