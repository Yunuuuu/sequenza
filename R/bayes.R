mufreq.dbinom <- function(mufreq, mufreq.model, depth.t, seq.errors = 0.01) {
   mufreq.model[mufreq.model == 0] <- seq.errors
   n.success       <- round(mufreq * depth.t, 0)
   dbinom( x = n.success, size = depth.t, prob = mufreq.model)
}

mufreq.dpois <- function(mufreq, mufreq.model, depth.t, seq.errors = 0.01) {
   mufreq.model[mufreq.model == 0] <- seq.errors
   n.success       <- round(mufreq * depth.t, 0)
   dpois( x = n.success, lambda = mufreq.model * depth.t)
}

baf.dbinom <- function(baf, baf.model, depth.t) {
   n.success       <- round(baf * depth.t, 0)
   dbinom( x = n.success, size = depth.t, prob = baf.model)
}

baf.dpois <- function(baf, baf.model, depth.t) {
   n.success       <- round(baf * depth.t, 0)
   dpois( x = n.success, lambda = baf.model * depth.t)
}

depth.ratio.dbinom <- function(size, depth.ratio, depth.ratio.model) {
   #n.success        <- round(depth.n * depth.ratio, 0)
   n.success        <- round(size * (depth.ratio/(1 + depth.ratio)), 0)
   prob             <- depth.ratio.model / (1 + depth.ratio.model)
   dbinom( x = n.success, size = size, prob = prob)
}

depth.ratio.dpois <- function(size, depth.ratio, depth.ratio.model) {
   #n.success        <- round(depth.n * depth.ratio, 0)
   n.success        <- round(size * (depth.ratio/(1 + depth.ratio)), 0)
   prob             <- depth.ratio.model / (1 + depth.ratio.model)
   dpois( x = n.success, lambda = prob * size)
}

mufreq.bayes <- function(mufreq, depth.ratio, cellularity, dna.content, avg.depth.ratio,
                         weight.mufreq = 10, weight.ratio = 10,
                         CNt.min = 1, CNt.max = 7, CNr = 2) {

   mufreq.tab <- data.frame(F = mufreq, ratio = depth.ratio,
                            weight.mufreq = weight.mufreq, weight.ratio = weight.ratio)
   types <- types.matrix(CNt.min = CNt.min, CNt.max = CNt.max, CNr = CNr)
   mufreq.depth.ratio <- cbind(types, model.points(cellularity = cellularity, dna.content = dna.content,
                                                   types = types, avg.depth.ratio = avg.depth.ratio))
   rows.x             <- 1:nrow(mufreq.tab)
   bayes.fit <- function (x, mat, model.pts) {
      test.ratio <- model.pts$depth.ratio
      test.mufrq <- model.pts$mufreqs
      min.offset <- 1e-323
      score.r    <- depth.ratio.dbinom(size = mat[x,]$weight.ratio, depth.ratio = mat[x,]$ratio, test.ratio)
      score.m    <- mufreq.dbinom(mufreq = mat[x,]$F, depth.t = mat[x,]$weight.mufreq, test.mufrq)
      
      priors <- rep(1, length(score.r))
      priors <- priors/sum(priors)

      score.r    <- score.r * priors
      score.m    <- score.m * priors

      post.model <- score.r * score.m

      post.model[post.model == 0] <- min.offset
      # max.lik <-  which.max(post.model)
      # max.post <- c(as.numeric(model.pts[max.lik,1:3]), log2(post.model[max.lik]))
      # max.post

      res.cn     <- model.pts$CNt[which.max(score.r)]
      idx.pts    <- model.pts$CNt == res.cn
      model.lik  <- cbind(model.pts[idx.pts, 1:3], log2(post.model[idx.pts]))
      if (is.null(dim(model.lik))) {
         max.post <- model.lik
      } else {
         max.post   <- model.lik[which.max(model.lik[,4]),]
      }

      max.post
   }
   types.L           <- mapply(FUN = bayes.fit, rows.x,
                         MoreArgs = list(mat = mufreq.tab, model.pts = mufreq.depth.ratio),
                                         SIMPLIFY = FALSE)
   types.L           <- do.call(rbind, types.L)
   colnames(types.L) <- c("CNr","CNt","Mt", "L")
   types.L
}

baf.bayes <- function(Bf, depth.ratio, cellularity, dna.content, avg.depth.ratio,
                      weight.Bf = 100, weight.ratio = 100, CNt.min = 0,
                      CNt.max = 7, CNr = 2) {
   
   mufreq.tab <- data.frame(Bf = Bf, ratio = depth.ratio,
                            weight.Bf = weight.Bf, weight.ratio = weight.ratio)
   mufreq.depth.ratio <- model.points(cellularity = cellularity, dna.content = dna.content, 
                                      types = cbind(CNr = CNr, CNt = CNt.min:CNt.max, Mt = 0),
                                      avg.depth.ratio = avg.depth.ratio)
   model.d.ratio      <- cbind(CNt = CNt.min:CNt.max, depth.ratio = mufreq.depth.ratio[, 2])
   model.baf          <- theoric.baf(CNr = CNr, CNt = CNt.max, cellularity = cellularity)
   if(CNt.min == 0) {
      model.baf          <- as.data.frame(rbind(c(0,0,0.5,0), model.baf))
   }
   # B-allele freq are never 0.5, always smaller. work around on this bias
   #model.baf$BAF[model.baf$BAF == 0.5] <- quantile(rep(mufreq.tab$Bf, times = mufreq.tab$weight.Bf, na.rm = T), probs = .75)
   model.pts          <- merge(model.baf, model.d.ratio)
   # model.pts          <- cbind(baf.type = apply(model.pts[, 1:3], 1, FUN = function(x) paste(x, collapse = "_")),
   #                             model.pts[, 4:5])
   rows.x             <- 1:nrow(mufreq.tab)

   bayes.fit <- function (x, mat, model.pts) {

      test.ratio <- model.pts$depth.ratio
      test.baf   <- model.pts$BAF
      min.offset <- 1e-323
      score.r    <- depth.ratio.dbinom(size = mat[x,]$weight.ratio, depth.ratio = mat[x,]$ratio, test.ratio)
      score.b    <- baf.dbinom(baf = mat[x,]$Bf, depth.t = mat[x,]$weight.Bf, test.baf)

      priors <- rep(1,length(score.b))
      priors.cn <- priors
      ## Some small priors on the diploid status
      priors.cn[model.pts$CNt == 2]  <- 3
      priors <- priors/sum(priors)
      priors.cn <- priors.cn/sum(priors.cn)

      score.r    <- score.r * priors.cn
      score.b    <- score.b * priors

      post.model <- score.r * score.b

      post.model[post.model == 0] <- min.offset
      max.lik <-  which.max(post.model)
      max.post <- c(as.numeric(model.pts[max.lik,1:3]), log2(post.model[max.lik]))
      max.post
   }
   bafs.L           <- mapply(FUN = bayes.fit, rows.x,
                         MoreArgs = list(mat = mufreq.tab, model.pts = model.pts),
                                         SIMPLIFY = FALSE)
   bafs.L           <- do.call(rbind, bafs.L)
   colnames(bafs.L) <- c("CNt", "A", "B", "L")
   bafs.L 
}

shannon.types <- function(types.mat) {
   data.types <- apply(types.mat, 1, FUN = function(x) paste(x, collapse = "_"))
   tab.types <- table(data.types)/length(data.types)
   tab.types <- tab.types * log2(tab.types)
   -sum(tab.types)
}

mufreq.model.fit <- function(mufreq, depth.ratio, weight.mufreq = 10, weight.ratio = 10, avg.depth.ratio,
                             cellularity.range = c(0.3,0.8), dna.content.range = c(0.5,4),
                             by.c = 0.01, by.p = 0.01, mc.cores = 2, CNt.max = 7, CNr = 2) {
   require(parallel)

   no.mut  <- 1
   c.range <- seq( from = min(cellularity.range), to = max(cellularity.range), by = by.c)
   p.range <- seq( from = min(dna.content.range), to = max(dna.content.range), by = by.p)
   pc.comb <- function(x, dna.content, cellularity) {
    cbind(rep(dna.content[x], length(cellularity)),
          cellularity)
   }
   C.P     <- mapply( x = 1:length(p.range), FUN = pc.comb,
                    MoreArgs=list(dna.content = p.range, cellularity = c.range),
                    SIMPLIFY = FALSE)
   C.P     <- do.call(rbind, C.P)
   types   <- types.matrix(CNt.min = 1, CNt.max = CNt.max, CNr = CNr)
   types   <- types[types[,3] >= no.mut, ]

   fit.cp <- function(x, C.P. = C.P, mufreq. = mufreq, depth.ratio. = depth.ratio,
                      weight.mufreq. = weight.mufreq, weight.ratio. = weight.ratio, 
                      avg.depth.ratio. = avg.depth.ratio, CNt.min. = 1, CNt.max. = CNt.max, CNr. = CNr) {
      dna.content. <- C.P.[x, 1]
      cellularity. <- C.P.[x, 2]
      L.model <- mufreq.bayes(mufreq = mufreq., depth.ratio = depth.ratio., weight.mufreq = weight.mufreq., 
                              weight.ratio = weight.ratio., cellularity = cellularity., dna.content = dna.content.,
                              avg.depth.ratio = avg.depth.ratio., CNt.min = CNt.min., CNt.max = CNt.max., CNr = CNr.)
      L.sum <- sum(L.model[,4])
      c(dna.content.,  cellularity., L.sum)  
   }

   bayes.res <- lapply_pb(X=1:nrow(C.P), FUN = fit.cp, mc.cores = mc.cores)
   #bayes.res <- lapply(X=1:nrow(C.P), FUN = fit.cp)

   bayes.res <- do.call(rbind, bayes.res)
   colnames(bayes.res) <- c("dna.content", "cellularity", "L")
   bayes.res
} 

baf.model.fit <- function(Bf, depth.ratio, weight.Bf = 10, weight.ratio = 10,
                          cellularity.range = c(0.3,1), dna.content.range = c(0.7,4),
                          by.c = 0.01, by.p = 0.01, avg.depth.ratio, mc.cores = 2,
                          CNt.max = 7, CNr = 2) {

   require(parallel)
   c.range <- seq( from = min(cellularity.range), to = max(cellularity.range), by = by.c)
   p.range <- seq( from = min(dna.content.range), to = max(dna.content.range), by = by.p)
   pc.comb <- function(x, dna.content, cellularity) {
    cbind(rep(dna.content[x], length(cellularity)),
          cellularity)
   }
   C.P   <- mapply( x = 1:length(p.range), FUN = pc.comb,
                    MoreArgs=list(dna.content = p.range, cellularity = c.range),
                    SIMPLIFY = FALSE)
   C.P   <-  do.call(rbind, C.P)
   
   fit.cp <- function(x, C.P. = C.P, Bf. = Bf, depth.ratio. = depth.ratio, weight.Bf. = weight.Bf, weight.ratio. = weight.ratio,
                      avg.depth.ratio. = avg.depth.ratio, CNt.min. = 0, CNt.max. = CNt.max, CNr. = CNr) {
      dna.content <- C.P.[x, 1]
      cellularity <- C.P.[x, 2]
      L.model <- baf.bayes(Bf = Bf., depth.ratio = depth.ratio., weight.Bf = weight.Bf., weight.ratio = weight.ratio.,
                           cellularity = cellularity, dna.content = dna.content,
                           avg.depth.ratio = avg.depth.ratio., CNt.min = CNt.min., 
                           CNt.max = CNt.max., CNr = CNr.)
      L.sum <- sum(L.model[,4])
      c(dna.content,  cellularity, L.sum)  
   }
   bayes.res <- lapply_pb(X=1:nrow(C.P), FUN = fit.cp, mc.cores = mc.cores)
   #bayes.res <- lapply(X=1:nrow(C.P), FUN = fit.cp)
   bayes.res <- do.call(rbind, bayes.res)
   colnames(bayes.res) <- c("dna.content", "cellularity", "L")
   bayes.res
}

