theoretical.depth.ratio <- function(cellularity = 0.5, dna.index = 1, copy.number.ratio = 1.5, avg.depth.ratio = 1) {
   cellu.copy.term <- (1 - cellularity) + (copy.number.ratio * cellularity / dna.index)
   avg.depth.ratio * cellu.copy.term
}

theoretical.mufreq <- function(cellularity, CNr = 2, CNt = 2, Mt = 1) {
   copy.ratio <- CNt / CNr
   cell.fract <- 1 / cellularity
   mut.factor <- Mt / CNr
   mut.factor/(copy.ratio + cell.fract - 1)
}

types.matrix <- function(CNt.min = 1, CNt.max = 7, CNr = 2) {
   cn.ratio.vect <- seq(from = CNt.min / CNr, to =  CNt.max / CNr, by = 1 / CNr)
   CNt           <- cn.ratio.vect * CNr
   mut.comb      <- lapply(CNt, FUN = function(x) seq(from = 0, to = x))
   times.muts    <- sapply(mut.comb, length)
   data.frame(CNr = CNr, CNt = rep(CNt, times = times.muts), 
              Mt = unlist(mut.comb))
}

model.points <- function(cellularity = 0.5, dna.index = 1,
                         types = cbind(CNr = 2, CNt = 2, Mt = 1),
                         avg.depth.ratio = 1) {
   mufreqs     <-  theoretical.mufreq(cellularity = cellularity , CNr = types[, 1], CNt = types[, 2], Mt = types[, 3])
   depth.ratio <-  theoretical.depth.ratio(cellularity = cellularity, dna.index = dna.index,
                                       copy.number.ratio = types[, 2] / types[, 1],
                                       avg.depth.ratio = avg.depth.ratio)
   cbind(mufreqs,depth.ratio)
}

# theoretical.baf <- function(cellularity = 0.5, CNt = 2, B = 1, CNr = 2){
#    B.tot <- ((B * cellularity)  + (1 - cellularity)) / 
#             ((CNt * cellularity) + CNr*(1 - cellularity))
#    B.tot
# }

theoretical.baf <- function(CNr, CNt, cellularity) {
   alleles       <- seq(from = 1, to = CNt, by = 1)
   max.b <- function(CNt) {
      max.b.alleles <- CNt / 2
      if (CNt %% 2 != 0 ) {
         max.b.alleles <- trunc(max.b.alleles) 
      }
      max.b.alleles
   }
   fract.normal.alleles <- (1 - cellularity)
   res                  <- list()
   for (i in 1:length(alleles)) {
      max.b.alleles <- max.b(alleles[i])
      max.a.alleles <- alleles[i] - max.b.alleles
      decrements.b  <- seq(from = max.b.alleles, to = 0, by = -1)
      res[[i]]     <- list()
      for (n in 1:length(decrements.b)) {
         A.i <- (max.a.alleles + decrements.b[n])
         B.i <- (max.b.alleles - decrements.b[n])
         BAF <- (fract.normal.alleles + (cellularity * B.i)) / ((alleles[i] * cellularity) + (CNr * fract.normal.alleles))
         res[[i]][[n]] <- cbind(A = A.i, B = B.i, BAF = BAF, CNt = alleles[i])
      }
   }
   for (i in 1:length(res)) {
      res[[i]] <- do.call(rbind,res[[i]]) 
   }
   as.data.frame(do.call(rbind,res))
}

baf.model.points <- function (cellularity, dna.index, avg.depth.ratio,
                               CNr = 2, CNt.min = 1, CNt.max = 4) {
   mufreq.depth.ratio <- model.points(cellularity = cellularity, dna.index = dna.index, 
                                      types = cbind(CNr = CNr, CNt = CNt.min:CNt.max, Mt = 0),
                                      avg.depth.ratio = avg.depth.ratio)
   model.d.ratio      <- cbind(CNt = CNt.min:CNt.max, depth.ratio = mufreq.depth.ratio[, 2])
   model.baf          <- theoretical.baf(CNr = CNr, CNt = CNt.max, cellularity = cellularity)
   if (CNt.min == 0) {
      model.baf <- as.data.frame(rbind(c(0,0,0.5,0), model.baf))
   }
   model.pts          <- merge(model.baf, model.d.ratio)
   model.pts
}
