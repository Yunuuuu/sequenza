subclonal.matrix <- function(mut.tab, cellularity = seq(0.1, 1, 0.05), dna.index, avg.depth.ratio, mc.cores = 2){
  
  mut.types.list <- lapply(X = 1:nrow(mut.tab),
                           FUN = function(x) { 
                             types.matrix(CNr = mut.tab[x, 'CNr'],
                                          CNt.min = mut.tab[x, 'CNt'],
                                          CNt.max = mut.tab[x, 'CNt'])})
  mut.cloanlity <- function(F, depth.t, types, cellularity, dna.index, avg.depth.ratio) {
    theorethic <- model.points(cellularity = cellularity,
                               dna.index = dna.index,
                               types = types,
                               avg.depth.ratio = avg.depth.ratio)
    max(mufreq.dpois(mufreq = F, mufreq.model = theorethic[, 1], depth.t = depth.t),na.rm = TRUE)
  }
  res <- mclapplyPb (mc.cores = mc.cores, X = 1:nrow(mut.tab),
                     FUN = function (i) {
                             sapply(X = cellularity, FUN = function(x) {
                                mut.cloanlity(F = mut.tab$F[i],
                                           depth.t = mut.tab$good.s.reads[i],
                                           types = mut.types.list[[i]],
                                           cellularity = x, dna.index = dna.index,
                                           avg.depth.ratio = avg.depth.ratio)
                           })
                         })
  res <- do.call(rbind, res)
  res/rowSums(res)
}
