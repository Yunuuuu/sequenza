mufreq.ccf <- function(cellularity, mufreq, Mt, CNt, CNn = 2, normal.ploidy = 2) {
  rel.freq    <- mufreq/cellularity
  tumor.comp  <- cellularity * CNt
  normal.comp <- normal.ploidy * (1 - cellularity)
  mutation.multiplicity <- rel.freq * (tumor.comp + normal.comp)
  mutation.multiplicity/Mt
}