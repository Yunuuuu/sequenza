mufreq.ccf <- function(cellularity, mufreq, Mt, CNt, CNn = 2, normal.ploidy = 2) {
  rel.freq    <- mufreq/cellularity
  tumor.comp  <- cellularity * CNt
  normal.comp <- normal.ploidy * (1 - cellularity)
  mutation.multiplicity <- rel.freq * (tumor.comp + normal.comp)
  mutation.multiplicity/Mt
}

sequenza.subclonal <- function(sequenza.extract, cp.table = NULL, sample.id, out.dir = getwd(),
                             cellularity = NULL, ploidy = NULL, female = TRUE, CNt.max = 20,
                             ratio.priority = FALSE, XY = c(X = "X", Y = "Y"), chromosome.list = 1:24){
  if(!file.exists(out.dir)) {
    dir.ok <- dir.create(path = out.dir, recursive = TRUE)
    if(!dir.ok) stop('Directory does not exist and cannot be created: ', out.dir)
  }
  makeFilename <- function(x) file.path(out.dir, paste(sample.id, x, sep = '_'))
  
  muts.file <- makeFilename("mutations_ccf.txt")
  segs.file <- makeFilename("segments_floats.txt")

  seg.tab     <- do.call(rbind, sequenza.extract$segments[chromosome.list])
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
  types <- types.matrix(CNt.min = 0, CNt.max = CNt.max, CNn = 2)
  model <- model.points(cellularity = cellularity, ploidy = ploidy, types = types, avg.depth.ratio = avg.depth.ratio)
  model <- cbind(types, model)
  model.ratio     <- unique(model[, c("CNn","CNt","depth.ratio")])
  CNt.to.ratio    <- setNames(model.ratio$depth.ratio, model.ratio$CNt)
  mean.ratio.step <- mean(diff(model.ratio$depth.ratio))
  
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
  model.vect <- CNt.to.ratio[as.character(seg.res$CNt)]
  float.vect <- (seg.res$depth.ratio - model.vect)/mean.ratio.step
  float.vect[float.vect < 0 && seg.res$CNt == 0] <- 0
  seg.res    <- cbind(seg.res, CNt.float = seg.res$CNt + round(float.vect, 3))
  if (!female){
    if (sum(segs.is.xy) >= 1) {
      cn.alleles  <- baf.bayes(Bf = NA, CNt.max = CNt.max,
                               depth.ratio = seg.tab$depth.ratio[segs.is.xy],
                               cellularity = cellularity, ploidy = ploidy,
                               avg.depth.ratio = avg.depth.ratio, sd.ratio = seg.tab$sd.ratio[segs.is.xy],
                               weight.ratio = seg.len[segs.is.xy], sd.Bf = NA,
                               weight.Bf = NA, ratio.priority = ratio.priority, CNn = 1)   
      
      seg.xy     <- cbind(seg.tab[segs.is.xy, ], cn.alleles)
      model.vect <- CNt.to.ratio[as.character(seg.xy$CNt)]
      float.vect <- (seg.xy$depth.ratio - model.vect)/mean.ratio.step
      float.vect[float.vect < 0] <- 0
      seg.xy     <- cbind(seg.xy, CNt.float = seg.xy$CNt + round(float.vect,3))      
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
    
    get.ccf <- function(x) {
      Mta <- as.numeric(x['Mt'])
      if (Mta == 0) {Mta = 1}
      mufreq.ccf(cellularity = cellularity,
                 mufreq = as.numeric(x['F']),
                 CNt = as.numeric(x['CNt']),
                 Mt = Mta,
                 CNn = as.numeric(x['CNn']))
    }
    
    ccfs <-  apply( mut.res, 1, function(x) get.ccf(x))
    mut.res <- cbind(mut.res, CCF = round(ccfs, 3))
    write.table(mut.res, file = muts.file,
                col.names = TRUE, row.names = FALSE, sep = "\t")   
  }
}


