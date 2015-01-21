exactModelSelection <- function(cost, model.complexity, peaks) 
{
  stopifnot(is.numeric(cost))
  stopifnot(is.numeric(model.complexity))
  stopifnot(diff(model.complexity) > 0)
  stopifnot(diff(cost) < 0)
  stopifnot(length(cost) == length(model.complexity))
  n.models <- length(cost)
  Kmax <- model.complexity[n.models]
  Kcurrent <- Kmax
  Lcurrent <- 0
  vK <- Kmax
  vL <- 0
  vP <- peaks[n.models]
  i <- 2
  min.complexity <- model.complexity[1]
  while(Kcurrent > min.complexity) {
    is.smaller <- model.complexity < Kcurrent
    is.current <- model.complexity == Kcurrent
    smallerK <- model.complexity[is.smaller]
    smallerPeaks <- peaks[is.smaller]
    cost.term <- cost[is.current] - cost[is.smaller]
    complexity.term <- smallerK - model.complexity[is.current]
    lambdaTransition <- cost.term/complexity.term
    next.i <- which.min(lambdaTransition)
    Kcurrent <- smallerK[next.i]
    Lcurrent <- min(lambdaTransition)
    vL[i] <- Lcurrent
    vK[i] <- Kcurrent
    vP[i] <- smallerPeaks[next.i]
    i <- i + 1
  }
  L <- log(vL)
  data.frame(min.log.lambda = L,
             max.log.lambda = c(L[-1], Inf),
             model.complexity = vK,
             peaks=vP,
             min.lambda = vL,
             max.lambda = c(vL[-1], Inf))
}


load("dp.peaks.matrices.RData")

oracle.optimal <- list()
for(set.name in names(dp.peaks.matrices)){
  chunk.list <- dp.peaks.matrices[[set.name]]
  for(chunk.name in names(chunk.list)){
    err.mat <- chunk.list[[chunk.name]]$PeakSeg
    model.file <- sprintf("data/%s/dp.model.RData", chunk.name)
    load(model.file)
    counts.file <- sprintf("data/%s/counts.RData", chunk.name)
    ## load(counts.file)
    ## counts.list <- split(counts, counts$sample.id)
    for(sample.id in rownames(err.mat)){
      cat(sprintf("%s %s\n", sample.id, chunk.name))
      sample.model <- dp.model[[sample.id]]
      ## With an unconstrained model the loss should always decrease
      ## with model complexity. However with our constrained model
      ## sometimes it increases, for example "H3K36me3_AM_immune/20"
      ## "McGill0004" in this case we consider only the less complex
      ## models.
      loss.df <- sample.model$error
      while(length(next.bigger <- which(diff(loss.df$error) > 0))){
        loss.df <- loss.df[-(next.bigger[1]+1), ]
      }
      bases <- with(sample.model$segments[1,], chromEnd-chromStart)
      in.sqrt <- 1.1 + log(bases / loss.df$segments)
      in.square <- 1 + 4 * sqrt(in.sqrt)
      complexity <- in.square * in.square * loss.df$segments
      exact.df <- exactModelSelection(loss.df$error, complexity, loss.df$peaks)
      if(FALSE){ # visual check of exact function with inexact grid search.
        log.lambda.seq <- 
          seq(exact.df$max.log.lambda[1]-1,
              exact.df$min.log.lambda[nrow(exact.df)]+1,
              l=200)
        grid.peaks <- sapply(log.lambda.seq, function(log.lambda){
          with(loss.df, {
            peaks[which.min(error + exp(log.lambda) * complexity)]
          })
        })
        grid.df <- data.frame(peaks=grid.peaks, log.lambda=log.lambda.seq)
        ggplot()+
          geom_point(aes(log.lambda, peaks), data=grid.df, pch=1, color="red")+
            geom_segment(aes(min.log.lambda, peaks,
                             xend=max.log.lambda, yend=peaks),
                         data=exact.df)
      }
      oracle.optimal[[set.name]][[chunk.name]][[sample.id]] <- exact.df
    }
  }
}

save(oracle.optimal, file="oracle.optimal.RData")
