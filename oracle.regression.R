works_with_R("3.1.1",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             dplyr="0.2")

load("oracle.intervals.RData")
load("dp.peaks.sets.RData")
load("oracle.optimal.RData")
load("dp.peaks.matrices.RData")

source("interval-regression.R")

g <- as.matrix(expand.grid(log.max.coverage=seq(2, 8, l=100),
                           log.total.weight=seq(6, 15, l=100)))
dp.peaks.grid.list <- list()
dp.peaks.polygon.list <- list()
dp.peaks.segment.list <- list()
oracle.regression <- NULL
dp.peaks.prediction.list <- list()
dp.peaks.roc <- list()
dp.peaks.roc.chosen <- list()
getPolygons <- function(df){
  indices <- with(df, chull(log.max.coverage, log.total.weight))
  df[indices, ]
}
for(set.name in names(oracle.intervals)){
  chunk.list <- oracle.intervals[[set.name]]
  train.sets <- dp.peaks.sets[[set.name]]
  for(set.i in seq_along(train.sets)){
    testSet <- paste(set.name, "split", set.i)
    train.chunks <- train.sets[[set.i]]
    train.features <- NULL
    train.intervals <- NULL
    for(chunk.name in train.chunks){
      data.list <- chunk.list[[chunk.name]]
      train.features <- rbind(train.features, data.list$features)
      train.intervals <- rbind(train.intervals, data.list$intervals)
    }
    ##train.features <- train.features[, "log.max.coverage", drop=FALSE]
    fit <- regression.funs$square(train.features, train.intervals)
    dp.peaks.grid.list[[testSet]] <-
      data.frame(set.name, set.i,
                 testSet,
                 g, log.lambda=fit$predict(g))
    dp.peaks.prediction.list[[paste(testSet, "train")]] <-
      data.frame(set.name, set.i,
                 test.chunk=NA,
                 sample.id=NA,
                 testSet, set="train",
                 train.features,
                 log.lambda=fit$train.predict,
                 peaks=NA,
                 test.errors=0)
    ## Now use the model on the test chunks.
    test.chunks <- names(chunk.list)[!names(chunk.list) %in% train.chunks]
    test.diffs <- test.complex <- list()
    for(test.chunk in test.chunks){
      error.list <- dp.peaks.matrices[[set.name]][[test.chunk]]
      fp.list <- dp.peaks.matrices.fp[[set.name]][[test.chunk]]
      tp.list <- dp.peaks.matrices.tp[[set.name]][[test.chunk]]
      data.list <- chunk.list[[test.chunk]]
      log.lambda.hat <- fit$predict(data.list$features)
      optimal.list <- oracle.optimal[[set.name]][[test.chunk]]
      this.test <- data.frame(set.name, set.i,
                              test.chunk,
                              sample.id=names(optimal.list),
                              testSet, set="test",
                              data.list$features,
                              log.lambda=log.lambda.hat,
                              peaks=NA,
                              test.errors=NA,
                              row.names=NULL)
      for(sample.i in seq_along(optimal.list)){
        optimal.df <- optimal.list[[sample.i]]
        sample.id <- names(optimal.list)[[sample.i]]
        fp <- fp.list$PeakSeg[sample.id, ]
        tp <- tp.list$PeakSeg[sample.id, ]
        l <- log.lambda.hat[sample.i, ]
        optimal.df$min.thresh <- optimal.df$min.log.lambda-l
        optimal.df$max.thresh <- optimal.df$max.log.lambda-l
        optimal.df$fp <- fp[as.character(optimal.df$peaks)]
        optimal.df$tp <- tp[as.character(optimal.df$peaks)]
        most.complex <- optimal.df[1,] %.%
          select(fp, tp) %.%
          mutate(possible.tp=tp.list$possible.tp[[sample.id]],
                 possible.fp=fp.list$possible.fp[[sample.id]])
        optimal.diffs <- with(optimal.df, {
          data.frame(thresh=min.thresh[-1], diff.fp=diff(fp), diff.tp=diff(tp))
        })
        pred.row <-
          subset(optimal.df, min.log.lambda < l & l < max.log.lambda)
        ## actually it is OK and it makes sense to have non-monotonic
        ## ROC curves, since peaks may disappear as model complexity increases.
        ##stopifnot(optimal.diffs$diff.fp <= 0) 
        ##stopifnot(optimal.diffs$diff.tp <= 0)
        test.diffs[[paste(test.chunk, sample.id)]] <-
          filter(optimal.diffs, diff.fp != 0 | diff.tp != 0)
        test.complex[[paste(test.chunk, sample.id)]] <- most.complex
        stopifnot(nrow(pred.row)==1)
        peaks <- this.test$peaks[[sample.i]] <-
          as.character(pred.row$peaks)
        errors <- this.test$test.errors[[sample.i]] <-
          error.list$PeakSeg[sample.i, peaks]
        regions <- error.list$regions[sample.i]
        oracle.regression <- rbind(oracle.regression, {
          data.frame(set.name, set.i, test.chunk,
                     sample.id,
                     peaks, errors, regions)
        })
      }
      if(nrow(this.test) == 2){
        one.seg <- this.test[1,]
        one.seg$log.max.coverage.end <- this.test[2, "log.max.coverage"]
        one.seg$log.total.weight.end <- this.test[2, "log.total.weight"]
        dp.peaks.segment.list[[paste(testSet, test.chunk)]] <-
          one.seg
      }else if(nrow(this.test) > 2){
        dp.peaks.polygon.list[[paste(testSet, test.chunk)]] <-
          getPolygons(this.test)
      }else{
        stop(nrow(this.test), "rows")
      }
      dp.peaks.prediction.list[[paste(testSet, test.chunk)]] <-
        this.test
    }#test.chunk
    all.complex <- do.call(rbind, test.complex)
    sum.complex <- colSums(all.complex)
    all.diffs <- do.call(rbind, test.diffs) %.%
      arrange(thresh) %.%
      mutate(cum.tp=cumsum(diff.tp),
             cum.fp=cumsum(diff.fp),
             tp=cum.tp+sum.complex[["tp"]],
             fp=cum.fp+sum.complex[["fp"]],
             TPR=tp/sum.complex[["possible.tp"]],
             FPR=fp/sum.complex[["possible.fp"]])
    last <- tail(all.diffs, 1)
    chosen <- all.diffs %.%
      filter(thresh < 0) %.%
      tail(1)
    if(nrow(chosen) == 0){
      chosen <- all.diffs[1,]
    }
    stopifnot(last[, c("tp", "fp")] == 0)
    ggplot()+
      coord_equal()+
      geom_point(aes(FPR, TPR), data=chosen, pch=1)+
      geom_path(aes(FPR, TPR), data=all.diffs)
    stopifnot(diff(all.diffs$thresh) != 0)
    dp.peaks.roc[[paste(set.name, set.i)]] <-
      data.frame(set.name, set.i, all.diffs)
    dp.peaks.roc.chosen[[paste(set.name, set.i)]] <-
      data.frame(set.name, set.i, chosen)
  }
}

dp.peaks.prediction <- do.call(rbind, dp.peaks.prediction.list)
dp.peaks.grid <- do.call(rbind, dp.peaks.grid.list)
dp.peaks.polygon <- do.call(rbind, dp.peaks.polygon.list)
dp.peaks.segment <- do.call(rbind, dp.peaks.segment.list)

save(oracle.regression,
     dp.peaks.roc,
     dp.peaks.roc.chosen,
     dp.peaks.prediction,
     dp.peaks.grid,
     dp.peaks.polygon,
     dp.peaks.segment,
     file="oracle.regression.RData")
