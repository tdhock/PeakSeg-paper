works_with_R("3.1.2",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             reshape2="1.2.2",
             data.table="1.9.4",
             directlabels="2014.6.13",
             dplyr="0.4.0")

load("dp.peaks.intervals.RData")
load("dp.peaks.sets.RData")
load("dp.peaks.optimal.RData")
load("dp.peaks.matrices.RData")

source("interval-regression.R")

g <- as.matrix(expand.grid(log.max.coverage=seq(2, 8, l=100),
                           log.total.weight=seq(6, 15, l=100)))
dp.peaks.grid.list <- list()
dp.peaks.polygon.list <- list()
dp.peaks.segment.list <- list()
dp.peaks.regression <- NULL
dp.peaks.prediction.list <- list()
dp.peaks.roc <- list()
dp.peaks.roc.chosen <- list()
getPolygons <- function(df){
  indices <- with(df, chull(log.max.coverage, log.total.weight))
  df[indices, ]
}
set.name <- "H3K4me3_PGP_immune"
chunk.list <- dp.peaks.intervals[[set.name]]
train.sets <- dp.peaks.sets[[set.name]]
set.i <- 1
testSet <- paste(set.name, "split", set.i)
train.validation <- train.sets[[set.i]]

e.list <- list()
w.list <- list()
for(seed in 1:4){
  print(seed)
  set.seed(seed)
  train.chunks <- sample(train.validation, length(train.validation)/2)
  set.chunks <-
    list(train=train.chunks,
         validation=train.validation[!train.validation %in% train.chunks])

  set.data <- list()
  for(tv in names(set.chunks)){
    f <- NULL
    i <- NULL
    chunk.names <- set.chunks[[tv]]
    for(chunk.name in chunk.names){
      data.list <- chunk.list[[chunk.name]]
      f <- rbind(f, data.list$all.features)
      i <- rbind(i, data.list$intervals)
    }
    set.data[[tv]] <- list(features=f, intervals=i)
  }

  dup.vars <- function(features){
    for(feature.i in 1:(ncol(features)-1)){
      feature.vec <- features[, feature.i]
      other.i <- (feature.i+1):ncol(features)
      other.mat <- features[, other.i, drop=FALSE]
      residual <- colSums(abs(feature.vec - other.mat))
      is.dup <- residual == 0
      if(any(is.dup)){
        return(colnames(other.mat)[is.dup])
      }
    }
    character(0)
  }
  train.features <- set.data$train$features
  feature.sd <- apply(train.features, 2, sd)
  keep <- is.finite(feature.sd) & feature.sd > 0
  train.mat <- train.features[, keep]
  while(length(dup.names <- dup.vars(train.mat))){
    is.dup <- colnames(train.mat) %in% dup.names
    train.mat <- train.mat[, !is.dup]
  }
  fit <-
    regularized.interval.regression(train.mat, set.data$train$intervals,
                                    calc.grad=calc.grad.list$square,
                                    calc.loss=calc.loss.list$square)
  ## Use these models on the train and validation sets.
  set.data <- list()
  for(tv in names(set.chunks)){
    chunk.names <- set.chunks[[tv]]
    error.list <- list()
    for(chunk.name in chunk.names){
      f.mat <- chunk.list[[chunk.name]]$all.features
      penalty.mat <- fit$predict(f.mat)
      colnames(penalty.mat) <- fit$gamma.seq
      names(dimnames(penalty.mat)) <- c("sample.id", "regularization")
      penalty.dt <- data.table(melt(penalty.mat, value.name="penalty")) %>%
        mutate(penalty2=penalty)
      setkey(penalty.dt, sample.id, penalty, penalty2)
      sample.list <- dp.peaks.optimal[[set.name]][[chunk.name]]
      optimal.dt <-
        do.call(rbind, lapply(names(sample.list), function(sample.id){
          data.table(sample.id, sample.list[[sample.id]]) %>%
            mutate(min.log.lambda=ifelse(min.log.lambda == -Inf,
                     -1e10, min.log.lambda),
                   max.log.lambda=ifelse(max.log.lambda == Inf,
                     1e10, max.log.lambda)) # bug in foverlaps when -Inf.
        }))
      setkey(optimal.dt, sample.id, min.log.lambda, max.log.lambda)
      one.id <- "McGill0011"
      optimal.one <- optimal.dt %>%
        filter(sample.id == one.id)
      setkey(optimal.one, min.log.lambda, max.log.lambda)
      penalty.one <- penalty.dt %>%
        filter(sample.id == one.id)
      setkey(penalty.one, penalty, penalty2)
      ## Why do we need this filter? Shouldn't foverlaps do this
      ## automatically?
      overlap.one <- foverlaps(penalty.one, optimal.one) %>%
        filter(min.log.lambda < penalty & penalty < max.log.lambda) %>%
          select(regularization, min.log.lambda, penalty, max.log.lambda) %>%
            arrange(regularization)
      n.regularization <- length(fit$gamma.seq)
      stopifnot(nrow(overlap.one) == n.regularization)
      overlap.dt <- foverlaps(penalty.dt, optimal.dt) %>%
        filter(min.log.lambda < penalty & penalty < max.log.lambda)
      overlap.counts <- overlap.dt %>%
        group_by(sample.id) %>%
          summarise(count=n())
      stopifnot(overlap.counts$count == n.regularization)
      err.mat <- dp.peaks.matrices[[set.name]][[chunk.name]]$PeakSeg
      names(dimnames(err.mat))[2] <- "model.complexity"
      err.dt <- data.table(melt(err.mat, value.name="errors"))
      join.dt <- inner_join(err.dt, overlap.dt)
      stopifnot(nrow(join.dt) == nrow(overlap.dt))
      error.list[[chunk.name]] <- data.table(set=tv, join.dt)
    }
    set.data[[tv]] <- do.call(rbind, error.list)
  }
  set.dt <- do.call(rbind, set.data)
  set.stats <- set.dt %>%
    group_by(set, regularization) %>%
      summarise(errors=sum(errors)) %>%
        mutate(what="error") %>%
          arrange(set, regularization)
  coef.mat <- fit$coefs
  colnames(coef.mat) <- fit$gamma.seq
  names(dimnames(coef.mat)) <- c("feature", "regularization")
  coef.tall <- melt(coef.mat[-1, ]) %>%
    group_by(regularization) %>%
      mutate(arclength=sum(abs(value)),
             what="weights")
  arclength <- unique(coef.tall[, c("regularization", "arclength")])
  join.stats <- inner_join(arclength, set.stats)

  selected <- join.stats %>%
    filter(set=="validation") %>%
    group_by() %>%
    filter(seq_along(errors) == which.min(errors)) %>%
    select(arclength, regularization, errors)

  best.coefs <- coef.tall %>%
    filter(regularization == selected$regularization)
  w.list[[seed]] <- data.table(seed, best.coefs)

  nonzero.coefs <- best.coefs %>%
    filter(value != 0)

  e.list[[seed]] <- data.table(seed, coef.tall)

  lasso <-
    ggplot()+
    geom_vline(aes(xintercept=arclength), data=selected, color="grey")+
    xlab("model complexity (L1 norm of weights)")+
  geom_line(aes(arclength, value, group=feature, color=feature),
            data=coef.tall)+
  geom_line(aes(arclength, errors, group=set, linetype=set),
            data=join.stats, show_guide=FALSE)+
  geom_dl(aes(arclength, errors, label=set),
          data=join.stats, method="lines2")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(what ~ ., scales="free")

  lasso.selected <- list(function(df, ...){
    subset(df, groups %in% nonzero.coefs$feature)
  }, "lasso.labels")
  lasso.dl <- direct.label(lasso, "lasso.selected")
  print(lasso.dl)
}
## All of this normalization code is already inside
## regularized.interval.regression.
apply(train.mat, 2, range)
train.names <- colnames(train.mat)
train.mean <- colMeans(train.mat)
train.sd <- feature.sd[keep]
train.mean.mat <-
  matrix(train.mean, nrow(train.mat), ncol(train.mat), byrow=TRUE)
train.sd.mat <-
  matrix(train.sd, nrow(train.mat), ncol(train.mat), byrow=TRUE)
train.norm <- (train.mat - train.mean.mat)/train.sd.mat
stopifnot(abs(colMeans(train.norm)) < 1e-5)
stopifnot(abs(1-apply(train.norm, 2, sd)) < 1e-5)

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
  optimal.list <- dp.peaks.optimal[[set.name]][[test.chunk]]
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
    optimal.df$fp <- fp[as.character(optimal.df$model.complexity)]
    optimal.df$tp <- tp[as.character(optimal.df$model.complexity)]
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
      as.character(pred.row$model.complexity)
    errors <- this.test$test.errors[[sample.i]] <-
      error.list$PeakSeg[sample.i, peaks]
    regions <- error.list$regions[sample.i]
    dp.peaks.regression <- rbind(dp.peaks.regression, {
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

dp.peaks.prediction <- do.call(rbind, dp.peaks.prediction.list)
dp.peaks.grid <- do.call(rbind, dp.peaks.grid.list)
dp.peaks.polygon <- do.call(rbind, dp.peaks.polygon.list)
dp.peaks.segment <- do.call(rbind, dp.peaks.segment.list)

save(dp.peaks.regression,
     dp.peaks.roc,
     dp.peaks.roc.chosen,
     dp.peaks.prediction,
     dp.peaks.grid,
     dp.peaks.polygon,
     dp.peaks.segment,
     file="dp.peaks.regression.RData")
