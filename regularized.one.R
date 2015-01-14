works_with_R("3.1.2",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             reshape2="1.2.2",
             "Rdatatable/data.table@31cf650e6fe54fe76b73086e71ad4815bc7eeb8c",
             directlabels="2014.6.13",
             dplyr="0.4.0")

load("dp.peaks.intervals.RData")
load("dp.peaks.sets.RData")
load("dp.peaks.optimal.RData")
load("dp.peaks.matrices.RData")

source("interval-regression.R")

set.name <- "H3K4me3_PGP_immune"
chunk.list <- dp.peaks.intervals[[set.name]]
bad.count.mat <- sapply(chunk.list, function(L){
  colSums(!is.finite(L$all.features))
})
feature.is.good <- apply(bad.count.mat==0, 1, all)
rownames(bad.count.mat)[feature.is.good]
train.sets <- dp.peaks.sets[[set.name]]
set.i <- 1
testSet <- paste(set.name, "split", set.i)
train.validation <- train.sets[[set.i]]

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

e.list <- list()
w.list <- list()
s.list <- list()
for(seed in 1:10){
  print(seed)
  set.seed(seed)
  ## Randomly divide the train set into 1/2 train, 1/2 validation.
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

  train.features <- set.data$train$features[, feature.is.good]
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
  set.data.error <- list()
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
      sample.ids <- unique(optimal.dt$sample.id)
      for(one.id in sample.ids){
        optimal.one <- optimal.dt %>%
          filter(sample.id == one.id)
        setkey(optimal.one, min.log.lambda, max.log.lambda)
        penalty.one <- penalty.dt %>%
          filter(sample.id == one.id)
        setkey(penalty.one, penalty, penalty2)
        ## Why do we need this filter? Shouldn't foverlaps do this
        ## automatically?
        overlap.one <- foverlaps(penalty.one, optimal.one, nomatch=0) %>%
          select(regularization, min.log.lambda, penalty, max.log.lambda) %>%
            arrange(regularization)
        n.regularization <- length(fit$gamma.seq)
        stopifnot(nrow(overlap.one) == n.regularization)
      }
      
      overlap.dt <- foverlaps(penalty.dt, optimal.dt, nomatch=0)
      overlap.counts <- overlap.dt %>%
        group_by(sample.id) %>%
        summarise(count=n())
      stopifnot(overlap.counts$count == n.regularization)
      err.mat <- dp.peaks.matrices[[set.name]][[chunk.name]]$PeakSeg
      regions <- dp.peaks.matrices[[set.name]][[chunk.name]]$regions
      names(dimnames(err.mat))[2] <- "model.complexity"
      err.dt <- data.table(melt(err.mat, value.name="errors"))
      join.dt <- inner_join(err.dt, overlap.dt)
      stopifnot(nrow(join.dt) == nrow(overlap.dt))
      join.dt$regions <- regions[as.character(join.dt$sample.id)]
      error.list[[chunk.name]] <- data.table(set=tv, join.dt)
    }
    set.data.error[[tv]] <- do.call(rbind, error.list)
  }
  set.dt <- do.call(rbind, set.data.error)
  set.stats <- set.dt %>%
    group_by(set, regularization) %>%
    summarise(errors=sum(errors),
              regions=sum(regions)) %>%
    mutate(what="percent error",
           percent=errors/regions*100) %>%
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
    arrange(errors, arclength) %>%
    select(arclength, regularization, errors) %>%
    head(1)

  s.list[[seed]] <- data.table(seed, selected)

  best.coefs <- coef.tall %>%
    filter(regularization == selected$regularization)
  e.list[[seed]] <- data.table(seed, join.stats)

  nonzero.coefs <- best.coefs %>%
    filter(value != 0)

  w.list[[seed]] <- data.table(seed, coef.tall)

  lasso <-
    ggplot()+
    geom_vline(aes(xintercept=arclength), data=selected, color="grey")+
    xlab("model complexity (L1 norm of weights)")+
  geom_line(aes(arclength, value, group=feature, color=feature),
            data=coef.tall)+
  geom_point(aes(arclength, value, group=feature, color=feature),
            data=subset(coef.tall, value != 0), pch=1)+
  geom_line(aes(arclength, percent, group=set, linetype=set),
            data=join.stats, show_guide=FALSE)+
  geom_dl(aes(arclength, percent, label=set),
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
selected.dt <- do.call(rbind, s.list)
r <- selected.dt$regularization
hist(log(r));points(log(r), rep(0, length(r)))
mean.log.reg <- mean(log(r))
abline(v=mean.log.reg)
hist(r);points(r, rep(0, length(r)))
mean.reg <- mean(r)
abline(v=mean.reg)
abline(v=exp(mean.log.reg), col="red")

## Fit model to entire train set.
set.chunks <-
  list(train=train.validation,
       test=names(chunk.list)[! names(chunk.list) %in% train.validation])
stopifnot(sum(sapply(set.chunks, length)) == length(unique(unlist(set.chunks))))
with(set.chunks, stopifnot(length(intersect(train, test)) == 0))

## actually we don't need to do the test set here...
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

train.features <- set.data$train$features[, feature.is.good]
feature.sd <- apply(train.features, 2, sd)
keep <- is.finite(feature.sd) & feature.sd > 0
train.mat <- train.features[, keep]
while(length(dup.names <- dup.vars(train.mat))){
  is.dup <- colnames(train.mat) %in% dup.names
  train.mat <- train.mat[, !is.dup]
}

small.mat <- train.mat[, c("log.bases", "log.unweighted.quartile.100%")]

fits <-
  list(L1.reg.log=smooth.interval.regression(train.mat, set.data$train$intervals,
                             gamma=exp(mean.log.reg),
                             calc.grad=calc.grad.list$square,
                             calc.loss=calc.loss.list$square),
       L1.reg=smooth.interval.regression(train.mat, set.data$train$intervals,
                             gamma=mean.reg,
                             calc.grad=calc.grad.list$square,
                             calc.loss=calc.loss.list$square),
       log.bases.log.max=smooth.interval.regression(small.mat,
         set.data$train$intervals,
         calc.grad=calc.grad.list$square,
         calc.loss=calc.loss.list$square))

error.list <- list()
for(chunk.name in set.chunks$test){
  f.mat <- chunk.list[[chunk.name]]$all.features
  for(model.name in names(fits)){
    fit <- fits[[model.name]]
    penalty.mat <- fit$predict(f.mat)
    stopifnot(ncol(penalty.mat) == 1)
    colnames(penalty.mat) <- "penalty"
    penalty.dt <- data.table(sample.id=rownames(penalty.mat), penalty.mat) %>%
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
    
    overlap.dt <- foverlaps(penalty.dt, optimal.dt, nomatch=0)
    stopifnot(nrow(overlap.dt) == nrow(penalty.dt))
    err.mat <- dp.peaks.matrices[[set.name]][[chunk.name]]$PeakSeg
    regions <- dp.peaks.matrices[[set.name]][[chunk.name]]$regions
    names(dimnames(err.mat))[2] <- "model.complexity"
    err.dt <- data.table(melt(err.mat, value.name="errors"))
    join.dt <- inner_join(err.dt, overlap.dt, c("sample.id", "model.complexity"))
    join.dt$regions <- regions[as.character(join.dt$sample.id)]
    stopifnot(nrow(join.dt) == nrow(overlap.dt))
    error.list[[paste(model.name, chunk.name)]] <-
      data.table(model.name, join.dt)
  }
}
test.err <- do.call(rbind, error.list)

test.err.stats <- test.err %>%
  group_by(model.name) %>%
  summarise(problems=n(),
            errors=sum(errors),
            regions=sum(regions)) %>%
  mutate(percent=errors/regions*100)
stopifnot(test.err$problems == nrow(set.data$test$features))

regularized.one <-
  list(weights=do.call(rbind, w.list),
       errors=do.call(rbind, e.list),
       testError.stats=test.err.stats,
       testError=test.err)

save(regularized.one, file="regularized.one.RData")
