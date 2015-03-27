works_with_R("3.1.3",
             directlabels="2014.6.13",
             data.table="1.9.4",
             "tdhock/PeakSegDP@c104b7681745851b286c0e89a7abed277277e415",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             dplyr="0.4.0")

load("oracle.intervals.RData")
load("oracle.optimal.RData")
load("dp.peaks.matrices.RData")
load("dp.peaks.sets.RData")

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

error.result.list <- list()
for(set.name in c("H3K36me3_AM_immune", "H3K4me3_TDH_other")){
  chunk.list <- oracle.intervals[[set.name]]
  bad.count.mat <- sapply(chunk.list, function(L){
    colSums(!is.finite(L$all.features))
  })
  feature.is.good <- apply(bad.count.mat==0, 1, all)
  all.chunks <- names(chunk.list)
  set.i <- 1
  all.train.validation <- dp.peaks.sets[[set.name]][[set.i]]
  test.chunks <- all.chunks[! all.chunks %in% all.train.validation]
  ##for(set.i in seq_along(train.sets)){
  testSet <- paste(set.name, "split", set.i)

  train.sizes <- seq(2, length(all.train.validation), by=2)
  for(train.size in train.sizes){
    train.validation <- all.train.validation[1:train.size]

    n.folds <- if(length(train.validation) == 2) 2 else 3
    fold.id <- rep(1:n.folds, l=length(train.validation))

    e.list <- list()
    w.list <- list()
    s.list <- list()
    for(validation.fold in 1:n.folds){
      is.validation <- fold.id == validation.fold
      ## Randomly divide the train set into 1/2 train, 1/2 validation.
      set.chunks <-
        list(train=train.validation[!is.validation],
             validation=train.validation[is.validation])

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
                                        L0=ncol(train.mat)*1.5,
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
          optimal.by.sample <- oracle.optimal[[set.name]][[chunk.name]]
          err.mat <- dp.peaks.matrices[[set.name]][[chunk.name]]$PeakSeg
          regions <- dp.peaks.matrices[[set.name]][[chunk.name]]$regions
          for(regularization.i in 1:ncol(penalty.mat)){
            regularization <- fit$gamma.seq[[regularization.i]]
            penalty.vec <- penalty.mat[, regularization.i]
            for(sample.id in names(penalty.vec)){
              penalty <- penalty.vec[[sample.id]]
              optimal <- optimal.by.sample[[sample.id]]
              optimal.peaks <- subset(optimal, {
                min.log.lambda < penalty & penalty < max.log.lambda
              })$peaks
              stopifnot(length(optimal.peaks) == 1)
              errors <- err.mat[sample.id, paste(optimal.peaks)]
              error.list[[paste(chunk.name, regularization.i, sample.id)]] <- 
                data.table(set=tv, chunk.name, regularization.i, regularization,
                           sample.id, errors, regions=regions[[sample.id]])
            }#sample.id
          }#regularization.i
        }#chunk.name
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
      arclength <-
        unique(data.frame(coef.tall)[, c("regularization", "arclength")])
      rownames(arclength) <- arclength$regularization
      set.stats$arclength <-
        arclength[paste(set.stats$regularization), "arclength"]

      selected <- set.stats %>%
        filter(set=="validation") %>%
        group_by() %>%
        arrange(errors, regularization) %>%
        select(regularization, errors) %>%
        head(1)

      s.list[[validation.fold]] <- data.table(validation.fold, selected)

      best.coefs <- coef.tall %>%
        filter(regularization == selected$regularization)
      e.list[[validation.fold]] <- data.table(validation.fold, set.stats)

      nonzero.coefs <- best.coefs %>%
        filter(value != 0)

      w.list[[validation.fold]] <- data.table(validation.fold, coef.tall)

      lasso <-
        ggplot()+
          ggtitle(paste(testSet, "validation set", validation.fold))+
          geom_vline(aes(xintercept=-log10(regularization)),
                     data=selected, color="grey")+
          xlab("model complexity (L1 norm of weights)")+
          geom_line(aes(-log10(regularization), value,
                        group=feature, color=feature),
                    data=coef.tall, show_guide=FALSE)+
          geom_point(aes(-log10(regularization),
                         value, group=feature, color=feature),
                     data=subset(coef.tall, value != 0),
                     pch=1, show_guide=FALSE)+
          geom_line(aes(-log10(regularization), percent, group=set, linetype=set),
                    data=set.stats, show_guide=FALSE)+
          geom_dl(aes(-log10(regularization), percent, label=set),
                  data=set.stats, method="lines2")+
          theme_bw()+
          theme(panel.margin=grid::unit(0, "cm"))+
          facet_grid(what ~ ., scales="free")
      print(lasso)
    }#validation.fold

    selected.dt <- do.call(rbind, s.list)
    r <- selected.dt$regularization
    mean.reg <- mean(r)

    ## Fit model to entire train set.
    set.chunks <-
      list(train=train.validation,
           test=names(chunk.list)[! names(chunk.list) %in% train.validation])
    stopifnot(sum(sapply(set.chunks, length)) ==
                length(unique(unlist(set.chunks))))
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
      list(L1.reg=smooth.interval.regression(train.mat,
             set.data$train$intervals,
             gamma=mean.reg,
             calc.grad=calc.grad.list$square,
             calc.loss=calc.loss.list$square),
           log.bases.log.max=smooth.interval.regression(small.mat,
             set.data$train$intervals,
             calc.grad=calc.grad.list$square,
             calc.loss=calc.loss.list$square))

    for(chunk.name in set.chunks$test){
      optimal.by.sample <- oracle.optimal[[set.name]][[chunk.name]]
      err.mat <- dp.peaks.matrices[[set.name]][[chunk.name]]$PeakSeg
      regions <- dp.peaks.matrices[[set.name]][[chunk.name]]$regions
      f.mat <- chunk.list[[chunk.name]]$all.features
      for(model.name in names(fits)){
        fit <- fits[[model.name]]
        penalty.mat <- fit$predict(f.mat)
        stopifnot(ncol(penalty.mat) == 1)
        for(sample.id in rownames(penalty.mat)){
          penalty <- penalty.mat[sample.id, ]
          optimal <- optimal.by.sample[[sample.id]]
          optimal.peaks <- subset(optimal, {
            min.log.lambda < penalty & penalty < max.log.lambda
          })$peaks
          stopifnot(length(optimal.peaks) == 1)
          errors <- err.mat[sample.id, paste(optimal.peaks)]
          error.result.list[[paste(chunk.name, model.name, sample.id)]] <- 
            data.table(set.name, train.size, chunk.name, model.name,
                       sample.id, errors, regions=regions[[sample.id]])
        }#sample.id
      }#model.name
    }#chunk.name in test set
  }#train.size
}#set.name

error.results <- do.call(rbind, error.result.list)
error.stats <- error.results %>%
  group_by(set.name, train.size, model.name) %>%
  summarise(errors=sum(errors),
            regions=sum(regions)) %>%
  mutate(percent=errors/regions*100)

ggplot()+
  geom_line(aes(train.size, percent, color=model.name),
            data=error.stats)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(set.name ~ .)

## TODO: average over 6 train/test splits.
