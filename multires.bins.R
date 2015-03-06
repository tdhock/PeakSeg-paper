works_with_R("3.1.2",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             reshape2="1.2.2",
             data.table="1.9.4",
             "tdhock/PeakSegDP@45e95f50b957d82965dc55b1892a0e1e2516b649",
             dplyr="0.4.0")

load("dp.peaks.sets.RData")

pick.best.index <- function
### Minimizer for local models, described in article section 2.3
### "Picking the optimal model"
(err
### Vector of errors to minimize.
 ){
  nparam <- length(err)
  candidates <- which(err==min(err))
  if(length(err)==1)return(candidates)
  st <- abs(median(candidates)-candidates)
  middle <- candidates[which.min(st)]
  if(all(diff(err)==0))return(middle)
  if(nparam %in% candidates && 1 %in% candidates){
    cat("Warning: strange error profile, picking something near the center\n")
    print(as.numeric(err))
    d <- diff(candidates)>1
    if(any(d)){
      which(d)[1]
    }else{
      middle
    }
  }else if(1 %in% candidates){
    max(candidates)
  }else if(nparam %in% candidates){
    min(candidates)
  }else {
    middle
  }
### Integer index of the minimal error.
}

all.best.list <- list()
for(set.name in names(dp.peaks.sets)){
  train.sets <- dp.peaks.sets[[set.name]]
  label.dir <- paste0("~/genomelabels/", set.name)
  RData.files <- Sys.glob(file.path(label.dir, "*.RData"))
  names(RData.files) <- sub(".RData$", "", basename(RData.files))
  RData.by.sample <- list()
  e <- new.env()
  fl.list <- list()
  error.list <- list()
  for(sample.id in names(RData.files)){
    RData.file <- RData.files[[sample.id]]
    objs <- load(RData.file, e)
    RData.by.sample[[sample.id]] <- as.list(e)
    error.list[[sample.id]] <- data.table(sample.id, e$errors)
    to.check <- e$errors %>%
      group_by(bases.per.bin) %>%
      summarise(total.weight=sum(total.weight))
    exp.weight <- rep(to.check$total.weight[1], nrow(to.check))
    stopifnot(all.equal(exp.weight,
                        to.check$total.weight))
    for(bases.per.bin in names(e$features.limits)){
      fl <- e$features.limits[[bases.per.bin]]
      chunk.list <- fl$limits
      for(chunk.id in names(chunk.list)){
        limits <- chunk.list[[chunk.id]]
        chunk.name <- paste0(set.name, "/", chunk.id)
        features <- fl$features[rownames(limits), , drop=FALSE ]
        rownames(features) <- rownames(limits) <-
          paste(chunk.name, bases.per.bin, rownames(limits))
        fl.list[[bases.per.bin]][[chunk.name]][[sample.id]] <-
          list(features=features, limits=limits)
      }
    }
  }
  all.errors <- do.call(rbind, error.list)
  all.errors[, chunk.name := paste0(set.name, "/", chunk.id)]
  weights.by.resolution <- all.errors %>%
    group_by(bases.per.bin) %>%
    summarise(total.weight=sum(total.weight))
  max.total.weight <- max(weights.by.resolution$total.weight)
  complete.resolutions <- weights.by.resolution %>%
    filter(total.weight == max.total.weight)
  for(set.i in seq_along(train.sets)){
    testSet <- paste(set.name, "split", set.i)
    train.validation <- train.sets[[set.i]]

    error.by.sample <- 
    all.errors %>%
      filter(bases.per.bin %in% complete.resolutions$bases.per.bin,
             chunk.name %in% train.validation) %>%
      group_by(bases.per.bin, sample.id) %>%
      summarise(weighted.error=sum(weighted.error),
                total.weight=sum(total.weight))

    each.sample <- 
    ggplot(error.by.sample,
           aes(bases.per.bin, weighted.error, group=sample.id))+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free")+
    geom_line()+
    geom_point()+
    scale_x_log10()

    error.wide <- error.by.sample %>%
      group_by(bases.per.bin) %>%
      summarise(weighted.error=sum(weighted.error),
                total.weight=sum(total.weight)) %>%
      mutate(percent=weighted.error/total.weight*100)
    print(error.wide)
    error.tall <-
      melt(error.wide, measure.vars=c("weighted.error", "percent"))

    best.error <- error.tall %>%
      group_by(variable) %>%
      filter(seq_along(value) == which.min(value))

    ggplot(NULL, aes(bases.per.bin, value))+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(variable ~ ., scales="free")+
    geom_line(data=error.tall)+
    geom_point(data=best.error)+
    scale_x_log10()

    all.best.list[[testSet]] <- data.table(set.name, testSet, best.error)

    best.bases.per.bin <- paste(best.error$bases.per.bin[1])
    best.fl <- fl.list[[best.bases.per.bin]]

    n.folds <- if(length(train.validation) == 2) 2 else 3
    fold.id <- rep(1:n.folds, l=length(train.validation))

    gamma.vals <- list()
    for(validation.fold in 1:n.folds){
      is.validation <- fold.id == validation.fold
      ## Randomly divide the train set into 1/2 train, 1/2 validation.
      set.chunks <- 
        list(train=train.validation[!is.validation],
             validation=train.validation[is.validation])
      set.data <- list()
      for(tv in c("train", "validation")){
        for(data.type in c("features", "limits")){
          data.list <- list()
          for(set.chunk in set.chunks[[tv]]){
            sample.list <- best.fl[[set.chunk]]
            for(sample.id in names(sample.list)){
              fl <- sample.list[[sample.id]]
              data.list[[paste(set.chunk, sample.id)]] <- fl[[data.type]]
            }
          }
          set.data[[tv]][[data.type]] <- do.call(rbind, data.list)
        }
      }
      fmat <- with(set.data$train, {
        stopifnot(rownames(features) ==
                  rownames(limits))
        keep <- apply(is.finite(features), 2, all)
        features[, keep, drop=FALSE]
      })
      fit <-
        regularized.interval.regression(fmat, set.data$train$limits,
                                        L0=ncol(fmat)*1.5,
                                        calc.grad=calc.grad.list$square,
                                        calc.loss=calc.loss.list$square)
      ## Use these models on the train and validation sets.
      set.data.error <- list()
      for(tv in c("train", "validation")){
        fl <- set.data[[tv]]
        pred.log.lambda <- fit$predict(fl$features)
        too.hi <- fl$limits[,2] < pred.log.lambda
        too.lo <- pred.log.lambda < fl$limits[,1]
        is.error <- too.hi | too.lo
        errors <- colSums(is.error)
        possible.errors <- nrow(is.error)
        set.data.error[[tv]] <- 
        data.table(set=tv, gamma=fit$gamma, errors, possible.errors)
      }#tv
      tv.error <- do.call(rbind, set.data.error)
      tv.error[, percent := errors/possible.errors*100]
      setkey(tv.error, set)
      validation <- tv.error["validation"]
      best.i <- pick.best.index(validation$errors)
      gamma.vals[[validation.fold]] <- validation$gamma[best.i]
      best.dot <- validation[best.i, ]
      ggplot()+
        geom_vline(aes(xintercept=-log10(gamma)), data=best.dot, color="grey")+
        geom_line(aes(-log10(gamma), percent, linetype=set),
                  data=tv.error)+
        geom_point(aes(-log10(gamma), percent),
                   data=best.dot, pch=1)
    }#validation.fold
    mean.gamma <- mean(unlist(gamma.vals))
    ## fit model on combined train/validation set.
    is.tv <- names(best.fl) %in% train.validation
    set.chunks <- 
      list(train.validation=names(best.fl)[!is.tv],
           test=names(best.fl)[is.tv])
    set.data <- list()
    for(tv in c("train.validation", "test")){
      for(data.type in c("features", "limits")){
        data.list <- list()
        for(set.chunk in set.chunks[[tv]]){
          sample.list <- best.fl[[set.chunk]]
          for(sample.id in names(sample.list)){
            fl <- sample.list[[sample.id]]
            data.list[[paste(set.chunk, sample.id)]] <- fl[[data.type]]
          }
        }
        set.data[[tv]][[data.type]] <- do.call(rbind, data.list)
      }
    }
    fmat <- with(set.data$train.validation, {
      stopifnot(rownames(features) ==
                rownames(limits))
      keep <- apply(is.finite(features), 2, all)
      features[, keep, drop=FALSE]
    })
    fit <-
      smooth.interval.regression(fmat, set.data$train.validation$limits,
                                 L0=ncol(fmat)*1.5,
                                 gamma=mean.gamma,
                                 calc.grad=calc.grad.list$square,
                                 calc.loss=calc.loss.list$square)

    names(fit$weights)[fit$weights != 0] #selected features.
    ## predict peaks on the test set.
    pred.log.lambda <- fit$predict(set.data$test$features)
    ## TODO: combine peaks on each chromosome.

    ## TODO: evaluate test error for each chrom.
    
  }#set.i
}#set.name

all.best <- do.call(rbind, all.best.list)

L <- split(all.best, all.best$set.name)
lapply(L, with, table(bases.per.bin, variable))
