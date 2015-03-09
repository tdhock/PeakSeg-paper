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

min2 <- function(x, y){
  ifelse(x < y, x, y)
}
max2 <- function(x, y){
  ifelse(y < x, x, y)
}
overlapsNext <- function(chromStart, chromEnd){
  c(chromStart[-1] < chromEnd[-length(chromEnd)],
    FALSE)
}  

all.best.list <- list()
for(set.name in names(dp.peaks.sets)){
  train.sets <- dp.peaks.sets[[set.name]]
  experiment <- sub("_.*", "", set.name)
  label.dir <- paste0("~/genomelabels/", set.name)
  RData.files <- Sys.glob(file.path(label.dir, "*.RData"))
  names(RData.files) <- sub(".RData$", "", basename(RData.files))
  RData.by.sample <- list()
  e <- new.env()
  fl.list <- list()
  error.list <- list()
  chrom.list <- list()
  region.list <- list()
  sample.ids <- names(RData.files)
  for(sample.id in sample.ids){
    RData.file <- RData.files[[sample.id]]
    objs <- load(RData.file, e)
    RData.by.sample[[sample.id]] <- as.list(e)
    error.list[[sample.id]] <- data.table(sample.id, e$errors)
    region.list[[sample.id]] <- data.table(sample.id, e$regions)
    to.check <- e$errors %>%
      group_by(bases.per.bin) %>%
      summarise(total.weight=sum(total.weight))
    exp.weight <- rep(to.check$total.weight[1], nrow(to.check))
    stopifnot(all.equal(exp.weight,
                        to.check$total.weight))
    for(bases.per.bin in names(e$features.limits)){
      fl <- e$features.limits[[bases.per.bin]]
      stopifnot(!is.null(fl$problems))
      chunk.list <- fl$limits
      for(chunk.id in names(chunk.list)){
        limits <- chunk.list[[chunk.id]]
        chrom <- sub(":.*", "", rownames(limits)[1])
        all.problems <- fl$problems[[chrom]]
        setkey(all.problems, problem.name)
        problems <- all.problems[rownames(limits)]
        modelSelection <- fl$modelSelection[rownames(limits)]
        peaks <- fl$peaks[rownames(limits)]
        chunk.name <- paste0(set.name, "/", chunk.id)
        chrom.list[[chrom]][[chunk.name]] <- chunk.name
        features <- fl$features[rownames(limits), , drop=FALSE ]
        ## TODO: immediately filter out non-finite features.
        
        ## features and limits are matrices[problems, ]
        ## rownames(features) <- rownames(limits) <-
        ##   paste(chunk.name, bases.per.bin, rownames(limits))
        ## The problem index is INSIDE on the matrix rows.
        if(!is.null(fl.list[[bases.per.bin]][[chunk.name]][[sample.id]]))
          stop("overwriting some data, this should never happen")
        fl.list[[bases.per.bin]][[chunk.name]][[sample.id]] <-
          list(features=features, limits=limits,
               modelSelection=modelSelection,
               peaks=peaks,
               problems=problems)
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
    for(chrom in names(chrom.list)){
      chrom.chunks <- names(chrom.list[[chrom]])
      for(sample.id in sample.ids){
        sample.regions <- region.list[[sample.id]]
        setkey(sample.regions, chrom)
        chrom.regions <- sample.regions[chrom]
        chrom.peak.list <- list()
        problem.list <- list()
        for(chunk.name in chrom.chunks){
          chunk.info <- best.fl[[chunk.name]][[sample.id]]
          pred.log.lambda <- fit$predict(chunk.info$features)
          for(problem.name in names(chunk.info$modelSelection)){
            l <- pred.log.lambda[problem.name, ]
            exact <- chunk.info$modelSelection[[problem.name]]
            interval <- subset(exact, min.log.lambda < l & l < max.log.lambda)
            stopifnot(nrow(interval) == 1)
            peaks.str <- paste(interval$peaks)
            peaks <- chunk.info$peaks[[problem.name]][[peaks.str]]
            if(nrow(peaks)){
              problem <- chunk.info$problems[problem.name]
              problem.list[[problem.name]] <- problem
              peaks$problemStart <- problem$chromStart
              peaks$problemPeakStart <- problem$peakStart
              peaks$problemPeakEnd <- problem$peakEnd
              peaks$problemEnd <- problem$chromEnd
              peaks$problem.name <- problem.name
              chrom.peak.list[[paste(chunk.name, problem.name)]] <- peaks
            }
          }#problem.name
        }#chunk.name
        problem.rects <- do.call(rbind, problem.list)
        setkey(problem.rects, chromStart, chromEnd)
        chrom.peaks <- do.call(rbind, chrom.peak.list) %>%
          arrange(chromStart, chromEnd) %>%
          mutate(is.before=chromEnd < problemPeakStart,
                 is.after=problemPeakEnd < chromStart,
                 is.overlap=(!is.before)&(!is.after),
                 status=ifelse(is.overlap, "keep", "discard"))
        chrom.peaks$peak.i <- 1:nrow(chrom.peaks)
        filtered.peaks <- chrom.peaks %>%
          filter(is.overlap) %>%
          mutate(overlaps.next.peak=overlapsNext(chromStart, chromEnd),
                 status=ifelse(overlaps.next.peak, "merge", "keep"),
                 overlaps.prev.peak=c(FALSE,
                   overlaps.next.peak[-length(overlaps.next.peak)]))
        p.next <- filtered.peaks[filtered.peaks$overlaps.next.peak, ]
        p.prev <- filtered.peaks[filtered.peaks$overlaps.prev.peak, ]
        edited.peaks <- filtered.peaks %>%
          filter(!overlaps.next.peak) %>%
          mutate(status=ifelse(overlaps.prev.peak, "edited", "same"))
        edited.peaks$chromStart[edited.peaks$overlaps.prev.peak] <-
          min2(p.next$chromStart, p.prev$chromStart)
        edited.peaks$chromEnd[edited.peaks$overlaps.prev.peak] <-
          max2(p.next$chromEnd, p.prev$chromEnd)
        overlap.peaks <- filtered.peaks %>%
          filter(overlaps.next.peak) %>%
          mutate(overlap.i=seq_along(overlaps.next.peak))
        overlap.edited <- edited.peaks %>%
          filter(overlaps.prev.peak)
        if(nrow(overlap.peaks) > 0){
          zoom.problem.list <- list()
          zoom.peak.list <- list()
          for(overlap.i in 1:nrow(overlap.peaks)){
            overlap.peak <- data.table(overlap.peaks[overlap.i, ]) %>%
              mutate(mid=(chromStart+chromEnd)/2,
                     zoom=overlap.i)
            edited.peak <- overlap.edited[overlap.i, ] %>%
              mutate(problem.name="edited")
            setkey(overlap.peak, chromStart, chromEnd)
            overlap.rects <- foverlaps(problem.rects, overlap.peak, nomatch=0L)
            zoom.problem.list[[overlap.i]] <- overlap.rects
            zoom.chromEnd <- max(overlap.rects$chromEnd)
            zoom.chromStart <- min(overlap.rects$chromStart)
            relevant.peaks <- rbind(filtered.peaks, edited.peak) %>%
              filter(chromStart < zoom.chromEnd,
                     zoom.chromStart < chromEnd) %>%
              mutate(mid=overlap.peak$mid,
                     zoom=overlap.peak$zoom)
            zoom.peak.list[[overlap.i]] <- relevant.peaks
          }
          zoom.problems <- do.call(rbind, zoom.problem.list)
          zoom.peaks <- do.call(rbind, zoom.peak.list)
        ggplot()+
          geom_segment(aes((i.chromStart-mid)/1e3, i.problem.name,
                           xend=(i.chromEnd-mid)/1e3, yend=i.problem.name),
                        data=zoom.problems, size=5, alpha=1/10)+
          geom_segment(aes((peakStart-mid)/1e3, i.problem.name,
                           xend=(peakEnd-mid)/1e3, yend=i.problem.name),
                        data=zoom.problems, size=5, alpha=1/10)+
          geom_segment(aes((chromStart-mid)/1e3, problem.name,
                           xend=(chromEnd-mid)/1e3, yend=problem.name,
                           color=status),
                       data=zoom.peaks, size=2)+
          theme_bw()+
          theme(panel.margin=grid::unit(0, "cm"))+
          facet_grid(zoom ~ ., scales="free")
        }
        if(FALSE){
        ggplot()+
          geom_segment(aes(chromStart/1e3, peak.i,
                           xend=chromEnd/1e3, yend=peak.i,
                           color=status),
                       data=chrom.peaks)
        ggplot()+
          geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3),
                        data=problem.rects, alpha=1/10)+
          geom_tallrect(aes(xmin=peakStart/1e3, xmax=peakEnd/1e3),
                        data=problem.rects, alpha=1/10)+
          geom_segment(aes(chromStart/1e3, peak.i,
                           xend=chromEnd/1e3, yend=peak.i,
                           color=status),
                       data=chrom.peaks)+
          theme_bw()+
          theme(panel.margin=grid::unit(0, "cm"))+
          facet_grid(problem.name ~ .)
        ggplot()+
          geom_segment(aes(chromStart/1e3, problem.name,
                           xend=chromEnd/1e3, yend=problem.name),
                        data=problem.rects, size=5, alpha=1/10)+
          geom_segment(aes(peakStart/1e3, problem.name,
                           xend=peakEnd/1e3, yend=problem.name),
                        data=problem.rects, size=5, alpha=1/10)+
          geom_segment(aes(chromStart/1e3, problem.name,
                           xend=chromEnd/1e3, yend=problem.name,
                           color=status),
                       data=chrom.peaks, size=2)
        ggplot()+
          geom_segment(aes(chromStart/1e3, problem.name,
                           xend=chromEnd/1e3, yend=problem.name),
                        data=problem.rects, size=5, alpha=1/10)+
          geom_segment(aes(peakStart/1e3, problem.name,
                           xend=peakEnd/1e3, yend=problem.name),
                        data=problem.rects, size=5, alpha=1/10)+
          geom_segment(aes(chromStart/1e3, problem.name,
                           xend=chromEnd/1e3, yend=problem.name,
                           color=status),
                       data=filtered.peaks, size=2)
        }
        sorted.peaks <- edited.peaks %>%
          arrange(chromStart, chromEnd) %>%
          mutate(overlaps.next.peak=overlapsNext(chromStart, chromEnd))
        stopifnot(all(!sorted.peaks$overlaps.next.peak))
        ##     for every chunk
        ##       evaluate peaks
      }#sample.id
    }#chrom
    stop(1)

    pred.log.lambda <- fit$predict(set.data$test$features)
    ## TODO: combine peaks on each chromosome.

    ## TODO: evaluate test error for each chrom.
    
  }#set.i
}#set.name

all.best <- do.call(rbind, all.best.list)

L <- split(all.best, all.best$set.name)
lapply(L, with, table(bases.per.bin, variable))
