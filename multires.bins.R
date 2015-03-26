works_with_R("3.1.3",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             reshape2="1.2.2",
             data.table="1.9.4",
             "tdhock/PeakSegDP@af1e3af69e886e9cbf45e03e283083c88fb9fa43",
             dplyr="0.4.0")

load("dp.peaks.sets.RData")
load("dp.peaks.error.RData")
load("oracle.regularized.RData")

ref.err <- subset(oracle.regularized$error, model.name==model.name[1])

ref.regions <- ref.err %>%
  group_by(set.name, set.i) %>%
  summarise(regions=sum(regions)) %>%
  mutate(testSet=paste(set.name, "split", set.i))
setkey(ref.regions, testSet)

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

overlapsNext <- function(chromStart, chromEnd){
  next.chromStart <- chromStart[-1]
  prev.chromEnd <- chromEnd[-length(chromEnd)]
  cummax.chromEnd <- cummax(prev.chromEnd)
  c(next.chromStart < cummax.chromEnd,
    FALSE)
}  

all.best.list <- list()
test.error.list <- list()
test.peak.list <- list()
for(set.name in names(dp.peaks.sets)){
  all.chunks <- grep(set.name, names(dp.peaks.error), value=TRUE)
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
  bad.list <- list()
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
        ## immediately filter out non-finite features.
        is.bad <- apply(!is.finite(features), 2, any)
        bad.list[[bases.per.bin]][[paste(sample.id, chunk.id)]] <-
          colnames(features)[is.bad]
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
    chunks.used <- rep(NA, length(all.chunks))
    names(chunks.used) <- all.chunks
    chunks.used[train.validation] <- "train.validation"

    ## Check to make sure there are the right number of test regions
    ## in this dats set.
    expected.regions <- ref.regions[testSet]$regions
    set.regions <- do.call(rbind, region.list) %>%
      mutate(chunk.name=paste0(set.name, "/", chunk.id))
    test.regions <- set.regions %>%
      filter(! chunk.name %in% train.validation)
    stopifnot(nrow(test.regions) == expected.regions)

    set.test.error.list <- list()

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
    best.bad <- unique(do.call(c, bad.list[[best.bases.per.bin]]))

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
        keep <- ! colnames(features) %in% best.bad
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
      list(train.validation=names(best.fl)[is.tv],
           test=names(best.fl)[!is.tv])
    chunks.used[set.chunks$test] <- "test"
    stopifnot(!is.na(chunks.used))
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
      keep <- ! colnames(features) %in% best.bad
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
    test.chunk.list <- list()
    for(chrom in names(chrom.list)){
      chrom.chunks <- names(chrom.list[[chrom]])
      is.test.name <- chrom.chunks %in% set.chunks$test
      test.chunks <- chrom.chunks[is.test.name]
      test.chunk.list[[chrom]] <- test.chunks
      if(length(test.chunks))for(sample.id in sample.ids){
        sample.regions <- region.list[[sample.id]] %>%
          mutate(chunk.name=paste0(set.name, "/", chunk.id))
        setkey(sample.regions, chunk.name)
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
            problem <- chunk.info$problems[problem.name]
            problem.list[[problem.name]] <- problem

            peaks.str <- paste(interval$peaks)
            peaks <- chunk.info$peaks[[problem.name]][[peaks.str]]
            if(nrow(peaks)){
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
        if(length(chrom.peak.list) == 0){
          sorted.peaks <- Peaks()
        }else{
          chrom.peaks <- do.call(rbind, chrom.peak.list) %>%
            arrange(chromStart, chromEnd) %>%
            mutate(is.before=chromEnd < problemPeakStart,
                 is.after=problemPeakEnd < chromStart,
                 is.overlap=(!is.before)&(!is.after),
                 status=ifelse(is.overlap, "keep", "discard"))
          chrom.peaks$peak.i <- 1:nrow(chrom.peaks)
          filtered.peaks <- chrom.peaks %>%
            filter(is.overlap) %>%
            arrange(chromStart, chromEnd) %>%
            mutate(overlaps.next.peak=overlapsNext(chromStart, chromEnd),
                 overlaps.prev.peak=c(FALSE,
                   overlaps.next.peak[-length(overlaps.next.peak)]),
                 status=ifelse(overlaps.next.peak|overlaps.prev.peak,
                   "merge", "keep"),
                 diff=c(ifelse(overlaps.next.peak[1], 1, 0),
                   diff(overlaps.next.peak)))
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
          cluster.starts <- which(filtered.peaks$diff == 1)
          cluster.ends <- which(filtered.peaks$diff == -1)
          stopifnot(length(cluster.starts) == length(cluster.ends))
          merged.peaks <- if(length(cluster.starts) > 0){
            zoom.problem.list <- list()
            zoom.peak.list <- list()
            for(cluster.i in seq_along(cluster.starts)){
              cluster.start <- cluster.starts[[cluster.i]]
              cluster.end <- cluster.ends[[cluster.i]]
              overlap.peak <- data.table(filtered.peaks[1,]) %>%
                mutate(chromStart=filtered.peaks$chromStart[cluster.start],
                       chromEnd=filtered.peaks$chromEnd[cluster.end],
                       problem.name="edited")
              problem.names <-
                filtered.peaks$problem.name[c(cluster.start, cluster.end)]
              mid <- filtered.peaks$problemPeakEnd[cluster.start]
              overlap.rects <- problem.rects %>%
                filter(problem.name %in% problem.names)
              zoom.problem.list[[cluster.i]] <-
                data.table(overlap.rects, cluster.i)
              problem.peaks <- filtered.peaks %>%
                filter(problem.name %in% problem.names)
              relevant.peaks <-
                data.table(rbind(problem.peaks, overlap.peak), mid, cluster.i)
              zoom.peak.list[[cluster.i]] <- relevant.peaks
            }
            zoom.problems <- do.call(rbind, zoom.problem.list)
            zoom.peaks <- do.call(rbind, zoom.peak.list)
            ggplot()+
              geom_segment(aes((chromStart-mid)/1e3, problem.name,
                           xend=(chromEnd-mid)/1e3, yend=problem.name),
                       data=zoom.problems, size=5, alpha=1/10)+
              geom_segment(aes((peakStart-mid)/1e3, problem.name,
                           xend=(peakEnd-mid)/1e3, yend=problem.name),
                        data=zoom.problems, size=5, alpha=1/10)+
              geom_segment(aes((chromStart-mid)/1e3, problem.name,
                           xend=(chromEnd-mid)/1e3, yend=problem.name,
                           color=status),
                       data=zoom.peaks, size=2)+
              theme_bw()+
              theme(panel.margin=grid::unit(0, "cm"))+
              facet_grid(cluster.i ~ ., scales="free")
            zoom.peaks %>%
              filter(problem.name == "edited")
          }
          edited.peaks <-
            rbind(if(!is.null(merged.peaks))merged.peaks %>%
                  select(chromStart, chromEnd),
                  filtered.peaks %>%
                  filter(status=="keep") %>%
                  select(chromStart, chromEnd))
          sorted.peaks <- edited.peaks %>%
            arrange(chromStart, chromEnd) %>%
            mutate(overlaps.next.peak=overlapsNext(chromStart, chromEnd))
          stopifnot(all(!sorted.peaks$overlaps.next.peak))
          ## TDH 26 Mar 2015. try a simpler strategy involving clusterPeaks.
          all.peaks <- do.call(rbind, chrom.peak.list) %>%
            arrange(chromStart, chromEnd) %>%
            mutate(is.before=chromEnd < problemPeakStart,
                 is.after=problemPeakEnd < chromStart,
                   not.in.region=is.before|is.after,
                 is.overlap=!not.in.region,
                 status=ifelse(is.overlap, "keep", "discard"))
          some.peaks <- all.peaks %>%
            filter(is.overlap)
          clustered.peaks <- clusterPeaks(some.peaks)
          table(clustered.peaks$cluster)
          clustered.peak.list <- split(clustered.peaks, clustered.peaks$cluster)
          reduced.peak.list <- list()
          for(cluster.str in names(clustered.peak.list)){
            peaks <- clustered.peak.list[[cluster.str]]
            reduced.peak.list[[cluster.str]] <- 
              data.frame(chromStart=peaks$chromStart[1],
                         chromEnd=peaks$chromEnd[nrow(peaks)])
          }
          sorted.peaks <- do.call(rbind, reduced.peak.list) %>%
            arrange(chromStart, chromEnd) %>%
            mutate(overlaps.next.peak=overlapsNext(chromStart, chromEnd))
          stopifnot(all(!sorted.peaks$overlaps.next.peak))
        }#if/else
        test.peak.list[[testSet]][[chrom]][[sample.id]] <- sorted.peaks
        for(chunk.name in test.chunks){
          sid <- sample.id
          chunk.ref <- do.call(rbind, dp.peaks.error[[chunk.name]]) %>%
            filter(param.name == param.name[1],
                   sample.id == sid)
          test.regions <- sample.regions[chunk.name]
          stopifnot(nrow(chunk.ref) == nrow(test.regions))
          error <- PeakErrorChrom(sorted.peaks, test.regions)
          set.test.error.list[[paste(sample.id, chunk.name)]] <-
          data.table(sample.id, chunk.name, chrom, set.name, testSet, set.i,
                     error)
        }
      }#sample.id
    }#chrom
    set.test.error <- do.call(rbind, set.test.error.list)
    stopifnot(nrow(set.test.error) == expected.regions)
    test.error.list[[testSet]] <- set.test.error
  }#set.i
}#set.name
test.error <- do.call(rbind, test.error.list)
all.best <- do.call(rbind, all.best.list)

L <- split(all.best, all.best$set.name)
lapply(L, with, table(bases.per.bin, variable))

multires.bins <-
  list(test.error=test.error,
       hyperparameters=all.best,
       test.peaks=test.peak.list)

save(multires.bins, file="multires.bins.RData")
