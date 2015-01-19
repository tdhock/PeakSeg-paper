load("unsupervised.RData")
load("dp.peaks.matrices.RData")
load("dp.peaks.sets.RData")

default.params <-
  c(macs.trained="1.30103",
    hmcan.broad.trained="2.30258509299405")

elist <- list()
for(set.name in names(dp.peaks.sets)){
  train.sets <- dp.peaks.sets[[set.name]]
  chunk.list <- dp.peaks.matrices[[set.name]]
  for(set.i in seq_along(train.sets)){
    testSet <- paste(set.name, "split", set.i)
    train.chunks <- train.sets[[set.i]]
    cat(sprintf("%d / %d %s\n", set.i, length(train.sets), set.name))
    test.chunks <- names(chunk.list)[! names(chunk.list) %in% train.chunks]
    for(test.chunk in test.chunks){
      test.info <- chunk.list[[test.chunk]]
      seg.mat <- unsupervised[[test.chunk]]
      err.mat <- test.info$PeakSeg
      regions <- test.info$regions
      stopifnot(rownames(seg.mat) == rownames(err.mat))
      for(algorithm in colnames(seg.mat)){
        segs <- seg.mat[, algorithm]
        peaks <- (segs-1)/2
        param.name <- as.character(peaks)
        sample.id <- names(segs)
        i.mat <- cbind(sample.id, param.name)
        errors <- err.mat[i.mat]
        elist[[paste(set.name, set.i, test.chunk, algorithm)]] <- 
          data.table(set.name, set.i, testSet, test.chunk, algorithm,
                     param.name, sample.id, errors, regions)
      }
      for(algorithm in names(default.params)){
        err.mat <- test.info[[algorithm]]
        param.name <- default.params[[algorithm]]
        errors <- err.mat[, param.name]
        sample.id <- names(errors)
        algorithm <- sub("trained", "default", algorithm)
        elist[[paste(set.name, set.i, test.chunk, algorithm)]] <- 
          data.table(set.name, set.i, testSet, test.chunk, algorithm,
                     param.name, sample.id, errors, regions)
      }
    }
  }
}

unsupervised.error <- do.call(rbind, elist)

save(unsupervised.error, file="unsupervised.error.RData")
