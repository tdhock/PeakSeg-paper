works_with_R("3.1.2",
             "tdhock/PeakSegDP@5bcee97f494dcbc01a69e0fe178863564e9985bc")

load("dp.peaks.optimal.RData")
load("dp.peaks.matrices.RData")
load("dp.peaks.features.RData")

dp.peaks.intervals <- list()
for(set.name in names(dp.peaks.matrices)){
  matrices <- dp.peaks.matrices[[set.name]]
  optimal <- dp.peaks.optimal[[set.name]]
  for(chunk.name in names(optimal)){
    chunk.mat <- dp.peaks.features[[chunk.name]]
    optimal.list <- optimal[[chunk.name]]
    error.mat <- matrices[[chunk.name]]$PeakSeg
    chunk.intervals <- NULL
    chunk.features <- NULL
    counts.file <- sprintf("data/%s/counts.RData", chunk.name)
    print(counts.file)
    load(counts.file)
    counts.list <- split(counts, counts$sample.id)
    for(sample.id in names(optimal.list)){
      sample.counts <- counts.list[[sample.id]]
      one.feature <- with(sample.counts, {
        cbind(log.max.coverage=log(max(coverage)),
              log.total.weight=log(sum(chromEnd-chromStart)))
      })
      optimal.df <- optimal.list[[sample.id]]
      param.name <- as.character(optimal.df$model.complexity)
      optimal.df$error <- error.mat[sample.id, param.name]
      if(any(is.na(optimal.df$error))){
        stop("NA error")
      }
      indices <- with(optimal.df, {
        largestContinuousMinimum(error, max.log.lambda-min.log.lambda)
      })
      one.interval <- 
        cbind(min.log.lambda=optimal.df$min.log.lambda[indices$start],
              max.log.lambda=optimal.df$max.log.lambda[indices$end])
      rownames(one.interval) <- rownames(one.feature) <-
        paste(chunk.name, sample.id)
      chunk.intervals <- rbind(chunk.intervals, one.interval)
      chunk.features <- rbind(chunk.features, one.feature)
    }
    all.equal.show <- function(x, y){
      x <- as.numeric(x)
      y <- as.numeric(y)
      stopifnot(all.equal(x, y))
      rbind(x, y)
    }
    ## This verifies that the new feature computation in
    ## dp.peaks.features is the same as the old feature computation in
    ## this script.
    all.equal.show(chunk.mat[, "log.bases"],
                   chunk.features[, "log.total.weight"])
    all.equal.show(chunk.mat[, "log.weighted.quartile.100%"],
                   chunk.features[, "log.max.coverage"])
    dp.peaks.intervals[[set.name]][[chunk.name]] <-
      list(features=chunk.features,
           all.features=chunk.mat,
           intervals=chunk.intervals)
  }
}

save(dp.peaks.intervals, file="dp.peaks.intervals.RData")
