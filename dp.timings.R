works_with_R("3.1.1",
             "tdhock/PeakSegDP@5bcee97f494dcbc01a69e0fe178863564e9985bc",
             data.table="1.9.3")

library(parallel)

prefix <- "http://cbio.ensmp.fr/~thocking/chip-seq-chunk-db/"
u <- paste0(prefix, "RData.count.signal.txt")
count.files <- read.table(u, header=TRUE)
maxPeaks <- 9L
dp.timings <- NULL

for(file.i in 1:nrow(count.files)){
  f <- count.files$file[[file.i]]
  chunk.id.dir <- dirname(paste(f))
  chunk.id <- basename(chunk.id.dir)
  set.dir <- dirname(chunk.id.dir)
  set.name <- basename(set.dir)
  u <- paste0(prefix, f)
  local.f <- file.path("data", f)
  print(local.f)
  if(!file.exists(local.f)){
    local.dir <- dirname(local.f)
    dir.create(local.dir, showWarnings=FALSE, recursive=TRUE)
    download.file(u, local.f)
  }
  time.f <- sub("counts", "dp.timing", local.f)
  model.f <- sub("counts", "dp.model", local.f)
  ## run DP based on if we have a model computed or not.
  if(file.exists(time.f)){
    load(time.f)
  }else{
    load(local.f) #counts
    counts$bases <- with(counts, chromEnd-chromStart)
    counts$count <- counts$coverage
    sample.list <- split(counts, counts$sample.id, drop=TRUE)
    sample.ids <- names(sample.list)
    sample.results <- mclapply(seq_along(sample.ids), function(sample.i){
      sample.id <- sample.ids[[sample.i]]
      compressed <- sample.list[[sample.id]]
      bases <- sum(compressed$bases)
      n.data <- nrow(compressed)
      cat(sprintf("%4d / %4d chunks %4d / %4d sample %s %d bases %d data\n",
                  file.i, nrow(count.files),
                  sample.i, length(sample.list),
                  sample.id,
                  bases, n.data))
      result <- list()
      seconds <- system.time({
        result$model <- PeakSegDP(compressed, maxPeaks=maxPeaks)
      })[["elapsed"]]
      result$timing <- 
        data.frame(set.name, chunk.id, sample.id,
                   seconds, data=n.data, bases)
      result
    })
    names(sample.results) <- sample.ids
    dp.model <- lapply(sample.results, "[[", "model")
    timing.list <- lapply(sample.results, "[[", "timing")
    dp.timing <- do.call(rbind, timing.list)
    ## ggplot()+
    ##   geom_point(aes(data, seconds),
    ##              data=dp.timing, pch=1)+
    ##   scale_x_log10()+
    ##   scale_y_log10()
    save(dp.model, file=model.f)
    save(dp.timing, file=time.f)
  }
  dp.timings <- rbind(dp.timings, dp.timing)
}

save(dp.timings, file="dp.timings.RData")
