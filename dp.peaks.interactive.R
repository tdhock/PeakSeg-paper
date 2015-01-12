works_with_R("3.1.1", dplyr="0.2", ggplot2="1.0")

load("dp.peaks.regression.dots.RData")
load("dp.peaks.regression.RData")
load("dp.peaks.baseline.RData")
load("dp.peaks.RData")
load("dp.peaks.error.RData")

profile.list <- list()

prefix <- "http://cbio.ensmp.fr/~thocking/chip-seq-chunk-db"

dp.peaks.regression %.%
  group_by(set.name, set.i) %.%
  summarise(samples=n(),
            test.chunks=length(unique(test.chunk)))

data.list <-
  list(splitErrors=dp.peaks.regression.dots,
       baselineParams=dp.peaks.baseline,
       regressionParams=dp.peaks.regression)

by.testSet <- list()
for(data.name in names(data.list)){
  data.list[[data.name]]$testSet <-
    with(data.list[[data.name]], {
      paste(set.name, "split", set.i)
    })
  by.testSet[[data.name]] <-
    split(data.list[[data.name]],
          data.list[[data.name]]$testSet)
}

## First download all the data that we need to show the test chunks.
test.chunks <-
  as.character(unique(data.list$regressionParams$test.chunk))
baseline.regions <- list()
baseline.peaks <- list()
## Note that there are 125 test.chunks (out of 129 chunks total in the
## benchmark database,
## http://cbio.ensmp.fr/~thocking/chip-seq-chunk-db/file-rows.html)
## since we randomly chose 1/2 train, 1/2 test, 6 times in
## dp.peaks.sets, there are 4 chunks which appear in no test set.
for(test.chunk.i in seq_along(test.chunks)){
  test.chunk <- test.chunks[[test.chunk.i]]
  cat(sprintf("%4d / %4d %s\n", test.chunk.i, length(test.chunks), test.chunk))
  for(algorithm in c("hmcan.broad.trained", "macs.trained")){
    u <- url(sprintf("%s/%s/error/%s.RData", prefix, test.chunk, algorithm))
    load(u)
    close(u)
    baseline.regions[[test.chunk]][[algorithm]] <-
      split(error, error$param.name)
    u <- url(sprintf("%s/%s/peaks/%s.RData", prefix, test.chunk, algorithm))
    load(u)
    close(u)
    baseline.peaks[[test.chunk]][[algorithm]] <- peaks
  }
}

testSets <- names(by.testSet[[1]])

load(url(sprintf("%s/chunks.RData", prefix)))

region.list <- lapply(dp.peaks.error, function(L){
  df <- do.call(rbind, L)
  lapply(split(df, df$sample.id), function(df) split(df, df$param.name))
})

dp.peaks.interactive <- data.list
dp.peaks.interactive$peaks <- dp.peaks.interactive$regionErrors <- NULL
dp.peaks.interactive$chunks <- chunks %.%
  mutate(test.chunk=paste0(set.name, "/", chunk.id),
         bases=expandEnd-expandStart) %.%
  filter(test.chunk %in% test.chunks)
rownames(dp.peaks.interactive$chunks) <- dp.peaks.interactive$chunks$test.chunk
all.peak.list <- list()
all.region.list <- list()
for(testSet.i in seq_along(testSets)){
  testSet <- testSets[[testSet.i]]
  cat(sprintf("%4d / %4d %s\n", testSet.i, length(testSets), testSet))
  baselines <- by.testSet$baselineParams[[testSet]]
  all.test.samples <- by.testSet$regressionParams[[testSet]]
  chunk.list <- split(all.test.samples, all.test.samples$test.chunk, drop=TRUE)
  for(test.chunk in names(chunk.list)){
    chunk <- dp.peaks.interactive$chunks[test.chunk, ]
    normalize <- function(b)(b-chunk$expandStart)/chunk$bases
    normDF <- function(...){
      data.frame(...) %.%
        mutate(chromStartNorm=normalize(chromStart),
               chromEndNorm=normalize(chromEnd))
    }
    peak.meta <- data.frame(testSet, test.chunk)
    test.samples <- chunk.list[[test.chunk]]
    test.sample.ids <- as.character(test.samples$sample.id)
    counts.file <- sprintf("data/%s/counts.RData", test.chunk)
    load(counts.file)
    counts.list <- split(counts, counts$sample.id)
    for(sample.i in 1:nrow(test.samples)){
      s <- test.samples[sample.i, ]
      sample.id <- as.character(s$sample.id)
      sample.coverage <- counts.list[[sample.id]]
      cov.cols <- c("chromStart", "chromEnd", "coverage")
      profile.list[[paste(test.chunk, sample.id)]] <-
        data.frame(test.chunk, sample.id, sample.coverage[, cov.cols])
      peaks.by.sample <- dp.peaks[[test.chunk]]
      dp.models <- peaks.by.sample[[sample.id]]
      param.name <- as.character(s$peaks)
      dp.param.peaks <- dp.models[[param.name]]
      this.id <- paste(testSet, test.chunk, sample.id, "PeakSeg")
      if(!is.null(dp.param.peaks) && nrow(dp.param.peaks)){
        all.peak.list[[this.id]] <- 
          normDF(peak.meta, algorithm="PeakSeg", sample.id,
                 dp.param.peaks[, c("chromStart", "chromEnd")])
      }
      selected.error <- region.list[[test.chunk]][[sample.id]][[param.name]]
      all.region.list[[this.id]] <- 
        normDF(peak.meta, algorithm="PeakSeg", selected.error)
    }
    for(algorithm.i in 1:nrow(baselines)){
      baseline <- baselines[algorithm.i, ]
      algorithm <- as.character(baseline$algorithm)
      param.name <- as.character(baseline$param.name)
      error.by.param <- baseline.regions[[test.chunk]][[algorithm]]
      selected.error <-
        subset(error.by.param[[param.name]], sample.id %in% test.sample.ids)
      this.id <- paste(testSet, test.chunk, sample.id, algorithm)
      all.region.list[[this.id]] <-
        normDF(peak.meta, algorithm, selected.error)
      peaks <- baseline.peaks[[test.chunk]][[algorithm]]
      param.peaks <- peaks[[param.name]]
      if(nrow(param.peaks)){
        selected.peaks <- subset(param.peaks, sample.id %in% test.sample.ids)
        all.peak.list[[this.id]] <- 
          normDF(peak.meta, algorithm, selected.peaks)
      }
    }#algorithm
  }#test chunk
}

dp.peaks.interactive$peaks <- do.call(rbind, all.peak.list)
dp.peaks.interactive$regionErrors <- do.call(rbind, all.region.list)

## We will only display profiles in a plot that is 1000 pixels wide,
## so there is no need to get 254,451 data points. 
range(sapply(profile.list, nrow))

n.pixels <- 1300
profile.out.list <- list()
profile.info.list <- list()
for(coverage.file.i in seq_along(profile.list)){
  cat(sprintf("%4d / %4d profiles\n", coverage.file.i, length(profile.list)))
  profile <- profile.list[[coverage.file.i]]
  test.chunk <- as.character(profile$test.chunk[1])
  small.profile <- if(nrow(profile) > n.pixels){
    base <- with(profile, {
      as.integer(seq(chromStart[1], chromStart[length(chromStart)], l=n.pixels))
    })
    coverage <- with(profile, approx(chromStart, coverage, base))$y
    data.frame(profile[1, 1:2], base, coverage, row.names=NULL)
  }else{
    with(profile[1, ], {
      data.frame(test.chunk, sample.id, base=chromStart, coverage)
    })
  }
  chunk <- dp.peaks.interactive$chunks[test.chunk, ]
  profile.info <- chunk %.%
    mutate(max.coverage=max(profile$coverage),
           test.chunk=test.chunk,
           sample.id=profile$sample.id[1])
  profile.info.list[[coverage.file.i]] <- profile.info
  norm.profile <- small.profile %.%
    mutate(baseNorm=normalize(base),
           coverageNorm=coverage/profile.info$max.coverage)
  ggplot()+
    geom_step(aes(baseNorm, coverageNorm),
              data=norm.profile)
  ggplot()+
    geom_step(aes(chromStart/1e3, coverage, color=what),
              data=data.frame(profile, what="original"))+
    geom_step(aes(base/1e3, coverage, color=what),
              data=data.frame(small.profile, what="small"))
  profile.out.list[[coverage.file.i]] <- norm.profile
}
dp.peaks.interactive$profiles <- do.call(rbind, profile.out.list)
dp.peaks.interactive$profile.info <- do.call(rbind, profile.info.list)

range(dp.peaks.interactive$profiles$baseNorm)
range(dp.peaks.interactive$peaks$chromStartNorm)
range(dp.peaks.interactive$peaks$chromEndNorm)
range(dp.peaks.interactive$regionErrors$chromStartNorm)
range(dp.peaks.interactive$regionErrors$chromEndNorm)

save(dp.peaks.interactive, file="dp.peaks.interactive.RData")
