works_with_R("3.1.1", dplyr="0.2")

load("dp.peaks.sets.RData")
load("dp.peaks.matrices.RData")

source("pick.best.index.R")

dp.peaks.baseline <- NULL
dp.peaks.baseline.roc <- list()
dp.peaks.baseline.chosen <- list()
for(set.name in names(dp.peaks.sets)){
  train.sets <- dp.peaks.sets[[set.name]]
  error.list <- dp.peaks.matrices[[set.name]]
  fp.list <- dp.peaks.matrices.fp[[set.name]]
  tp.list <- dp.peaks.matrices.tp[[set.name]]
  for(set.i in seq_along(train.sets)){
    train.chunks <- train.sets[[set.i]]
    test.chunks <- names(error.list)[!names(error.list) %in% train.chunks]
    test.possible.tp <- 0
    test.possible.fp <- 0
    for(test.chunk in test.chunks){
      test.possible.tp <-
        test.possible.tp + sum(tp.list[[test.chunk]]$possible.tp)
      test.possible.fp <-
        test.possible.fp + sum(fp.list[[test.chunk]]$possible.fp)
    }
    for(algorithm in c("hmcan.broad.trained", "macs.trained")){
      train.mat <- NULL
      for(train.chunk in train.chunks){
        train.mat <- rbind(train.mat, error.list[[train.chunk]][[algorithm]])
      }
      error.curve <- colSums(train.mat)
      picked <- pick.best.index(error.curve)
      test.mat <- NULL
      test.regions <- NULL
      test.fp <- NULL
      test.tp <- NULL
      for(test.chunk in test.chunks){
        test.fp <- rbind(test.fp, fp.list[[test.chunk]][[algorithm]])
        test.tp <- rbind(test.tp, tp.list[[test.chunk]][[algorithm]])
        test.mat <- rbind(test.mat, error.list[[test.chunk]][[algorithm]])
        test.regions <- c(test.regions, error.list[[test.chunk]]$regions)
      }
      param.name <- colnames(test.tp)
      tp.fp <-
        data.frame(set.name, set.i,
                   algorithm, param.name,
                   tp=colSums(test.tp),
                   fp=colSums(test.fp)) %>%
        mutate(TPR=tp/test.possible.tp,
               FPR=fp/test.possible.fp)
      tp.fp.chosen <- tp.fp[picked, ]
      algo.code <- paste(set.name, set.i, algorithm)
      dp.peaks.baseline.roc[[algo.code]] <- tp.fp
      dp.peaks.baseline.chosen[[algo.code]] <- tp.fp.chosen
      errors <- sum(test.mat[,picked])
      regions <- sum(test.regions)
      dp.peaks.baseline <- rbind(dp.peaks.baseline, {
        data.frame(set.name, set.i, algorithm,
                   param.name=param.name[picked],
                   errors, regions)
      })
    }
  }
}

save(dp.peaks.baseline,
     dp.peaks.baseline.roc,
     dp.peaks.baseline.chosen,
     file="dp.peaks.baseline.RData")
