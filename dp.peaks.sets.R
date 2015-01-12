load("dp.peaks.matrices.RData")

sapply(dp.peaks.matrices, length)

set.seed(1)

dp.peaks.sets <- list()
for(set.name in names(dp.peaks.matrices)){
  set <- dp.peaks.matrices[[set.name]]
  chunk.names <- names(set)
  train.sets <- list()
  for(i in 1:3){
    for(j in (i+1):4){
      if(length(chunk.names)==4){
        train.i <- c(i, j)
        train.sets[[length(train.sets)+1]] <- chunk.names[train.i]
      }else{
        train.sets[[length(train.sets)+1]] <-
          sample(chunk.names, length(chunk.names)/2)
      }
    }
  }
  dp.peaks.sets[[set.name]] <- train.sets
}

save(dp.peaks.sets, file="dp.peaks.sets.RData")
