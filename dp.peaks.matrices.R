works_with_R("3.1.1", dplyr="0.2")

load("dp.peaks.error.RData")

prefix <- "http://cbio.ensmp.fr/~thocking/chip-seq-chunk-db"

set.names <- unique(sub("/.*", "", names(dp.peaks.error)))

dp.peaks.matrices <- list()
dp.peaks.matrices.fp <- list()
dp.peaks.matrices.tp <- list()
for(set.name in set.names){
  chunk.names <- grep(set.name, names(dp.peaks.error), value=TRUE)
  for(chunk.name in chunk.names){
    print(chunk.name)
    chunk.list <- dp.peaks.error[[chunk.name]]
    chunk.df <- do.call(rbind, chunk.list)
    long <- chunk.df %.%
      mutate(param.num=as.numeric(as.character(param.name))) %.%
      group_by(sample.id, param.num) %.%
      summarise(errors=sum(fp+fn),
                fp=sum(fp),
                tp=sum(tp),
                possible.fp=sum(possible.fp),
                possible.tp=sum(possible.tp),
                regions=n())
    long.list <- split(long, long$sample.id, drop=TRUE)
    err.mat <- fp.mat <- tp.mat <-
      matrix(NA, length(long.list), 10,
             dimnames=list(sample.id=names(long.list),
               param.name=0:9))
    for(row.i in seq_along(long.list)){
      sample.df <- long.list[[row.i]]
      param.name <- as.character(sample.df$param.num)
      err.mat[row.i, param.name] <- sample.df$errors
      fp.mat[row.i, param.name] <- sample.df$fp
      tp.mat[row.i, param.name] <- sample.df$tp
    }
    err.list <-
      list(PeakSeg=err.mat,
           regions=sapply(long.list, function(x)x$regions[[1]]))
    fp.list <-
      list(PeakSeg=fp.mat,
           possible.fp=sapply(long.list, function(x)x$possible.fp[[1]]))
    tp.list <-
      list(PeakSeg=tp.mat,
           possible.tp=sapply(long.list, function(x)x$possible.tp[[1]]))
    for(algorithm in c("macs.trained", "hmcan.broad.trained")){
      u <- url(sprintf("%s/%s/error/%s.RData", prefix, chunk.name, algorithm))
      load(u)
      close(u)
      a.df <- error %.%
        filter(sample.id %in% rownames(err.mat)) %.%
        mutate(param.num=as.numeric(as.character(param.name))) %.%
        group_by(sample.id, param.num) %.%
        summarise(errors=sum(fp+fn),
                  fp=sum(fp),
                  tp=sum(tp),
                  possible.fp=sum(possible.fp),
                  possible.tp=sum(possible.tp),
                  regions=n())
      alist <- split(a.df, a.df$sample.id, drop=TRUE)
      get.mat <- function(var.name){
        alg.mat <- t(sapply(alist, "[[", var.name))
        colnames(alg.mat) <- alist[[1]]$param.num
        alg.mat
      }
      err.list[[algorithm]] <- get.mat("errors")
      fp.list[[algorithm]] <- get.mat("fp")
      tp.list[[algorithm]] <- get.mat("tp")
    }
    dp.peaks.matrices[[set.name]][[chunk.name]] <- err.list
    dp.peaks.matrices.fp[[set.name]][[chunk.name]] <- fp.list
    dp.peaks.matrices.tp[[set.name]][[chunk.name]] <- tp.list
  }
}

save(dp.peaks.matrices,
     dp.peaks.matrices.fp, dp.peaks.matrices.tp,
     file="dp.peaks.matrices.RData")
