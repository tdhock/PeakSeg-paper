works_with_R("3.1.1", dplyr="0.3.0.2")

load("dp.peaks.error.RData")

prefix <- "http://cbio.ensmp.fr/~thocking/chip-seq-chunk-db"

dp.peaks.error[["H3K4me3_TDH_immune/5"]]$tcell %>%
  group_by(param.name) %>%
  summarise(errors=sum(fp+fn), regions=n())

dp.peaks.error[["H3K4me3_TDH_immune/5"]]$bcell %>%
  group_by(param.name) %>%
  summarise(errors=sum(fp+fn), regions=n())

dp.peaks.error[["H3K4me3_TDH_immune/5"]]$monocyte %>%
  group_by(param.name) %>%
  summarise(errors=sum(fp+fn), regions=n())

dp.peaks.train <- NULL
for(group.i in seq_along(dp.peaks.error)){
  chunk.name <- names(dp.peaks.error)[[group.i]]
  cat(sprintf("%4d / %4d %s\n", group.i, length(dp.peaks.error), chunk.name))
  type.list <- dp.peaks.error[[group.i]]
  prev.error <- NULL
  for(algorithm in c("macs.trained", "hmcan.broad.trained")){
    u <- url(sprintf("%s/%s/error/%s.RData", prefix, chunk.name, algorithm))
    load(u)
    close(u)
    prev.error <- rbind(prev.error, data.frame(algorithm, error))
  }

  for(cell.type in names(type.list)){
    dp.error <- type.list[[cell.type]]
    dp.error %>%
      group_by(sample.id, param.name) %>%
      summarise(regions=n())
    sample.ids <- as.character(unique(dp.error$sample.id))
    all.error <-
      rbind(prev.error,
            data.frame(algorithm="PeakSeg", dp.error))
    param.err <- all.error %>%
      filter(sample.id %in% sample.ids) %>%
      mutate(param.num=as.numeric(as.character(param.name))) %>%
      group_by(algorithm, param.num) %>%
      summarise(fp=sum(fp),
                fn=sum(fn),
                errors=sum(fp+fn),
                regions=n())
    ## There may be some PeakSeg models which were not feasible. In
    ## that case there are not the same number of regions in each of
    ## these rows.
    weird <- any(param.err$regions != param.err$regions[1])
    if(weird){
      counts <- with(dp.error, table(sample.id, param.name))
      if(!all(counts %in% c(0, max(counts)))){
        print(counts)
        stop("some regions/models/samples missing")
      }
    }
    dp.peaks.train <- rbind(dp.peaks.train, {
      data.frame(chunk.name, cell.type, param.err)
    })
  }
}

save(dp.peaks.train, file="dp.peaks.train.RData")
