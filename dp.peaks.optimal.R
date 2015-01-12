works_with_R("3.1.1",
             "tdhock/PeakSegDP@5bcee97f494dcbc01a69e0fe178863564e9985bc")

load("dp.peaks.matrices.RData")

dp.peaks.optimal <- list()
for(set.name in names(dp.peaks.matrices)){
  chunk.list <- dp.peaks.matrices[[set.name]]
  for(chunk.name in names(chunk.list)){
    err.mat <- chunk.list[[chunk.name]]$PeakSeg
    model.file <- sprintf("data/%s/dp.model.RData", chunk.name)
    load(model.file)
    counts.file <- sprintf("data/%s/counts.RData", chunk.name)
    ## load(counts.file)
    ## counts.list <- split(counts, counts$sample.id)
    for(sample.id in rownames(err.mat)){
      cat(sprintf("%s %s\n", sample.id, chunk.name))
      sample.model <- dp.model[[sample.id]]
      ##sample.counts <- counts.list[[sample.id]]
      ## ggplot()+
      ##   theme_bw()+
      ##   theme(panel.margin=grid::unit(0, "cm"))+
      ##   facet_grid(peaks ~ ., scales="free", labeller=function(var, val){
      ##     s <- ifelse(val==1, "", "s")
      ##     paste0(val, " peak", s)
      ##   })+
      ##   geom_step(aes(chromStart/1e3, coverage),
      ##             data=sample.counts)+
      ##   geom_segment(aes((chromStart-1/2)/1e3, mean,
      ##                    xend=(chromEnd+1/2)/1e3, yend=mean),
      ##                data=sample.model$segments,
      ##                color="green", alpha=3/4, size=1)+
      ##   geom_vline(aes(xintercept=(chromEnd+1/2)/1e3),
      ##              data=sample.model$breaks,
      ##              color="green")
      
      ## With an unconstrained model the loss should always decrease
      ## with model complexity. However with out constrained model
      ## sometimes it increases, for example "H3K36me3_AM_immune/20"
      ## "McGill0004" in this case we consider only the less complex
      ## models.
      loss.df <- sample.model$error
      while(length(next.bigger <- which(diff(loss.df$error) > 0))){
        loss.df <- loss.df[-(next.bigger[1]+1), ]
      }
      dp.peaks.optimal[[set.name]][[chunk.name]][[sample.id]] <-
        with(loss.df, exactModelSelection(error, peaks))
    }
  }
}

save(dp.peaks.optimal, file="dp.peaks.optimal.RData")
