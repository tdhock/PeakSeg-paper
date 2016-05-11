library(data.table)
library(dplyr)

load("dp.peaks.matrices.RData")

model.info <- list()
for(set.name in names(dp.peaks.matrices)){
  chunk.list <- dp.peaks.matrices[[set.name]]
  for(chunk.name in names(chunk.list)){
    chunk.info <- chunk.list[[chunk.name]]
    regions <- chunk.info$regions
    na.mat <- is.na(chunk.info$PeakSeg)
    models <- rowSums(na.mat)
    min.err <- apply(chunk.info$PeakSeg, 1, min, na.rm=TRUE)
    max.peaks <- apply(chunk.info$PeakSeg, 1, paste, collapse=" ")
    model.info[[paste(set.name, chunk.name)]] <- 
      data.table(set.name, chunk.name,
                 sample.id=names(regions), regions, models, min.err, max.peaks)
  }
}
models.tall <- do.call(rbind, model.info)
models.tall %>%
  summarise(regions=sum(regions),
            models=sum(models),
            problems=n(),
            all10=sum(models == 0))
## 14 problems with less than 10 models... there are 3 for which there
## is some error.
models.tall %>%
  filter(models != 0,
         min.err != 0)
dp.peaks.matrices$H3K4me3_XJ_immune[["H3K4me3_XJ_immune/2"]]$PeakSeg["McGill0028", ]
14/2752*100 # 0.5% of problems for which we did not compute all 10 models.
3/2752*100 # 0.1% of problems for which this actually caused an error.
models.tall %>%
  filter(models != 0,
         !grepl("NA$", max.peaks))
5/2752*100 # 0.2% could compute 9 peaks model but not a smaller model.
## TODO: plot these data/models!

load("dp.timings.RData")

dp.timings %>%
  summarise(seconds=sum(seconds)) %>%
  mutate(minutes=seconds/60,
         hours=minutes/60,
         days=hours/24)

dp.peaks.NA <- models.tall %>%
  filter(models != 0) %>%
  data.frame
save(dp.peaks.NA, file="dp.peaks.NA.RData")
