works_with_R("3.1.1",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             "tdhock/PeakSegDP@5bcee97f494dcbc01a69e0fe178863564e9985bc",
             data.table="1.9.4",
             reshape2="1.2.2",
             dplyr="0.3.0.2")

load("dp.peaks.train.RData")
load("dp.peaks.error.RData")
load("dp.peaks.RData")

groups <- dp.peaks.train %>%
  mutate(group=paste(chunk.name, cell.type))
good.group.df <- groups %>%
  group_by(group) %>%
  summarise(region.values=length(unique(regions))) %>%
  filter(region.values == 1)
good.groups <- good.group.df$group
min.error <- groups %>%
  filter(group %in% good.groups) %>%
  group_by(chunk.name, cell.type, algorithm, regions) %>%
  summarise(min=min(errors)) %>%
  mutate(set.name=sub("/.*", "", chunk.name),
         experiment=sub("_.*", "", set.name))
wide <-
  dcast(min.error,
        chunk.name + cell.type + experiment ~ algorithm,
        value.var="min") %>%
  mutate(baseline=ifelse(experiment=="H3K36me3",
           hmcan.broad.trained, macs.trained),
         advantage=baseline-PeakSeg) %>%
  arrange(advantage)
zero.error <- data.table(wide) %>%
  filter(PeakSeg==0)
biggest <- zero.error %>%
  group_by(experiment) %>%
  mutate(rank=rank(-advantage)) %>%
  filter(rank==1)
data.frame(biggest)

## We will make a plot for the window for which we have the biggest
## advantage, for each mark type.
prefix <- "http://cbio.ensmp.fr/~thocking/chip-seq-chunk-db"
chunks.file <- paste0(prefix, "/chunks.RData")
u <- url(chunks.file)
load(u)
close(u)
rownames(chunks) <- with(chunks, paste0(set.name, "/", chunk.id))

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

for(experiment.i in 1:nrow(biggest)){
  chunk.info <- biggest[experiment.i, ]
  experiment <- as.character(chunk.info$experiment)
  chunk.name <- as.character(chunk.info$chunk.name)
  more.info <- chunks[chunk.name, ]
  chunkChrom <- as.character(more.info$chunkChrom)
  cell.type <- as.character(chunk.info$cell.type)
  other.algo <-
    ifelse(experiment=="H3K4me3", "macs.trained", "hmcan.broad.trained")
  algorithms <- c("PeakSeg", other.algo)
  param.err <- dp.peaks.train %>%
    inner_join(chunk.info) %>%
    mutate(param.name=as.character(param.num))
  min.params <- param.err %>%
    filter(algorithm %in% algorithms) %>%
    group_by(algorithm) %>%
    filter(seq_along(errors) == which.min(errors))
  default.param.val <-
    ifelse(other.algo=="macs.trained", "1.30103", "2.30258509299405")
  default.param <- param.err %>%
    filter(algorithm==other.algo,
           param.name==default.param.val) %>%
    mutate(algorithm=sub("trained", "default", algorithm))
  show.params <- rbind(default.param, min.params)
  rownames(show.params) <- show.params$algorithm
  ## TODO: download peaks and error regions for baseline, plot them
  ## alongside PeakSeg model.
  dp.error <- dp.peaks.error[[chunk.name]][[cell.type]]
  dp.param <- show.params["PeakSeg", "param.name"]
  dp.regions <- subset(dp.error, param.name==dp.param)
  sample.ids <- as.character(unique(dp.error$sample.id))
  ## Try to show only a subset of samples.
  sample.ids <- sprintf("McGill%04d", c(4, 5, 24, 26, 29, 107))
  dp.peaks.samples <- dp.peaks[[chunk.name]]
  dp.peak.list <- list()
  for(sample.id in sample.ids){
    dp.peak.list[[sample.id]] <-
      data.frame(sample.id, dp.peaks.samples[[sample.id]][[dp.param]]) %>%
        select(sample.id, chromStart, chromEnd)
  }
  ## Download count signal data.
  counts.file <- file.path("data", chunk.name, "counts.RData")
  load(counts.file)
  sample.counts <- counts %>%
    filter(sample.id %in% sample.ids)
  tit <- with(chunk.info, paste({
    ifelse(experiment=="H3K4me3",
           "sharp peak pattern,", "broad peak pattern,")
  }, cell.type, chunk.name))
  sample.max.df <- sample.counts %>%
    group_by(sample.id) %>%
    summarise(count=max(coverage))
  sample.max <- sample.max.df$count
  names(sample.max) <- as.character(sample.max.df$sample.id)

  other.params <- subset(show.params, algorithm != "PeakSeg")
  trained.param <- subset(other.params, grepl("trained", algorithm))
  trained.algo <- as.character(trained.param$algorithm)

  u <- url(sprintf("%s/%s/error/%s.RData", prefix, chunk.name, trained.algo))
  load(u)
  close(u)
  u <- url(sprintf("%s/%s/peaks/%s.RData", prefix, chunk.name, trained.algo))
  load(u)
  close(u)

  show.peak.list <- list(PeakSeg=do.call(rbind, dp.peak.list))
  show.region.list <- list(PeakSeg=dp.regions)
  for(param.i in 2:1){
    other.param <- other.params[param.i, ]
    other.param.name <- other.param$param.name
    algorithm <- as.character(other.param$algorithm)
    show.region.list[[algorithm]] <- error %>%
      filter(sample.id %in% sample.ids,
             param.name == other.param.name)
    show.peak.list[[algorithm]] <- peaks[[other.param.name]] %>%
      filter(sample.id %in% sample.ids)
  }

  param.desc <-
    c(PeakSeg="maxPeaks",
      macs="log(qvalue)",
      hmcan="log(finalThreshold)")
  compare.region.list <- list()
  compare.peak.list <- list()
  label.sample <- "McGill0107"
  compare.label.list <- list()
  for(algorithm.i in seq_along(show.peak.list)){
    peak.df <- show.peak.list[[algorithm.i]]
    sample.id <- as.character(peak.df$sample.id)
    max.count <- sample.max[sample.id]
    algorithm <- names(show.peak.list)[[algorithm.i]]
    short.algo <- sub("[.].*", "", algorithm)
    this.desc <- param.desc[[short.algo]]
    ## if(algorithm.i==1){
    ##   algorithm.i <- 1
    ## }else{
    ##   algorithm.i <- algorithm.i +1
    ## }
    y.mid <- -algorithm.i*4*max.count/10
    compare.peak.list[[algorithm]] <-
      data.frame(algorithm, y.mid, peak.df)
    ## Also make regions.
    height <- 1
    region.df <- show.region.list[[algorithm]]
    label.i <- which(peak.df$sample.id==label.sample)[1]
    first <- peak.df[label.i, ]
    this.param <- show.params[algorithm, ]
    compare.label.list[[algorithm]] <- with(region.df, {
      data.frame(first, fp=sum(fp), fn=sum(fn),
                 algorithm,
                 y.mid=y.mid[label.i],
                 param.desc=this.desc,
                 param.name=this.param$param.name)
    })
    sample.id <- as.character(region.df$sample.id)
    max.count <- sample.max[sample.id]
    y.min <- (-algorithm.i*4-height)*max.count/10
    y.max <- (-algorithm.i*4+height)*max.count/10
    compare.region.list[[algorithm]] <-
      data.frame(algorithm, y.min, y.max, region.df) %>%
        select(sample.id, y.min, y.max,
               chromStart, chromEnd, annotation, status)
  }
  compare.regions <- do.call(rbind, compare.region.list)
  compare.peaks <- do.call(rbind, compare.peak.list)
  first <- dp.regions %>%
    filter(sample.id==label.sample,
           annotation=="peakStart")
  compare.labels <- do.call(rbind, compare.label.list) %>%
    mutate(chromStart=first$chromStart[2],
           sample.id=label.sample)

  algo.colors <-
    c(macs.default="#A6CEE3", macs.trained="#1F78B4", #lite dark blue
      hmcan.broad.default="#A6CEE3", hmcan.broad.trained="#1F78B4", #lite dark blue
      "#B2DF8A", "#33A02C", #green
      "#FB9A99", "#E31A1C", #red
      "#FDBF6F", "#FF7F00", #orange
      "#CAB2D6", PeakSeg="#6A3D9A", #purple
      "#FFFF99", "#B15928") #yellow/brown
  

  selectedPlot <- 
  ggplot()+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    fill=annotation),
                data=subset(dp.regions, sample.id %in% sample.ids),
                color="grey",
                alpha=1/2)+
  geom_step(aes(chromStart/1e3, coverage),
            data=sample.counts, color="grey50")+
  geom_text(aes(chromStart/1e3, y.mid,
                label=sprintf("%s, %s=%s, %2d FP, %2d FN ",
                  algorithm, param.desc,
                  substr(param.name, 1, 5),
                  fp, fn)),
            data=compare.labels, 
            ##vjust=0.25, size=2,
            size=2.5,
            hjust=1)+
  geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                ymin=y.min, ymax=y.max,
                linetype=status),
            data=subset(compare.regions, sample.id %in% sample.ids),
            fill=NA, color="black", size=0.5)+
  scale_linetype_manual("error type",
                        values=c(correct=0,
                          "false negative"=3,
                          "false positive"=1))+
  geom_point(aes(chromStart/1e3, y.mid, color=algorithm),
             data=compare.peaks,
             pch=1, size=2)+
  geom_segment(aes(chromStart/1e3, y.mid,
                   xend=chromEnd/1e3, yend=y.mid,
                   color=algorithm),
               data=compare.peaks, size=1)+
  scale_color_manual(values=algo.colors)+
  theme_bw()+
  coord_cartesian(xlim=with(more.info, c(expandStart, expandEnd)/1e3))+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
    sub("McGill0", "", val)
    paste(val)
  })+
  scale_y_continuous("aligned read coverage",
                     labels=function(x){
                       sprintf("%.1f", x)
                     },
                     breaks=function(limits){
                       c(0, limits[2])
                     })+
  xlab(paste("position on", chunkChrom, "(kilo base pairs)"))+
  scale_fill_manual("annotation", values=ann.colors,
                    breaks=names(ann.colors))+
  ggtitle(tit)

  png.file <- sprintf("figure-dp-peaks-train-%d.png", experiment.i)
  png(png.file,
      units="in", res=200, width=8, height=length(sample.ids))
  print(selectedPlot)
  dev.off()
  ##system(paste("firefox", png.file))
}
