works_with_R("3.1.2",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/ggplot2@aac38b6c48c016c88123208d497d896864e74bd7",
             "tdhock/PeakSegDP@5bcee97f494dcbc01a69e0fe178863564e9985bc",
             "Rdatatable/data.table@200b5b40dd3b05112688c3a9ca2dd41319c2bbae",
             reshape2="1.2.2",
             dplyr="0.4.0")

load("dp.peaks.train.RData")
load("dp.peaks.error.RData")
load("dp.peaks.RData")

chunk.name <- "H3K36me3_AM_immune/8"
chunk.name <- "H3K4me3_PGP_immune/7"

counts.file <- file.path("data", chunk.name, "counts.RData")
load(counts.file)

counts.list <- split(counts, counts$sample.id)
sample.id <- "McGill0026"
sample.counts <- counts.list[[sample.id]]
cell.type <- as.character(sample.counts$cell.type[1])
errors <- dp.peaks.error[[chunk.name]][[cell.type]]
error.list <- split(errors, errors$sample.id)
sample.errors <- error.list[[sample.id]]
sample.error.list <- split(sample.errors, sample.errors$param.name)
regions <- sample.error.list[[1]]

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

sample.counts$weight <- with(sample.counts, chromEnd-chromStart)
n <- length(count)
l <- 800
seg2.starts <- as.integer(seq(1, n, l=l)[-c(1, l)])
loss.list <- list()
for(model.i in seq_along(seg2.starts)){
  seg2.start <- seg2.starts[[model.i]]
  seg2.chromStart <- sample.counts$chromStart[seg2.start]
  seg.starts <- c(1, seg2.start)
  seg.ends <- c(seg2.start-1, n)
  for(seg.i in seq_along(seg.starts)){
    seg.start <- seg.starts[[seg.i]]
    seg.end <- seg.ends[[seg.i]]
    seg.data <- sample.counts[seg.start:seg.end, ]
    seg.mean <- with(seg.data, sum(coverage * weight)/sum(weight))
    seg.loss <- with(seg.data, PoissonLoss(coverage, seg.mean, weight))
    loss.list[[paste(model.i, seg.i)]] <-
      data.table(model.i, seg.i, seg2.chromStart,
                 seg2.start, seg.start, seg.end, seg.mean, seg.loss)
  }
}
loss.dt <- do.call(rbind, loss.list)
model.dt <- loss.dt %>%
  group_by(model.i, seg2.chromStart) %>%
  summarise(loss=sum(seg.loss))

selectedPlot <- 
  ggplot()+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    fill=annotation),
                data=data.table(regions, what="profile"),
                color="grey",
                alpha=1/2)+
  geom_step(aes(chromStart/1e3, coverage),
            data=data.table(sample.counts, what="profile"),
            color="grey50")+
  geom_line(aes(seg2.chromStart/1e3, loss),
            data=data.table(model.dt, what="loss"))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(what ~ ., scales="free")+
  xlab(paste("position on chromosome (kb = kilo bases)"))+
  scale_fill_manual("label", values=ann.colors,
                    breaks=names(ann.colors))

png("figure-dp.png",
    units="in", res=200, width=8, height=2.5)
print(selectedPlot)
dev.off()


