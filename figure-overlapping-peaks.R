works_with_R("3.1.3", reshape2="1.2.2", ggplot2="1.0",
             "tdhock/PeakSegDP@854678f21f6886ce99400fac8b81900d16051eb9",
             data.table="1.9.4",
             dplyr="0.4.0",
             xtable="1.7.4")

load("multires.bins.RData")
load("../chip-seq-paper/chunks.RData")

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

hyper.dt <- 
multires.bins$hyperparameters %>%
  filter(variable=="percent")
setkey(hyper.dt, testSet)

set <- "H3K4me3_TDH_immune"
split.id <- 1
testSet <- paste(set, "split", split.id)
bases.per.bin <- hyper.dt[testSet]$bases.per.bin
png.list <- list()
errors.by.testSet <-
  split(multires.bins$test.error, multires.bins$test.error$testSet)
set.errors <- errors.by.testSet[[testSet]]
errors.by.chunk <- split(set.errors, set.errors$chunk.name)
peaks.by.chrom <- multires.bins$test.peaks[[testSet]]
set.chunks <- chunks %>%
  filter(set.name == set)
chunks.by.chrom <- split(set.chunks, set.chunks$chunkChrom, drop=TRUE)

##chrom <- "chr11"
for(chrom in names(chunks.by.chrom)){

chrom.chunks <- data.table(chunks.by.chrom[[chrom]]) %>%
  mutate(chunk.name=paste0(set.name, "/", chunk.id))
setkey(chrom.chunks, expandStart, expandEnd)
peaks.by.sample <- peaks.by.chrom[[chrom]]
chrom.peak.list <- list()
for(sample.id in names(peaks.by.sample)){
  chrom.peak.list[[sample.id]] <-
    data.table(sample.id, peaks.by.sample[[sample.id]])
}
chrom.peaks <- do.call(rbind, chrom.peak.list) %>%
  select(sample.id, chromStart, chromEnd)
clustered.peaks <- clusterPeaks(chrom.peaks)
ggplot()+
  geom_segment(aes(chromStart/1e3, sample.id,
                   xend=chromEnd/1e3, yend=sample.id),
               data=clustered.peaks)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_wrap("cluster", scales="free_x")
setkey(clustered.peaks, chromStart, chromEnd)
peaks.and.chunks <- foverlaps(clustered.peaks, chrom.chunks, nomatch=0L)
peak.chunk.list <- split(peaks.and.chunks, peaks.and.chunks$chunk.name)
valid.chunks <- intersect(names(errors.by.chunk), names(peak.chunk.list))
for(chunk.name in valid.chunks){
  chunk.peaks <- peak.chunk.list[[chunk.name]]
  chunk.errors <- errors.by.chunk[[chunk.name]]
  counts.file <- paste0("data/", chunk.name, "/counts.RData")
  load(counts.file)
  sample.ids <- unique(counts$sample.id)

  cluster.list <- split(chunk.peaks, chunk.peaks$cluster)
  for(cluster.id in names(cluster.list)){
    cluster.peaks <- cluster.list[[cluster.id]]
    peakStart <- min(cluster.peaks$chromStart)
    peakEnd <- max(cluster.peaks$chromEnd)
    expand.bases <- 1
    clusterStart <- peakStart - expand.bases
    clusterEnd <- peakEnd + expand.bases
    cluster.bases <- clusterEnd-clusterStart
  }

  P <- 
ggplot()+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    fill=annotation),
                data=chunk.errors,
                color="grey",
                alpha=1/2)+
  geom_step(aes(chromStart/1e3, coverage),
            data=counts, color="grey50")+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    linetype=status),
                data=chunk.errors,
                fill=NA, color="black", size=1)+
  scale_linetype_manual("error type",
                        values=c(correct=0,
                          "false negative"=3,
                          "false positive"=1))+
  geom_segment(aes(chromStart/1e3, 0,
                   color=factor(cluster),
                   xend=chromEnd/1e3, yend=0),
               data=chunk.peaks,
               size=2)+
  theme_bw()+
  coord_cartesian(xlim=with(chunk.peaks[1, ], c(expandStart, expandEnd)/1e3))+
  scale_y_continuous("aligned read coverage signal",
                         labels=function(x){
                           sprintf("%.1f", x)
                         },
                         breaks=function(limits){
                           limits[2]
                         })+
  scale_fill_manual("annotation", values=ann.colors,
                    breaks=names(ann.colors))+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ ., scales="free_y", labeller=function(var, val){
    sub("McGill0", "", val)
  })

  png.file <-
    paste0("overlapping-peaks/",
           split.id, "/", chunk.name,
           ".png")
  png.dir <- dirname(png.file)
  dir.create(png.dir, recursive = TRUE, showWarnings=FALSE)
  png.list[[png.file]] <- png.file
  cat(png.file, "\n")
  if(!file.exists(png.file)){
  png(png.file,
      units="in", res=200, width=8, height=10)
  print(P)
  dev.off()
  }
}
}
include.lines <-
  paste0("\\begin{figure*}
\\includegraphics[width=\\textwidth]{", unlist(png.list), "}
\\end{figure*}")
cat(include.lines, file="figure-overlapping-peaks.tex")
