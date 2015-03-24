works_with_R("3.1.3",
             "tdhock/PeakSegDP@f36943548670d51932cbaf235690ff7c01eec74f",
             data.table="1.9.4",
             ggplot2="1.0",
             dplyr="0.4")

## The goal of this script is to make a data set that is small enough
## to be distributed along with the PeakSegDP package as an example. I
## chose 3 chunks on chr11 from the database.

set.name <- "H3K4me3_TDH_immune"
set.dir <- file.path("..", "chip-seq-paper", "chunks", set.name)
chunk.ids <- c(1, 5, 6)
type.list <-
  list(bcell=c(91, 322),
       tcell=c(95, 107))
sample.num <- unlist(type.list)
sample.ids <- sprintf("McGill%04d", sample.num)
region.list <- list()
for(chunk.id in chunk.ids){
  chunk.dir <- file.path(set.dir, chunk.id)
  regions.file <- file.path(chunk.dir, "regions.RData")
  load(regions.file)
  region.list[[chunk.id]] <- regions %>%
    filter(sample.id %in% sample.ids) %>%
    mutate(chunk.id=as.integer(chunk.id))
}
regions <- do.call(rbind, region.list)

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

window <- data.frame(chromStart=117100000, chromEnd=118250000)

h <- 0.3
ggplot()+
  geom_segment(aes(chromStart/1e3, 0,
                   xend=chromEnd/1e3, yend=0),
               data=window)+
  scale_fill_manual(values=ann.colors)+
  geom_rect(aes(xmin=chromStart/1e3, ymin=chunk.id-h,
                xmax=chromEnd/1e3, ymax=chunk.id+h,
                fill=annotation),
            color="black",
            data=regions)

samples.by.region <- split(regions, regions$sample.id)
experiment <- sub("_.*", "", set.name)
count.list <- list()
for(sample.id in sample.ids){
  sample.regions <- samples.by.region[[sample.id]]
  bg.file <- sprintf("~/genomecov/%s/%s.bedGraph", experiment, sample.id)
  print(bg.file)
  one <- fread(bg.file)
  setnames(one, c("chrom", "chromStart", "chromEnd", "count"))
  setkey(one, chrom, chromStart, chromEnd)
  chr11 <- one["chr11"]
  is.before <- chr11$chromEnd < window$chromStart
  is.after <- window$chromEnd < chr11$chromStart
  in.region <- !(is.before | is.after)
  some <- chr11[in.region, ]
  onePlot <- 
  ggplot()+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      fill=annotation),
                  color="grey",
                  alpha=0.5,
                  data=sample.regions)+
    scale_fill_manual(values=ann.colors)+
    geom_step(aes(chromStart/1e3, count),
              data=some,
              color="grey50")
  cell.type <- paste(sample.regions$cell.type[1])
  count.list[[sample.id]] <- some
}

for(sample.id in names(count.list)){
  print(sample.id)
  sample.regions <- samples.by.region[[sample.id]]
  some <- count.list[[sample.id]]
  cell.type <- paste(sample.regions$cell.type[1])
  out.bg <-
    sprintf("~/R/PeakSegDP/inst/exampleData/%s/%s.bedGraph",
            cell.type, sample.id)
  write.table(some, out.bg, quote=FALSE, sep="\t",
              row.names=FALSE, col.names=FALSE)
}
