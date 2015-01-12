works_with_R("3.1.1", 
             dplyr="0.3.0.2",
             "tdhock/PeakSegDP@5bcee97f494dcbc01a69e0fe178863564e9985bc",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732")

load("dp.peaks.error.RData")
load("dp.peaks.RData")

### Define what chunk and peaks to download.
set.name <- "H3K4me3_TDH_immune"
chunk.id <- 5
algorithm <- "macs.trained"
prefix <- "http://cbio.ensmp.fr/~thocking/chip-seq-chunk-db"
### Download the data.
## load(url(sprintf("%s/%s/%s/error/%s.RData", prefix, set.name, chunk.id, algorithm)))
## load(url(sprintf("%s/%s/%s/peaks/%s.RData", prefix, set.name, chunk.id, algorithm)))
## load(url(sprintf("%s/%s/%s/regions.RData", prefix, set.name, chunk.id)))

## load(url(sprintf("%s/chunks.RData", prefix)))
chunk.name <- paste0(set.name, "/", chunk.id)
regions <- do.call(rbind, dp.peaks.error[[chunk.name]])
counts.file <- file.path("data", set.name, chunk.id, "counts.RData")
load(counts.file)
model.file <- file.path("data", set.name, chunk.id, "dp.model.RData")
load(model.file)
region.list <- split(regions, regions$sample.id)
sample.ids <- c("McGill0322", "McGill0091", "McGill0002", "McGill0004")

chunk.peaks <- dp.peaks[[chunk.name]]
chunk.peak.list <- list()
chunk.region.list <- list()
chunk.loss.list <- list()
for(sample.id in sample.ids){
  model.list <- chunk.peaks[[sample.id]]
  model.info <- dp.model[[sample.id]]
  chunk.loss.list[[sample.id]] <- data.frame(sample.id, model.info$error)
  sample.regions <- region.list[[sample.id]]
  peak.region.list <- split(sample.regions, sample.regions$param.name)
  for(n.peaks in names(model.list)){
    model.name <- paste0("PeakSeg", n.peaks)
    peak.df <- model.list[[n.peaks]]
    region.df <- peak.region.list[[n.peaks]]
    peak.df <- if(nrow(peak.df)==0){
      data.frame(sample.id=character(), peak.df)
    }else{
      data.frame(sample.id, peak.df)
    }
    chunk.peak.list[[model.name]] <-
      rbind(chunk.peak.list[[model.name]], peak.df)
    chunk.region.list[[model.name]] <-
      rbind(chunk.region.list[[model.name]], region.df)
  }
}
only <- function(x)subset(x, sample.id %in% sample.ids)

PeakSeg4samples <-
  list(peaks=chunk.peak.list,
       regions=chunk.region.list,
       signal=only(counts),
       loss=do.call(rbind, chunk.loss.list))

save(PeakSeg4samples, file="PeakSeg4samples.RData")

