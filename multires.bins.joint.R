works_with_R("3.1.3", reshape2="1.2.2", ggplot2="1.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegDP@8601b78b937c609b2dcddef79a63ed1899fc0c2b",
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

all.peak.list <- list()
all.error.list <- list()
for(testSet in hyper.dt$testSet){
  print(testSet)
  bases.per.bin <- hyper.dt[testSet]$bases.per.bin
  errors.by.testSet <-
    split(multires.bins$test.error, multires.bins$test.error$testSet)
  set.errors <- errors.by.testSet[[testSet]]
  errors.by.chunk <- split(set.errors, set.errors$chunk.name)
  peaks.by.chrom <- multires.bins$test.peaks[[testSet]]
  setName <- sub(" .*", "", testSet)
  set.i <- sub(".* ", "", testSet)
  set.chunks <- chunks %>%
    filter(set.name == setName)
  chunks.by.chrom <- split(set.chunks, set.chunks$chunkChrom, drop=TRUE)

  chroms <- intersect(names(chunks.by.chrom), names(peaks.by.chrom))
  for(chrom in chroms){

    chrom.chunks <- data.table(chunks.by.chrom[[chrom]]) %>%
      mutate(chunk.name=paste0(set.name, "/", chunk.id),
             split=ifelse(chunk.name %in% names(errors.by.chunk),
               "test", "train"))
    setkey(chrom.chunks, expandStart, expandEnd)
    peaks.by.sample <- peaks.by.chrom[[chrom]]
    chrom.peak.list <- list()
    for(sample.id in names(peaks.by.sample)){
      sample.peaks <- peaks.by.sample[[sample.id]]
      if(nrow(sample.peaks)){ # possible to have no peaks.
        chrom.peak.list[[sample.id]] <-
          data.table(sample.id, sample.peaks)
      }
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

    chromPlot <- 
    ggplot()+
      geom_tallrect(aes(xmin=expandStart/1e3,
                        xmax=expandEnd/1e3,
                        linetype=split),
                    data=chrom.chunks,
                    fill="grey",
                    color="black",
                    alpha=0.5)+
      scale_linetype_manual(values=c(test=1, train=0))+
      geom_segment(aes(chromStart/1e3, sample.id,
                       xend=chromEnd/1e3, yend=sample.id),
                   data=clustered.peaks)
    ##print(chromPlot)
    
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

      cluster.list <- split(chunk.peaks, chunk.peaks$cluster, drop=TRUE)
      joint.peak.list <- list()
      joint.sample.peak.list <- list()
      tit <- paste0(set.i, "/", chunk.name)

      for(cluster.id in names(cluster.list)){
        cluster.peaks <- cluster.list[[cluster.id]]
        peakStart <- min(cluster.peaks$chromStart)
        peakEnd <- max(cluster.peaks$chromEnd)
        expand.bases <- peakEnd-peakStart
        next.id <- paste(as.integer(cluster.id) + 1)
        prev.id <- paste(as.integer(cluster.id) - 1)
        next.chromStart <- min(cluster.list[[next.id]]$chromStart)
        prev.chromEnd <- max(cluster.list[[prev.id]]$chromEnd)
        clusterStart <- peakStart - expand.bases
        clusterEnd <- peakEnd + expand.bases
        fillStart <- min(peakStart-1L, clusterStart)
        fillEnd <- max(peakEnd+1L, clusterEnd)
        if(clusterStart < prev.chromEnd){
          clusterStart <- prev.chromEnd
        }
        if(next.chromStart < clusterEnd){
          clusterEnd <- next.chromStart
        }
        sample.ids <- unique(paste(cluster.peaks$sample.id))
        cluster.counts.unfilled <-
          subset(counts,
                 sample.id %in% sample.ids &
                 clusterStart < chromStart &
                 chromEnd < clusterEnd)
        if(nrow(cluster.counts.unfilled)){
          cluster.counts.unfilled$count <- cluster.counts.unfilled$coverage
          unfilled.list <-
            split(cluster.counts.unfilled,
                  cluster.counts.unfilled$sample.id,
                  drop=TRUE)
          filled.list <- list()
          filled.cols <- c("sample.id", "chromStart", "chromEnd", "count")
          for(sample.id in names(unfilled.list)){
            unfilled <- unfilled.list[[sample.id]]
            filled.list[[sample.id]] <-
              rbind(data.frame(sample.id,
                               chromStart=fillStart,
                               chromEnd=unfilled$chromStart[1],
                               count=unfilled$count[1]),
                    unfilled[, filled.cols],
                    data.frame(sample.id,
                               chromStart=unfilled$chromEnd[nrow(unfilled)],
                               chromEnd=fillEnd,
                               count=unfilled$count[nrow(unfilled)]))
          }
          cluster.counts <- do.call(rbind, filled.list)
          cluster.bases <- clusterEnd-clusterStart
          peak <- multiSampleSegHeuristic(cluster.counts, 2)
          min.peakStart <- as.integer((prev.chromEnd+peakStart)/2+1)
          if(is.finite(min.peakStart) && peak$chromStart < min.peakStart){
            peak$chromStart <- min.peakStart
          }
          max.peakEnd <- as.integer((next.chromStart+peakEnd)/2-1)
          if(is.finite(max.peakEnd) && max.peakEnd < peak$chromEnd){
            peak$chromEnd <- max.peakEnd
          }
          if(with(peak, chromStart < chromEnd)){
            joint.peak.list[[cluster.id]] <-
              data.frame(sample.id="joint", cluster.id, peak)
            joint.sample.peak.list[[cluster.id]] <-
              data.frame(sample.id=sample.ids, cluster.id, peak)
          }
          ggplot()+
            ggtitle(paste(tit, "cluster", cluster.id))+
            xlab(paste("position on", cluster.peaks$chunkChrom[1],
                       "(kilobases = kb)"))+
            theme_bw()+
            scale_y_continuous("aligned read coverage signal",
                               labels=function(x){
                                 sprintf("%.1f", x)
                               },
                               breaks=function(limits){
                                 limits[2]
                               })+
            theme(panel.margin=grid::unit(0, "cm"))+
            facet_grid(sample.id ~ ., labeller=function(var, val){
              sub("McGill0", "", val)
            }, scales="free_y")+
            geom_segment(aes(chromStart/1e3, 0,
                             xend=chromEnd/1e3, yend=0),
                         data=cluster.peaks,
                         size=2)+
            coord_cartesian(xlim=c(fillStart, fillEnd)/1e3)+
            geom_segment(aes(chromStart/1e3, 0,
                             xend=chromEnd/1e3, yend=0),
                         data=peak,
                         color="green",
                         size=1)+
            geom_step(aes(chromStart/1e3, count),
                      data=cluster.counts, color="grey50")
        }
      }
      joint.sample.peaks <- do.call(rbind, joint.sample.peak.list)
      joint.peaks <- do.call(rbind, joint.peak.list)
      errors.by.sample <- split(chunk.errors, chunk.errors$sample.id)
      joint.peaks.by.sample <-
        split(joint.sample.peaks, joint.sample.peaks$sample.id)
      for(sample.id in names(errors.by.sample)){
        sample.errors <- errors.by.sample[[sample.id]]
        sample.peaks <- joint.peaks.by.sample[[sample.id]]
        if(is.null(sample.peaks)){
          sample.peaks <- Peaks()
        }
        new.errors <- PeakErrorChrom(sample.peaks, sample.errors)
        all.error.list[[paste(testSet, chunk.name, sample.id)]] <- 
          data.frame(testSet, chunk.name, sample.id, new.errors)
        if(nrow(sample.peaks)){
          all.peak.list[[testSet]][[chunk.name]][[sample.id]] <- 
            data.frame(testSet, chunk.name, sample.peaks)
        }
      }

      without.joint <- 
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
          coord_cartesian(xlim=with(chunk.peaks[1, ], {
            c(expandStart, expandEnd)/1e3
          }))+
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
          facet_grid(sample.id ~ ., labeller=function(var, val){
            sub("McGill0", "", val)
          }, scales="free_y")

      with.joint <-
        without.joint +
          ggtitle(tit)+
          xlab(paste("position on", cluster.peaks$chunkChrom[1],
                     "(kilobases = kb)"))+
          guides(color="none")+
          geom_segment(aes(chromStart/1e3, 0,
                           color=cluster.id,
                           xend=chromEnd/1e3, yend=0),
                       data=joint.peaks,
                       size=10)+
          geom_text(aes((chromStart+chromEnd)/2/1e3, 0, label=cluster.id),
                    data=joint.peaks)
    }#chunk.name
  }#chrom
}#testSet

multires.bins.joint <-
  list(test.error=do.call(rbind, all.error.list),
       test.peaks=all.peak.list)

save(multires.bins.joint, file="multires.bins.joint.RData")
