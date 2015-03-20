works_with_R("3.1.3", reshape2="1.2.2", ggplot2="1.0",
             "tdhock/PeakError@d9196abd9ba51ad1b8f165d49870039593b94732",
             "tdhock/PeakSegDP@dd9d52d1341d2824f9bb18d2c9fa84f6e12e6317",
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
  set.chunks <- chunks %>%
    filter(set.name == set)
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
    print(chromPlot)
    ## TODO: why are there no peaks/errors for all data sets except
    ## H3K4me3_TDH_immune ?
    
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
      tit <- paste0(split.id, "/", chunk.name)
              stop(1)

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
        if(clusterStart < prev.chromEnd){
          clusterStart <- prev.chromEnd
        }
        clusterEnd <- peakEnd + expand.bases
        if(next.chromStart < clusterEnd){
          clusterEnd <- next.chromStart
        }
        sample.ids <- unique(paste(cluster.peaks$sample.id))
        cluster.counts <-
          subset(counts,
                 sample.id %in% sample.ids &
                 clusterStart < chromStart &
                 chromEnd < clusterEnd)
        cluster.bases <- clusterEnd-clusterStart
        cluster.counts$count <- cluster.counts$coverage
        ##peak <- multiSampleSegOptimal(cluster.counts)
        peak <- multiSampleSegHeuristic(cluster.counts, 2)
        joint.peak.list[[cluster.id]] <-
          data.frame(sample.id="joint", cluster.id, peak)
        joint.sample.peak.list[[cluster.id]] <-
          data.frame(sample.id=sample.ids, cluster.id, peak)

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
          facet_grid(sample.id ~ ., scales="free_y", labeller=function(var, val){
            sub("McGill0", "", val)
          })+
          geom_segment(aes(chromStart/1e3, 0,
                           xend=chromEnd/1e3, yend=0),
                       data=cluster.peaks,
                       size=2)+
          geom_segment(aes(chromStart/1e3, 0,
                           xend=chromEnd/1e3, yend=0),
                       data=peak,
                       color="green",
                       size=1)+
          geom_step(aes(chromStart/1e3, coverage),
                    data=cluster.counts, color="grey50")
      }
      joint.sample.peaks <- do.call(rbind, joint.sample.peak.list)
      joint.peaks <- do.call(rbind, joint.peak.list)
      errors.by.sample <- split(chunk.errors, chunk.errors$sample.id)
      joint.peaks.by.sample <-
        split(joint.sample.peaks, joint.sample.peaks$sample.id)
        stop(1)
      for(sample.id in names(errors.by.sample)){
        sample.errors <- errors.by.sample[[sample.id]]
        sample.peaks <- joint.peaks.by.sample[[sample.id]]
        if(is.null(sample.peaks)){
          sample.peaks <- Peaks()
        }
        new.errors <- PeakErrorChrom(sample.peaks, sample.errors)
        stop(1)
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
          facet_grid(sample.id ~ ., scales="free_y", labeller=function(var, val){
            sub("McGill0", "", val)
          })

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
